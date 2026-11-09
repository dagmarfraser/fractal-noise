"""Stage-1 estimator runner, Python arm (v003). Hub-neutral: reads the manifest + series CSVs, writes long-format estimates.
Usage: python run_stage1_v003.py <manifest.csv> <series_dir> <out.csv> <fs> <band_label> <n_jobs> <all|spectral>
The fitting band is per cell, from manifest columns f_lo_<band_label>, f_hi_<band_label> (Hz).
'spectral' runs only the band-dependent estimators (DFA is band-free), so a second band costs less.
Output columns: cell_id,rep,estimator,language,alpha_hat,seconds,status,k_used,band (band = 'none' for DFA).
Welch/IRASA window = min(1024, N/2) for every tool. yasa runs twice: package default (detrend='constant',
collapses at steep alpha) and detrend='linear'. v003 vs v002: per-cell band, band column, spectral-only mode."""
import sys, time, warnings, csv
import numpy as np
from joblib import Parallel, delayed
warnings.filterwarnings("ignore")
man, sdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
fs, blab, njobs, subset = float(sys.argv[4]), sys.argv[5], int(sys.argv[6]), sys.argv[7]
if subset not in ("all", "spectral"): sys.exit(f"subset must be all|spectral, got {subset}")

def slope(f, p, b):
    m = (f >= b[0]) & (f <= b[1]) & (p > 0)
    return -np.polyfit(np.log10(f[m]), np.log10(p[m]), 1)[0]
def nps(x): return int(min(1024, len(x) // 2))
def e_ndsp(x, b):
    from neurodsp.aperiodic import compute_irasa
    f, pa, _ = compute_irasa(x, fs, f_range=b, nperseg=nps(x)); return slope(f, pa, b)
def yasa_fit(x, b, **kw):
    import yasa, logging; logging.getLogger("yasa").setLevel(logging.ERROR)
    _, _, _, fit = yasa.irasa(x, sf=fs, band=b, win_sec=nps(x) / fs, return_fit=True, verbose=False, **kw)
    return -float(fit["Slope"].iloc[0])
def e_yasa_def(x, b): return yasa_fit(x, b)
def e_yasa_lin(x, b): return yasa_fit(x, b, kwargs_welch=dict(average="median", window="hann", detrend="linear"))
def e_pyrasa(x, b):
    from pyrasa.irasa import irasa
    o = irasa(x, fs=int(fs), band=b, nperseg=nps(x)); return float(o.fit_aperiodic_model(fit_func="fixed").aperiodic_params["Exponent"].iloc[0])
def e_fooof(x, b):
    from fooof import FOOOF; from scipy.signal import welch
    f, p = welch(x, fs, nperseg=nps(x)); fm = FOOOF(verbose=False); fm.fit(f, p, list(b)); return fm.aperiodic_params_[-1]
def dfa(order):
    def e(x, b):
        import nolds; nv = np.unique(np.logspace(1, np.log10(len(x) // 10), 15).astype(int))
        return 2 * nolds.dfa(x, nvals=nv, order=order) - 1
    return e
SPECTRAL = [("neurodsp_irasa", e_ndsp), ("yasa_irasa_default", e_yasa_def), ("yasa_irasa_lindetrend", e_yasa_lin),
            ("pyrasa_irasa", e_pyrasa), ("fooof_fixed", e_fooof)]
BANDFREE = [("nolds_dfa1", dfa(1)), ("nolds_dfa2", dfa(2)), ("nolds_dfa3", dfa(3))]
EST = [(n, f, blab) for n, f in SPECTRAL] + ([(n, f, "none") for n, f in BANDFREE] if subset == "all" else [])

def one(cid, r, x, b):
    rows = []
    for name, fn, lab in EST:
        t0 = time.perf_counter()
        try: a, st = float(fn(x, b)), "ok"
        except Exception as e: a, st = float("nan"), f"FAIL {type(e).__name__}: {str(e)[:120]}"
        if st == "ok" and not np.isfinite(a): st = "FAIL non-finite"
        rows.append([cid, r, name, "python", f"{a:.6g}", f"{time.perf_counter() - t0:.4g}", st, "NA", lab])
    return rows

jobs = []
with open(man) as fm_:
    for row in csv.DictReader(fm_):
        X = np.loadtxt(f"{sdir}/{row['series_file']}", delimiter=","); X = X[:, None] if X.ndim == 1 else X
        b = (float(row[f"f_lo_{blab}"]), float(row[f"f_hi_{blab}"]))
        jobs += [(row["cell_id"], r + 1, X[:, r].copy(), b) for r in range(X.shape[1])]
res = Parallel(n_jobs=njobs)(delayed(one)(*j) for j in jobs)
with open(out, "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["cell_id", "rep", "estimator", "language", "alpha_hat", "seconds", "status", "k_used", "band"])
    for rows in res: w.writerows(rows)
