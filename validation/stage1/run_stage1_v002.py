"""Stage-1 estimator runner, Python arm (v002). Hub-neutral: reads the manifest + series CSVs, writes long-format estimates.
Usage: python run_stage1_v002.py <manifest.csv> <series_dir> <out.csv> <fs> <f_lo> <f_hi> [n_jobs]
Manifest columns: cell_id,alpha,N,k,H,nrep,seed,series_file[,alpha_eff]. Series file: N x nrep CSV.
v002 vs v001: Welch/IRASA window = min(1024, N/2) for every tool; yasa run twice (package default detrend='constant',
which collapses at steep alpha, and detrend='linear'); DFA orders 1-3; joblib over series; k_used column (NA here)."""
import sys, time, warnings, csv
import numpy as np
from joblib import Parallel, delayed
warnings.filterwarnings("ignore")
man, sdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
fs, flo, fhi = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])
njobs = int(sys.argv[7]) if len(sys.argv) > 7 else 1

def slope(f, p):
    m = (f >= flo) & (f <= fhi) & (p > 0)
    return -np.polyfit(np.log10(f[m]), np.log10(p[m]), 1)[0]
def nps(x): return int(min(1024, len(x) // 2))
def e_ndsp(x):
    from neurodsp.aperiodic import compute_irasa
    f, pa, _ = compute_irasa(x, fs, f_range=(flo, fhi), nperseg=nps(x)); return slope(f, pa)
def yasa_fit(x, **kw):
    import yasa, logging; logging.getLogger("yasa").setLevel(logging.ERROR)
    _, _, _, fit = yasa.irasa(x, sf=fs, band=(flo, fhi), win_sec=nps(x) / fs, return_fit=True, verbose=False, **kw)
    return -float(fit["Slope"].iloc[0])
def e_yasa_def(x): return yasa_fit(x)
def e_yasa_lin(x): return yasa_fit(x, kwargs_welch=dict(average="median", window="hann", detrend="linear"))
def e_pyrasa(x):
    from pyrasa.irasa import irasa
    o = irasa(x, fs=int(fs), band=(flo, fhi), nperseg=nps(x)); return float(o.fit_aperiodic_model(fit_func="fixed").aperiodic_params["Exponent"].iloc[0])
def e_fooof(x):
    from fooof import FOOOF; from scipy.signal import welch
    f, p = welch(x, fs, nperseg=nps(x)); fm = FOOOF(verbose=False); fm.fit(f, p, [flo, fhi]); return fm.aperiodic_params_[-1]
def dfa(order):
    def e(x):
        import nolds; nv = np.unique(np.logspace(1, np.log10(len(x) // 10), 15).astype(int))
        return 2 * nolds.dfa(x, nvals=nv, order=order) - 1
    return e
EST = [("neurodsp_irasa", e_ndsp), ("yasa_irasa_default", e_yasa_def), ("yasa_irasa_lindetrend", e_yasa_lin),
       ("pyrasa_irasa", e_pyrasa), ("fooof_fixed", e_fooof), ("nolds_dfa1", dfa(1)), ("nolds_dfa2", dfa(2)), ("nolds_dfa3", dfa(3))]

def one(cid, r, x):
    rows = []
    for name, fn in EST:
        t0 = time.perf_counter()
        try: a, st = float(fn(x)), "ok"
        except Exception as e: a, st = float("nan"), f"FAIL {type(e).__name__}: {str(e)[:120]}"
        if st == "ok" and not np.isfinite(a): st = "FAIL non-finite"
        rows.append([cid, r, name, "python", f"{a:.6g}", f"{time.perf_counter() - t0:.4g}", st, "NA"])
    return rows

jobs = []
with open(man) as fm_:
    for row in csv.DictReader(fm_):
        X = np.loadtxt(f"{sdir}/{row['series_file']}", delimiter=","); X = X[:, None] if X.ndim == 1 else X
        jobs += [(row["cell_id"], r + 1, X[:, r].copy()) for r in range(X.shape[1])]
res = Parallel(n_jobs=njobs)(delayed(one)(*j) for j in jobs)
with open(out, "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["cell_id", "rep", "estimator", "language", "alpha_hat", "seconds", "status", "k_used"])
    for rows in res: w.writerows(rows)
