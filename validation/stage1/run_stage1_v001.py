"""Stage-1 estimator runner, Python arm (v001). Hub-neutral: reads the manifest + series CSVs, writes long-format estimates.
Usage: python run_stage1_v001.py <manifest.csv> <series_dir> <out.csv> <fs> <f_lo> <f_hi>
Manifest columns: cell_id,alpha,N,k,H,nrep,seed,series_file. Series file: N x nrep CSV."""
import sys, time, warnings, csv
import numpy as np
warnings.filterwarnings("ignore")
man, sdir, out = sys.argv[1], sys.argv[2], sys.argv[3]
fs, flo, fhi = float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])

def slope(f, p):
    m = (f >= flo) & (f <= fhi) & (p > 0)
    return -np.polyfit(np.log10(f[m]), np.log10(p[m]), 1)[0]
def e_ndsp(x):
    from neurodsp.aperiodic import compute_irasa
    f, pa, _ = compute_irasa(x, fs, f_range=(flo, fhi)); return slope(f, pa)
def e_yasa(x):
    import yasa; _, _, _, fit = yasa.irasa(x, sf=fs, band=(flo, fhi), return_fit=True); return -float(fit["Slope"].iloc[0])
def e_pyrasa(x):
    from pyrasa.irasa import irasa
    o = irasa(x, fs=int(fs), band=(flo, fhi), nperseg=1024); return float(o.fit_aperiodic_model(fit_func="fixed").aperiodic_params["Exponent"].iloc[0])
def e_fooof(x):
    from fooof import FOOOF; from scipy.signal import welch
    f, p = welch(x, fs, nperseg=1024); fm = FOOOF(verbose=False); fm.fit(f, p, [flo, fhi]); return fm.aperiodic_params_[-1]
def e_dfa(x):
    import nolds; return 2 * nolds.dfa(x, nvals=np.unique(np.logspace(1, np.log10(400), 15).astype(int))) - 1
EST = [("neurodsp_irasa", e_ndsp), ("yasa_irasa", e_yasa), ("pyrasa_irasa", e_pyrasa), ("fooof_fixed", e_fooof), ("nolds_dfa1", e_dfa)]

with open(man) as fm_, open(out, "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["cell_id", "rep", "estimator", "language", "alpha_hat", "seconds", "status"])
    for row in csv.DictReader(fm_):
        X = np.loadtxt(f"{sdir}/{row['series_file']}", delimiter=","); X = X[:, None] if X.ndim == 1 else X
        for r in range(X.shape[1]):
            for name, fn in EST:
                t0 = time.perf_counter()
                try: a, st = float(fn(X[:, r])), "ok"
                except Exception as e: a, st = float("nan"), f"FAIL {type(e).__name__}: {str(e)[:120]}"
                w.writerow([row["cell_id"], r + 1, name, "python", f"{a:.6g}", f"{time.perf_counter() - t0:.4g}", st])
