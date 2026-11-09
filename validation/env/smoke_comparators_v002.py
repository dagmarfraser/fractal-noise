"""Multi-seed smoke test (v002): every Python comparator on K shared exact-fGn series (alpha = 0.5).
Usage: python smoke_comparators_v002.py <series.csv (N x K)> <fs> <seeds.csv (K)>.
Prints RES|tool|k|value|status lines. Generators are seeded per k from seeds.csv."""
import sys, warnings, traceback
import numpy as np
warnings.filterwarnings("ignore")
X = np.loadtxt(sys.argv[1], delimiter=","); X = X[:, None] if X.ndim == 1 else X
fs = float(sys.argv[2]); seeds = np.atleast_1d(np.loadtxt(sys.argv[3], dtype=np.int64))
def slope(f, p, lo, hi):
    m = (f >= lo) & (f <= hi) & (p > 0)
    return -np.polyfit(np.log10(f[m]), np.log10(p[m]), 1)[0]
def run(name, fn, k):
    try: print(f"RES|{name}|{k}|{float(fn()):.6g}|ok")
    except Exception as e: print(f"RES|{name}|{k}|nan|FAIL {type(e).__name__}: {str(e)[:150]}")
def _cn():
    import colorednoise as cn; y = cn.powerlaw_psd_gaussian(2.0, 4096, random_state=int(seed))
    from scipy.signal import welch; f, p = welch(y, 1.0, nperseg=1024); return slope(f, p, 0.005, 0.2)
def _ndsp_sim():
    from neurodsp.sim import sim_powerlaw; from scipy.signal import welch
    y = sim_powerlaw(40.96, 100.0, exponent=-2.0); f, p = welch(y, 100.0, nperseg=1024); return slope(f, p, 1, 20)
def _ndsp_irasa():
    from neurodsp.aperiodic import compute_irasa; f, pa, _ = compute_irasa(x, fs, f_range=(1, 20)); return slope(f, pa, 1, 20)
def _yasa():
    import yasa; f, pa, po, fit = yasa.irasa(x, sf=fs, band=(1, 20), return_fit=True); return -float(fit["Slope"].iloc[0])
def _pyrasa():
    from pyrasa.irasa import irasa
    out = irasa(x, fs=int(fs), band=(1, 20), nperseg=1024); ap = out.fit_aperiodic_model(fit_func="fixed")
    return float(ap.aperiodic_params["Exponent"].iloc[0])
def _fooof():
    from fooof import FOOOF; from scipy.signal import welch
    f, p = welch(x, fs, nperseg=1024); fm = FOOOF(verbose=False); fm.fit(f, p, [1, 20]); return fm.aperiodic_params_[-1]
def _nolds():
    import nolds; return 2 * nolds.dfa(x) - 1
def _nolitsa():
    from nolitsa import surrogates; from scipy.signal import welch
    y = surrogates.iaaft(x)[0]; f, p = welch(y, fs, nperseg=1024); return slope(f, p, 1, 20)
def _fbm():
    from fbm import FBM; from scipy.signal import welch
    y = FBM(n=4096, hurst=0.75, method="daviesharte").fgn(); f, p = welch(y, 1.0, nperseg=1024); return slope(f, p, 0.005, 0.2)
for k in range(X.shape[1]):
    x = X[:, k]; seed = seeds[k]; np.random.seed(int(seed) % 2**32)   # neurodsp / nolitsa / fbm use the global RNG
    for n, fn in [("colorednoise_gen_a2", _cn), ("neurodsp_sim_a2", _ndsp_sim), ("neurodsp_irasa", _ndsp_irasa),
                  ("yasa_irasa", _yasa), ("pyrasa_irasa", _pyrasa), ("fooof", _fooof), ("nolds_dfa", _nolds),
                  ("nolitsa_iaaft_slope", _nolitsa), ("fbm_dh_gen_a05", _fbm)]:
        run(n, fn, k + 1)
