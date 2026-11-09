"""V1L: estimator leakage floor at steep alpha. Whole-record multitaper slope on exact V0 series, with manipulations:
  M0 raw (NW=4, K=7, eigenspectra averaged)   M0a Thomson adaptive weights (as MATLAB pmtm default 'adapt')        M1 prewhiten: diff k, multitaper, recolour by |1-e^{-iw}|^{-2k}
  M1e prewhiten with k estimated from the data (DFA-3 pre-estimate, k = floor((a+1)/2), as run_stage1_v003.R)
  M2 more tapers (NW=8, K=14)          M3 linear detrend      M4 end-match (subtract line x1->xN)   M5 single Hann taper
Scored against the band-effective alpha of the analytic spectrum on the same frequency grid.
Usage: python v1l_leakage_v001.py <out.csv> [nrep] [seed]"""
import sys, csv, numpy as np
from scipy.signal.windows import dpss, hann
from exact_ref_v001 import exact_ref_series

FS, F_HI = 100.0, 20.0
def mt_psd(x, NW=4, K=7, taper=None):
    N = len(x)
    W = dpss(N, NW, K) if taper is None else taper[None, :] / np.sqrt(np.sum(taper ** 2))
    P = np.mean(np.abs(np.fft.rfft(W * x[None, :], axis=1)) ** 2, 0)
    return np.fft.rfftfreq(N, 1 / FS), P
def mt_psd_adapt(x, NW=4, K=7, iters=50):
    """Thomson (1982) adaptive weighting, as in MATLAB pmtm('adapt'): down-weights tapers whose leakage dominates."""
    N = len(x); W, lam = dpss(N, NW, K, return_ratios=True)
    Sk = np.abs(np.fft.rfft(W * x[None, :], axis=1)) ** 2
    sig2 = np.var(x); S = np.mean(Sk[:2], 0)
    for _ in range(iters):
        b = S[None, :] / (lam[:, None] * S[None, :] + (1 - lam[:, None]) * sig2)
        w = b ** 2 * lam[:, None]
        Sn = np.sum(w * Sk, 0) / np.sum(w, 0)
        if np.max(np.abs(Sn - S) / S) < 1e-6: S = Sn; break
        S = Sn
    return np.fft.rfftfreq(N, 1 / FS), S
def slope(f, P, band):
    m = (f >= band[0]) & (f <= band[1]) & (P > 0)
    return -np.polyfit(np.log10(f[m]), np.log10(P[m]), 1)[0]
def dfa_alpha(x, order=3):
    y = np.cumsum(x - x.mean()); n = len(y)
    s = np.unique(np.round(np.exp(np.linspace(np.log(10), np.log(n // 10), 15))).astype(int)); F = []
    for w in s:
        nb = n // w; Y = y[:nb * w].reshape(nb, w).T; t = np.arange(w)
        V = np.vander((t - t.mean()) / w, order + 1); Q, _ = np.linalg.qr(V)
        F.append(np.sqrt(np.mean((Y - Q @ (Q.T @ Y)) ** 2)))
    return 2 * np.polyfit(np.log(s), np.log(F), 1)[0] - 1
def prewhiten(x, k, band):
    d = np.diff(x, k) if k > 0 else x
    f, P = mt_psd(d); lam = 2 * np.pi * f / FS
    with np.errstate(divide="ignore"): Pr = P / (2 - 2 * np.cos(lam)) ** k
    return slope(f, Pr, band)
def methods(x, k, band):
    N = len(x); t = np.arange(N); out = {}
    out["M0_raw"] = slope(*mt_psd(x), band)
    out["M0a_adaptive"] = slope(*mt_psd_adapt(x), band)
    out["M1_prewhiten_oraclek"] = prewhiten(x, k, band)
    out["M1e_prewhiten_estk"] = prewhiten(x, int(min(max(np.floor((dfa_alpha(x) + 1) / 2), 0), 3)), band)
    out["M2_more_tapers"] = slope(*mt_psd(x, NW=8, K=14), band)
    out["M3_lin_detrend"] = slope(*mt_psd(x - np.polyval(np.polyfit(t, x, 1), t)), band)
    out["M4_end_match"] = slope(*mt_psd(x - (x[0] + (x[-1] - x[0]) * t / (N - 1))), band)
    out["M5_hann"] = slope(*mt_psd(x, taper=hann(N, sym=False)), band)
    return out

if __name__ == "__main__":
    out = sys.argv[1]; nrep = int(sys.argv[2]) if len(sys.argv) > 2 else 30; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261108
    rng = np.random.default_rng(seed)
    with open(out, "w", newline="") as fo:
        w = csv.writer(fo); w.writerow(["alpha", "N", "band", "alpha_eff", "rep", "method", "alpha_hat"])
        for alpha in [0.5, 2.5, 4.0, 4.5, 5.5, 6.5]:
            for N in [256, 1024, 4096]:
                X, info = exact_ref_series(N, alpha, nrep, rng)
                f = np.fft.rfftfreq(N, 1 / FS)
                for bname, band in [("record", (20 * FS / N, F_HI)), ("fixed", (1.0, F_HI))]:
                    m = (f >= band[0]) & (f <= band[1])
                    a_eff = -np.polyfit(np.log10(f[m]), np.log10(info["S"](2 * np.pi * f[m] / FS)), 1)[0]
                    for r in range(nrep):
                        for name, a in methods(X[:, r], info["k"], band).items():
                            w.writerow([alpha, N, bname, f"{a_eff:.4f}", r + 1, name, f"{a:.4f}"])
