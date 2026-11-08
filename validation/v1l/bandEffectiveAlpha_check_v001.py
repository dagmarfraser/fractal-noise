"""Band-effective alpha of exact V0 spectra (OLS of log S on log f over 400 LINEARLY spaced points in the band, fs = 100).
Source of the first in-band numbers quoted 2026-11-08 ~09:50 (alpha 4 -> 3.906, 0.5 -> 0.530, 6.5 -> 6.389 over 1-20 Hz)
and the differencing local slope 2x cot x at 20 Hz (1.73). The stage-1 runs use bandEffectiveAlpha_v001.m instead
(Welch grid, nperseg = min(1024, N/2)), which gives e.g. 3.905 at alpha 4, N 4096. Usage: python bandEffectiveAlpha_check_v001.py"""
import numpy as np
from exact_ref_v001 import exact_ref_series, alpha_to_kH, fgn_spectrum
fs = 100
for alpha in [0.5, 2.5, 4, 4.5, 6.5]:
    k, H = alpha_to_kH(alpha)
    for band in [(1, 20), (1, 10), (1, 5)]:
        f = np.linspace(band[0], band[1], 400); lam = 2 * np.pi * f / fs
        S = fgn_spectrum(lam, H) / (2 - 2 * np.cos(lam)) ** k
        print(f"alpha {alpha} band {band}: effective {-np.polyfit(np.log(f), np.log(S), 1)[0]:.3f}")
x = np.pi * 20 / fs; print("local slope of one differencing step at 20 Hz:", round(2 * x / np.tan(x), 4))
