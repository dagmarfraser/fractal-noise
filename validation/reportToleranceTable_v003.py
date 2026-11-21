# reportToleranceTable_v003.py  Measured variance factor of estimateWhittleAlpha against the cosine-bell constant of Velasco (1999).
# Run:  python3 -I reportToleranceTable_v003.py   (from any folder). Reads out_toleranceTable_v001/tolerance_table.csv (written by reportToleranceTable_v001.py); writes out_toleranceTable_v003/report.txt.
#
# Supersedes the note in v002, which said Velasco's paper was not in the library. Source read, with page images:
# Velasco, C. (1999). Gaussian semiparametric estimation of non-stationary time series. Journal of Time Series Analysis, 20(1), 87-127.
#   Theorem 6 (p. 100): m^(1/2) (dHat - d0) -> N(0, p Phi / 4). Equation (10) (p. 101) defines Phi.
#   p. 101, after equation (10): for the full cosine window h_t = (1 - cos(2 pi t / n)) / 2, with mu = 0 and d < 3/2, Phi = 1 if the taper is treated as order p = 3 (every third frequency),
#   and Phi = 35/18 if all Fourier frequencies from lambda_2 to lambda_m are used. p. 97 says the cosine bell is not treated explicitly, so this is a stated remark, not a proved theorem.
# Hurvich & Chen (2000, p. 159) print this as "35/18m". Equation (10) with p = 1 gives Phi = 35/18, so the variance of dHat is (35/18) / (4m); the printed form omits the 4.
#
# What differs from Velasco's setting (the comparison below is empirical, not a consequence of his theorem):
#   the taper is applied after differencing k times and the periodogram is then recoloured; k is chosen from the data; the band starts at 1 Hz, not at lambda_2;
#   the band is fixed in Hz, so m / n does not tend to 0 as condition (9) requires.
#
# Part 3 (derived here, not from a paper): if periodogram ordinates at bins j and l have correlation rho(|j - l|)^2, with rho = 1, 2/3, 1/6 for the cosine bell, the first-order expansion of the estimate gives
#   factor over sdBand^2 = sum_jl nu_j nu_l rho(|j - l|)^2 / sum_j nu_j^2, nu = log(lambda_j) - mean(log(lambda_j)). Over a band starting at the lowest frequency with large m this tends to 35/18.
import math
import os
import sys

import numpy as np
import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
out = os.path.join(here, 'out_toleranceTable_v003')
os.makedirs(out, exist_ok=True)
T = pd.read_csv(os.path.join(here, 'out_toleranceTable_v001', 'tolerance_table.csv'))
PHI = 35 / 18
lines = []


def phi10(n, p):
    """Equation (10) of Velasco (1999) at finite n, for the taper estimateWhittleAlpha uses."""
    h2 = (0.5 * (1 - np.cos(2 * np.pi * np.arange(n) / n))) ** 2
    C = np.real(np.fft.fft(h2))
    return float((C[np.arange(0, n - p + 1, p)] ** 2).sum() / h2.sum() ** 2)


def bandFactor(fs, N, sdBand):
    """Finite-band factor (part 3). k is recovered by matching the sdBand stored by v001; no match is an error."""
    for k in range(4):
        M = N - k
        f = np.arange(M) * fs / M
        f = f[(f >= 1.0) & (f <= 20.0)]
        nu = np.log(f) - np.log(f).mean()
        if abs(1 / math.sqrt((nu ** 2).sum()) - sdBand) < 1e-12:
            rho2 = {0: 1.0, 1: (2 / 3) ** 2, 2: (1 / 6) ** 2}
            return sum(r * (2 if lag else 1) * float((nu[:len(nu) - lag] * nu[lag:]).sum()) for lag, r in rho2.items()) / float((nu ** 2).sum())
    raise RuntimeError(f'No differencing order reproduces sdBand for fs {fs}, N {N}.')


lines.append('1. Equation (10) evaluated for the estimateWhittleAlpha taper')
for n in (1000, 3000):
    lines.append(f'   n = {n}: p = 1 gives {phi10(n, 1):.6f} (35/18 = {PHI:.6f}); p = 3 gives {phi10(n - n % 3, 3):.6f}')

T['measured'] = (T.sdObserved / T.sdBand) ** 2
T['bandFactor'] = [bandFactor(r.fs, r.N, r.sdBand) for r in T.itertuples()]
T['approxSE'] = T.measured * 1.2533 * math.sqrt(2 / 199) / np.sqrt(T.cells)  # variance of 200 series per alpha, median over the alpha cells in the row; approximate
z = (T.measured - PHI) / T.approxSE
lines.append(f'\n2. Measured variance factor over sdBand^2, exact references, {len(T)} rows (fs x N x alpha range)')
lines.append(f'   median {T.measured.median():.3f}, range {T.measured.min():.3f} to {T.measured.max():.3f}; 35/18 = {PHI:.3f}; median measured / (35/18) = {(T.measured / PHI).median():.3f}')
lines.append(f'   approximate standard error per row {T.approxSE.median():.3f}; rows within 2 approximate standard errors of 35/18: {int((z.abs() <= 2).sum())} of {len(T)}')
lines.append('\n3. Finite-band factor (derived here) for the same rows')
lines.append(f'   median {T.bandFactor.median():.3f}, range {T.bandFactor.min():.3f} to {T.bandFactor.max():.3f}; median measured / finite-band factor = {(T.measured / T.bandFactor).median():.3f}')
text = '\n'.join(lines) + '\n'
T[['fs', 'N', 'alphaRange', 'cells', 'mBins', 'sdBand', 'sdObserved', 'measured', 'approxSE', 'bandFactor']].to_csv(os.path.join(out, 'factor_table.csv'), index=False)
with open(os.path.join(out, 'report.txt'), 'w') as f:
    f.write(text)
print(text)
