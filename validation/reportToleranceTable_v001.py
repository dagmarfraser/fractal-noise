# reportToleranceTable_v001.py  How large is AlphaTol = 0.3 relative to the sampling error of estimateWhittleAlpha, on exact references only?
# Run:  python3 -I reportToleranceTable_v001.py   (from any folder). Writes out_toleranceTable_v001/tolerance_table.csv and report_tables.txt, and prints the same text.
# Inputs (existing runs, no new simulation): out_lowAlpha_v002/cells.csv and out_highAlpha_v003/cells.csv (exact references, 200 series per cell),
# out_drawCalibration_v001/series.csv (per-draw gaps, surrogate alpha minus target).
#
# Predicted SD of the exponent for the band's m Hann periodogram bins (series differenced k = floor((alpha + 1) / 2) times first, as the estimator does on average):
#   sdRobinson = 1 / sqrt(m). Local Whittle: sqrt(m) (dHat - d) -> N(0, 1/4), so alpha = 2d has SD 1 / sqrt(m). Expected source Robinson (1995, Ann. Stat. 23); NOT in the Zotero library, so VERIFY before citing.
#                It assumes the band starts at the lowest frequency.
#   sdBand     = 1 / sqrt(sum(nu^2)), nu = log(lambda_j) - mean(log(lambda_j)). Fisher information for the slope when periodogram ordinates are independent exponentials about the spectrum (derived here,
#                not taken from a paper). It allows for a band that starts above the lowest frequency.
#   sdHann     = sqrt(1.5) * sdBand. The factor 1.5 is the variance inflation usually quoted for a Hann-type taper (Hurvich & Chen 2000; VERIFY before citing). The observed/sdHann column tests it.
# Standard derived from the design: if the gap is centred with SD s, one draw is accepted with probability p = 2 Phi(tol / s) - 1, and 20 draws succeed with 1 - (1 - p)^20.
import math
import os
import sys

import numpy as np
import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
out = os.path.join(here, 'out_toleranceTable_v001')
os.makedirs(out, exist_ok=True)
pd.set_option('display.width', 240, 'display.max_columns', 40, 'display.max_rows', 200)
V2 = pd.read_csv(os.path.join(here, 'out_lowAlpha_v002', 'cells.csv'))
V2 = V2[V2.generator == 'exact'][['fs', 'N', 'alpha', 'sdAlpha', 'nEstimated']]
V3 = pd.read_csv(os.path.join(here, 'out_highAlpha_v003', 'cells.csv'))[['fs', 'N', 'alpha', 'sdAlpha', 'nEstimated']]
A = pd.concat([V2, V3[~V3.alpha.isin([2.95, 3.05])]], ignore_index=True)
ranges = [(-1.01, 0.99, '-0.75 to 0.95'), (0.99, 2.99, '1.05 to 2.95'), (2.99, 4.99, '3.05 to 4.95'), (4.99, 6.8, '5.05 to 6.75')]
lines = []


def show(title, frame):
    lines.append('\n=== ' + title + ' ===\n' + frame.to_string())


def predicted(fs, N, k, band=(1.0, 20.0)):
    M = N - k
    f = np.arange(M) * fs / M
    f = f[(f >= band[0]) & (f <= band[1])]
    nu = np.log(f) - np.log(f).mean()
    return len(f), 1 / math.sqrt(len(f)), 1 / math.sqrt((nu ** 2).sum())


rows = []
for (fs, N), q in A.groupby(['fs', 'N']):
    for lo, hi, name in ranges:
        c = q[(q.alpha > lo) & (q.alpha <= hi)]
        k = int(math.floor((c.alpha.median() + 1) / 2))
        m, sdR, sdB = predicted(fs, N, k)
        obs = c.sdAlpha.median()
        rows.append({'fs': fs, 'N': N, 'alphaRange': name, 'cells': len(c), 'mBins': m, 'sdRobinson': sdR, 'sdBand': sdB, 'sdHann': math.sqrt(1.5) * sdB, 'sdObserved': obs, 'sdObservedMax': c.sdAlpha.max(),
                     'observedOverHann': obs / (math.sqrt(1.5) * sdB), 'k_for_0.3': 0.3 / obs, 'k_at_worst_SD': 0.3 / c.sdAlpha.max()})
T = pd.DataFrame(rows)
T.to_csv(os.path.join(out, 'tolerance_table.csv'), index=False)
show('SD of estimateWhittleAlpha on exact references, and k = 0.3 / SD', T.round(3))
show('Summary per fs x N (all alpha ranges)', T.groupby(['fs', 'N']).agg(mBins=('mBins', 'median'), sdObserved=('sdObserved', 'median'), sdHann=('sdHann', 'median'), k=('k_for_0.3', 'median'), kWorst=('k_at_worst_SD', 'min')).round(3))
show('Overall', pd.DataFrame({'value': {'median k': T['k_for_0.3'].median(), 'smallest k (worst cell)': T['k_at_worst_SD'].min(), 'largest k': T['k_for_0.3'].max(), 'median observed / sdHann': T.observedOverHann.median(),
                                        'min observed / sdHann': T.observedOverHann.min(), 'max observed / sdHann': T.observedOverHann.max()}}).round(3))


# What the design implies: acceptance for a centred gap with SD s at tol 0.3
def phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


rows = []
for s in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.75, 1.0, 1.13, 1.5]:
    p = 2 * phi(0.3 / s) - 1
    rows.append({'gapSD': s, 'perDrawP': p, 'acceptedAt20': 1 - (1 - p) ** 20})
show('Design: tolerance 0.3, centred normal gap', pd.DataFrame(rows).round(4))
pReq = {c: 1 - (1 - c) ** (1 / 20) for c in [0.95, 0.99]}
show('Per-draw p needed for a given acceptance at 20 draws', pd.Series(pReq).round(4).to_frame('p'))

# Per-draw gaps on exact references that passed the segment rule, and on real residuals: is low acceptance a wide gap or an off-centre gap?
S = pd.read_csv(os.path.join(here, 'out_drawCalibration_v001', 'series.csv'))
S = S[S.segPass.astype(bool)].copy()
gc = [c for c in S.columns if c.startswith('gap_')]
G = S[gc].to_numpy(dtype=float)
S['gapMean'] = np.nanmean(G, axis=1)
S['gapSD'] = np.nanstd(G, axis=1, ddof=1)
S['p'] = np.nansum(np.abs(G) <= 0.3, axis=1) / np.maximum(np.isfinite(G).sum(1), 1)
S['name'] = np.where(S.kind == 'synthetic', S.group + ' alpha ' + S.alpha.astype(str), S.group)


def gaps(q):
    low = q[q.p < 0.21]
    return pd.Series({'series': len(q), 'medianGapSD': q.gapSD.median(), 'p90GapSD': q.gapSD.quantile(.9), 'medianAbsGapMean': q.gapMean.abs().median(), 'shareP<0.21': 100 * len(low) / len(q),
                      'lowP_medianAbsGapMean': low.gapMean.abs().median(), 'lowP_medianGapSD': low.gapSD.median()})


show('Per-draw gap, surrogate alpha minus target: spread and centring. lowP = series with per-draw p below 0.21, the value needed for 99% at 20 draws', S.groupby('name', sort=False).apply(gaps, include_groups=False).round(3))
text = '\n'.join(lines) + '\n'
with open(os.path.join(out, 'report_tables.txt'), 'w') as f:
    f.write(text)
print(text)
