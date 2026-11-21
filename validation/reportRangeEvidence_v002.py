# reportRangeEvidence_v002.py  reportRangeEvidence_v001 on the reruns with the shapeNoise MinLength default of 120 (2026-11-21); v001 kept, unedited.
# Only the input and output folders differ from v001. The numbers behind the RANGE statements in the shapeNoise help, README and CHANGELOG.
# Run:  python3 -I reportRangeEvidence_v001.py   (from any folder). Writes out_rangeEvidence_v002/report_tables.txt and prints the same text.
# Reads only existing run outputs in this folder:
#   out_lowAlpha_minL120_v001/cells.csv  (the v002 grid, MinLength 120) exact references and noiseXu, nominal alpha -1 to 3.05, fs 120/133/240/300 x N 1000/3000, 200 series per cell
#   out_highAlpha_minL120_v001/cells.csv (the v003 grid, MinLength 120) exact references, nominal alpha 2.95 to 6.75, same fs x N, 200 series per cell
#   out_layers_minL120_v001/layers_cells.csv (the v004 layer grid, MinLength 120) exact, noiseXu Phi 0.99 and Phi 1, nominal alpha -4 to 10, same fs x N, 100 series per cell
#   out_realResiduals_minL120_v001/axes.csv (MinLength 120) real residuals, stripping K = 1 and K = 4 harmonics
# Exact references are binned by band-effective alpha: the in-band slope of the analytic spectrum over 1-20 Hz, which is what the estimator targets. A failure is an error (noSegment, segmentAlpha) or no accepted draw.
import math
import os
import sys

import numpy as np
import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
out = os.path.join(here, 'out_rangeEvidence_v002')
os.makedirs(out, exist_ok=True)
pd.set_option('display.width', 220, 'display.max_columns', 40, 'display.max_rows', 200)
lines = []


def show(title, frame):
    lines.append('\n=== ' + title + ' ===\n' + frame.to_string())


V2 = pd.read_csv(os.path.join(here, 'out_lowAlpha_minL120_v001', 'cells.csv'))
V2 = V2[V2.generator == 'exact'].copy()
V2['fail'] = V2.errShapeNoise + V2.nNoAcceptedDraw
V2['noSegSegAlpha'] = V2.errShapeNoise
V3 = pd.read_csv(os.path.join(here, 'out_highAlpha_minL120_v001', 'cells.csv'))
V3 = V3[V3.alpha > 3.05].copy()
V3['fail'] = V3.nRep - V3.okDefault
V3['noSegSegAlpha'] = V3.noSegmentDefault + V3.segmentAlphaDefault
cols = ['fs', 'N', 'alpha', 'alphaEffective', 'nRep', 'fail', 'noSegSegAlpha', 'biasEffective', 'sdAlpha', 'errEstimateRef']
A = pd.concat([V2[cols], V3[cols]], ignore_index=True)
A['cfg'] = A.fs.astype(str) + '/' + A.N.astype(str)
edges = [-1.01, 0, 1, 2, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 5.5, 7]
names = ['<0', '0-1', '1-2', '2-2.5', '2.5-2.75', '2.75-3', '3-3.5', '3.5-4', '4-4.5', '4.5-5', '5-5.5', '>=5.5']
A['effBin'] = pd.cut(A.alphaEffective, edges, right=False, labels=names)


def pooled(q):
    c = q.groupby('cfg').apply(lambda x: 100 * x.fail.sum() / x.nRep.sum(), include_groups=False)
    return pd.Series({'cells': len(q), 'series': q.nRep.sum(), 'failPct': 100 * q.fail.sum() / q.nRep.sum(), 'lowestCfgFailPct': c.min(), 'highestCfgFailPct': c.max(), 'worstCellFailPct': (100 * q.fail / q.nRep).max(),
                      'noSegOrSegAlphaPct': 100 * q.noSegSegAlpha.sum() / q.nRep.sum(), 'worstAbsBias': q.biasEffective.abs().max(), 'estimatorErrors': q.errEstimateRef.sum()})


show('Exact references by band-effective alpha (v002 for nominal alpha up to 3.05, v003 above)', A.groupby('effBin', observed=True).apply(pooled, include_groups=False).round(3))
show('Exact references by nominal alpha (failure percent, pooled over the 8 configurations)', A.groupby('alpha').apply(lambda q: pd.Series({'alphaEffectiveMin': q.alphaEffective.min(), 'alphaEffectiveMax': q.alphaEffective.max(), 'failPct': 100 * q.fail.sum() / q.nRep.sum(),
                                                                                                                       'highestCfgFailPct': q.groupby('cfg').apply(lambda x: 100 * x.fail.sum() / x.nRep.sum(), include_groups=False).max()}), include_groups=False).round(2))
L = pd.read_csv(os.path.join(here, 'out_layers_minL120_v001', 'layers_cells.csv'))
L = L[L.alphaEffective >= 2].copy()
L['effBin'] = pd.cut(L.alphaEffective, [2, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 11], right=False, labels=['2-3', '3-3.5', '3.5-4', '4-4.5', '4.5-5', '5-5.5', '5.5-6', '6-7', '>=7'])
show('Accepted percent by band-effective alpha, v004 layer grid (100 series per cell)', L.groupby(['generator', 'effBin'], observed=True).apply(lambda x: 100 * x.okShape.sum() / x.nGenerated.sum(), include_groups=False).unstack(0).round(1))
L2 = pd.read_csv(os.path.join(here, 'out_layers_minL120_v001', 'layers_cells.csv'))
E = L2[L2.generator == 'exact']
show('Estimator on exact references by nominal alpha, v004 (percent of series with an estimator error; worst absolute bias against band-effective alpha over configurations)',
     E.groupby('alpha').apply(lambda x: pd.Series({'estimatorErrorPct': 100 * (1 - x.nEstimated.sum() / x.nGenerated.sum()), 'worstAbsBias': x.biasEffective.abs().max()}), include_groups=False).round(3))
q = E[(E.alpha >= -1.5) & (E.alpha <= 9.5)]
show('Estimator summary, nominal alpha -1.5 to 9.5', pd.DataFrame({'value': {'worstAbsBias': q.biasEffective.abs().max(), 'highestEstimatorErrorPctByAlpha': (100 * (1 - q.groupby('alpha').nEstimated.sum() / q.groupby('alpha').nGenerated.sum())).max()}}).round(3))
R = pd.read_csv(os.path.join(here, 'out_realResiduals_minL120_v001', 'axes.csv'))
R = R[R.N > 0].copy()
R['bin'] = pd.cut(R.alphaResidual, [-np.inf, 3, 4, 5, np.inf], right=False, labels=['<3', '3-4', '4-5', '>=5'])
show('Real residuals, each axis alone', R.groupby('K').apply(lambda x: pd.Series({'axes': len(x), 'acceptedPct': 100 * (x.cls == 'ok').mean(), 'noSegment': (x.cls == 'noSegment').sum(), 'segmentAlpha': (x.cls == 'segmentAlpha').sum(), 'noAccepted': (x.cls == 'noAccepted').sum()}), include_groups=False).round(2))
show('Real residuals, accepted percent by residual alpha', R.pivot_table(index='K', columns='bin', values='cls', aggfunc=lambda s: 100 * (s == 'ok').mean(), observed=True).round(1))


def lowerTail(x, n, q=0.99):
    # exact P(X <= x), X ~ Binomial(n, q), on the log scale; independent of summariseShapeNoiseBatch
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(q) + (n - i) * math.log1p(-q) for i in range(x + 1)]
    return min(1.0, sum(math.exp(v) for v in logs))


rows = []
for (ds, K), q in R.groupby(['dataset', 'K']):
    n, x = len(q), int((q.cls == 'ok').sum())
    p95 = lowerTail(x, n, 0.95)
    p99 = lowerTail(x, n, 0.99)
    rows.append({'dataset': ds, 'K': K, 'channels': n, 'accepted': x, 'acceptedPct': 100 * x / n, 'pAgainst95': p95, 'wouldWarn': (n >= 90) and (p95 < 0.01), 'tested': n >= 90, 'pAgainst99': p99, 'belowStandard99': p99 < 0.01})
show('What the batch rule would do on each real dataset (warning: margin 0.95, level 0.01, at least 90 channels; 0.99 reported for information; each axis alone)', pd.DataFrame(rows).set_index(['K', 'dataset']).sort_index().round({'acceptedPct': 2}))
R['L'] = R.Lfrac * R.N
rows = []
for K, q in R[R.L.notna()].groupby('K'):
    rows.append({'K': K, 'axesWithSegment': len(q), 'L<100': int((q.L < 99.5).sum()), 'L 100-109': int(((q.L >= 99.5) & (q.L < 109.5)).sum()), 'L 110-119': int(((q.L >= 109.5) & (q.L < 119.5)).sum()), 'L<200': int((q.L < 199.5).sum()), 'shortestL': int(round(q.L.min()))})
show('Real residuals: chosen segment length L, each axis alone', pd.DataFrame(rows).set_index('K'))
text = '\n'.join(lines) + '\n'
with open(os.path.join(out, 'report_tables.txt'), 'w') as f:
    f.write(text)
print(text)
