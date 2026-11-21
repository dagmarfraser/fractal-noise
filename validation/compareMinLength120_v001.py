# compareMinLength120_v001.py  Does MinLength = 120 (default 32) change the v002 (low alpha) and v003 (high alpha) shapeNoise results?
# Run:  python3 -I compareMinLength120_v001.py   (from any folder). Writes out_compareMinLength120_v001/report.txt and prints it. Stops if that file exists.
# Reads: out_lowAlpha_v002/ and out_highAlpha_v003/ (18 Nov, MinLength 32) and out_lowAlpha_minL120_v001/ and out_highAlpha_minL120_v001/ (MinLength 120),
# each cells.csv and errors.csv. Same grids, seeds and references, so every difference comes from MinLength.
# A failure is an error (noSegment, segmentAlpha, other) or no accepted draw, as in reportRangeEvidence_v001.py.
# Check first: the reference-estimator columns do not involve shapeNoise, so they must be identical; if not, the script stops.
import os
import sys

import numpy as np
import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
out = os.path.join(here, 'out_compareMinLength120_v001')
if os.path.exists(os.path.join(out, 'report.txt')):
    raise SystemExit('report exists; use a new version')
os.makedirs(out, exist_ok=True)
pd.set_option('display.width', 230, 'display.max_columns', 40, 'display.max_rows', 300)
lines = []


def show(title, frame):
    lines.append('\n=== ' + title + ' ===\n' + (frame.to_string() if hasattr(frame, 'to_string') else str(frame)))


def load(folder):
    return (pd.read_csv(os.path.join(here, folder, 'cells.csv')), pd.read_csv(os.path.join(here, folder, 'errors.csv')))


L0, LE0 = load('out_lowAlpha_v002')
L1, LE1 = load('out_lowAlpha_minL120_v001')
H0, HE0 = load('out_highAlpha_v003')
H1, HE1 = load('out_highAlpha_minL120_v001')
keyL, keyH = ['fs', 'N', 'generator', 'alpha'], ['fs', 'N', 'alpha']
ref = ['nEstimated', 'biasNominal', 'biasEffective', 'sdAlpha', 'rmseEffective', 'errReference', 'errEstimateRef']
for name, a, b, key in [('low', L0, L1, keyL), ('high', H0, H1, keyH)]:
    m = a.merge(b, on=key, suffixes=('_old', '_new'), validate='one_to_one')
    if len(m) != len(a) or len(m) != len(b):
        raise SystemExit(f'{name}: grids differ ({len(a)} vs {len(b)} cells)')
    for c in ref:
        x, y = m[c + '_old'].to_numpy(float), m[c + '_new'].to_numpy(float)
        if not np.all((np.isnan(x) & np.isnan(y)) | (np.abs(x - y) <= 1e-9 * np.maximum(1, np.abs(x)))):
            raise SystemExit(f'{name}: reference column {c} differs, so the grids are not the same runs')
lines.append('Check passed: reference-estimator columns identical in both grids (the references and seeds are the same).')

# failure per cell
L0['fail'] = L0.errShapeNoise + L0.nNoAcceptedDraw
L1['fail'] = L1.errShapeNoise + L1.nNoAcceptedDraw
for H in (H0, H1):
    H['fail'] = H.nRep - H.okDefault
    H['generator'] = 'exact'


def pooled(T, by):
    g = T.groupby(by)
    return pd.DataFrame({'series': g.nRep.sum(), 'failPct': 100 * g.fail.sum() / g.nRep.sum()})


def compare(old, new, by, title):
    a, b = pooled(old, by), pooled(new, by)
    c = a.join(b, lsuffix='_32', rsuffix='_120')
    c['change'] = c.failPct_120 - c.failPct_32
    show(title, c[['series_32', 'failPct_32', 'failPct_120', 'change']].round(2))
    return c


cl = compare(L0, L1, ['generator', 'alpha'], 'Low alpha (v002 grid): failure percent by generator and nominal alpha, pooled over 8 fs x N; MinLength 32 vs 120')
ch = compare(H0, H1, ['alpha'], 'High alpha (v003 grid, exact references): failure percent by nominal alpha, pooled over 8 fs x N; MinLength 32 vs 120')
both = pd.concat([L0[L0.generator == 'exact'].assign(src='v002'), H0[H0.alpha > 3.05].assign(src='v003')])
bothN = pd.concat([L1[L1.generator == 'exact'].assign(src='v002'), H1[H1.alpha > 3.05].assign(src='v003')])
compare(both, bothN, ['N', 'alpha'], 'Exact references (v002 up to 3.05, v003 above): failure percent by N and nominal alpha')

# failure classes on the high grid
cls = ['noSegmentDefault', 'segmentAlphaDefault', 'noAcceptedDefault', 'otherDefault']
k = H0.groupby('alpha')[cls].sum().join(H1.groupby('alpha')[cls].sum(), lsuffix='_32', rsuffix='_120')
show('High alpha: failure classes (counts of series, default call), MinLength 32 vs 120', k[[c + s for c in cls for s in ('_32', '_120')]])
lc = ['errShapeNoise', 'nNoAcceptedDraw']
show('Low alpha: shapeNoise errors and no-accepted-draw counts, MinLength 32 vs 120 (all generators)',
     L0.groupby('generator')[lc].sum().join(L1.groupby('generator')[lc].sum(), lsuffix='_32', rsuffix='_120'))


def ids(E, tag):
    s = E.groupby('identifier')['count'].sum().rename(tag)
    return s


show('Error identifiers, summed over the grid (high: default and checkOff calls)', pd.concat([ids(HE0, 'high_32'), ids(HE1, 'high_120')], axis=1).fillna(0).astype(int))
show('Error identifiers, low grid', pd.concat([ids(LE0, 'low_32'), ids(LE1, 'low_120')], axis=1).fillna(0).astype(int))

# surrogate quality on the series that succeed
q = ['surrBiasEffective', 'gapSelfMedian', 'LfracMedian', 'LfracMin', 'd1Median', 'd2Median', 'nDrawMedian']
rows = []
for name, a, b, key in [('low', L0, L1, keyL), ('high', H0, H1, keyH)]:
    m = a.merge(b, on=key, suffixes=('_32', '_120'))
    for c in q:
        x, y = m[c + '_32'], m[c + '_120']
        both_ok = x.notna() & y.notna()
        rows.append({'grid': name, 'column': c, 'cellsBoth': int(both_ok.sum()), 'median_32': x.median(), 'median_120': y.median(),
                     'maxAbsChange': (y - x)[both_ok].abs().max(), 'cellsChanged>1e-9': int(((y - x)[both_ok].abs() > 1e-9).sum()),
                     'cellsNaN_32': int(x.isna().sum()), 'cellsNaN_120': int(y.isna().sum())})
show('Surrogate scores per cell (medians over cells; changes on cells scored in both runs)', pd.DataFrame(rows).round(4))

# same bands as rangeEvidence_v001.md
edges = [-1.01, 2.3, 2.6, 2.8, 4.01, 4.6, 4.8, 7]
names = ['-0.75 to 2.25', '2.50', '2.75', '2.95 to 4', '4.25 to 4.5', '4.75', '4.95 to 6.75']
for T in (both, bothN):
    T['band'] = pd.cut(T.alpha, edges, labels=names)
r = pd.DataFrame({'failPct_32': both.groupby('band', observed=True).apply(lambda x: 100 * x.fail.sum() / x.nRep.sum(), include_groups=False),
                  'failPct_120': bothN.groupby('band', observed=True).apply(lambda x: 100 * x.fail.sum() / x.nRep.sum(), include_groups=False)})
r['change'] = r.failPct_120 - r.failPct_32
show('Exact references in the rangeEvidence_v001.md bands (pooled over fs x N)', r.round(2))
text = '\n'.join(lines) + '\n'
with open(os.path.join(out, 'report.txt'), 'w') as f:
    f.write(text)
print(text)
