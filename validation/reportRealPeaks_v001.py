# reportRealPeaks_v001.py  Aggregate counts (no trial codes) of the estimateWhittleAlpha peak check on the real residuals.
# Run:  python3 -I reportRealPeaks_v001.py   (from any folder). Reads out_realPeaks_v001/axes.csv (private, written by countRealResidualPeaks_v001.m);
# writes out_realPeaks_v001/report_counts.txt (aggregates only). A hit is peakRatio > peakThreshold (PeakFalseAlarm 1e-4 per call).
# "Near a harmonic" means the peak frequency is within 0.05 of an integer multiple of the trial's f0 (peakOverF0 within 0.05 of an integer >= 1).
import os
import sys

import numpy as np
import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
out = os.path.join(here, 'out_realPeaks_v001')
pd.set_option('display.width', 220, 'display.max_columns', 30, 'display.max_rows', 200)
T = pd.read_csv(os.path.join(out, 'axes.csv'), keep_default_na=True)
T = T[T.axis.notna()].copy()
T['scored'] = T.peakRatio.notna()
T['hit'] = T.peakInBand.astype(str).str.lower().isin(['true', '1'])
h = T.peakOverF0
T['nearHarmonic'] = T.hit & (np.abs(h - np.round(h)) <= 0.05) & (np.round(h) >= 1)
T['harmonicNumber'] = np.where(T.nearHarmonic, np.round(h), np.nan)
g = T.groupby(['K', 'dataset'])
R = pd.DataFrame({'axesScored': g.scored.sum(), 'hits': g.hit.sum(), 'hitPct': 100 * g.hit.sum() / g.scored.sum(),
                  'hitsNearHarmonic': g.nearHarmonic.sum(), 'medianPeakHz': g.apply(lambda q: q.peakHz[q.hit].median(), include_groups=False)})
lines = ['Peak check on real residuals, each axis alone, band 1-20 Hz (counts only)', R.round(2).to_string(), '',
         'All datasets: ' + ', '.join(f'K = {k}: {int(q.hit.sum())} hits of {int(q.scored.sum())} axes ({100 * q.hit.sum() / q.scored.sum():.2f}%), {int(q.nearHarmonic.sum())} near a harmonic of f0'
                                      for k, q in T.groupby('K')),
         '', 'Harmonic number of hits near a harmonic (count by K):',
         T[T.nearHarmonic].groupby(['K', 'harmonicNumber']).size().unstack(fill_value=0).to_string() if T.nearHarmonic.any() else 'none']
text = '\n'.join(lines) + '\n'
with open(os.path.join(out, 'report_counts.txt'), 'w') as f:
    f.write(text)
print(text)
