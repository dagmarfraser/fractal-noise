# compareRealResidualsMinL120_v001.py  Which real-residual axes change outcome between shapeNoise MinLength 32 and 120? Aggregates only, no trial codes.
# Run:  python3 -I compareRealResidualsMinL120_v001.py   (from any folder). Reads out_realResiduals_v001/axes.csv (MinLength 32) and
# out_realResiduals_minL120_v001/axes.csv (MinLength 120), both private (trial level). Writes out_compareRealResidualsMinL120_v001/report.txt (counts only).
import os
import sys

import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
out = os.path.join(here, 'out_compareRealResidualsMinL120_v001')
if os.path.exists(os.path.join(out, 'report.txt')):
    raise SystemExit('report exists; use a new version')
os.makedirs(out, exist_ok=True)
k = ['dataset', 'trial', 'K', 'axis']
a = pd.read_csv(os.path.join(here, 'out_realResiduals_v001', 'axes.csv'), low_memory=False)
b = pd.read_csv(os.path.join(here, 'out_realResiduals_minL120_v001', 'axes.csv'), low_memory=False)
m = a.merge(b, on=k, suffixes=('_32', '_120'), validate='one_to_one')
if len(m) != len(a) or len(m) != len(b):
    raise SystemExit('the two runs do not cover the same axes')
ch = m[m.cls_32 != m.cls_120]
lines = [f'Rows compared: {len(m)} (every row of axes.csv: nine datasets, K = 1 and 4, one row per axis; rows with N = 0 included).',
         f'Axes that moved between accepted and failed: {int(((m.cls_32 == "ok") != (m.cls_120 == "ok")).sum())}.',
         f'Axes whose failure class changed: {len(ch)}.', '',
         ch.groupby(['K', 'cls_32', 'cls_120']).size().rename('axes').to_frame().to_string(), '',
         'Record length N of the changed axes (samples), counted: ' + ch.N_120.value_counts().sort_index().to_string().replace('\n', '; ')]
text = '\n'.join(lines) + '\n'
with open(os.path.join(out, 'report.txt'), 'w') as f:
    f.write(text)
print(text)
