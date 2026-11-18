# reportToleranceTable_v002.py  Observed sampling SD of estimateWhittleAlpha against the theoretical scales, with the sources checked.
# Run:  python3 -I reportToleranceTable_v002.py   (from any folder). Reads out_toleranceTable_v001/tolerance_table.csv (written by reportToleranceTable_v001.py); writes out_toleranceTable_v002/report.txt.
#
# Correction to v001. v001 used sqrt(1.5) as "the Hann factor". Hurvich & Chen (2000, J. Time Ser. Anal. 21, 155-180) give 1.5/(4m) for their own complex-valued taper applied to differenced data, not for a Hann
# (cosine bell) taper. For the cosine bell they quote Velasco (1999b) as "35/18m" (p. 159); Velasco's paper is not in the library, so that value is not used here.
# Sources read: Robinson (1995, Ann. Stat. 23, 1630-1661), Theorem 2: m^(1/2) (H_hat - H_0) -> N(0, 1/4). With alpha = 2H - 1 this gives SD(alpha) = 1 / sqrt(m) for an untapered estimate.
# sdBand = 1 / sqrt(sum(nu^2)), nu = log(lambda_j) - mean(log(lambda_j)), allows for a band that does not start at the lowest frequency (derived in v001, not taken from a paper).
import os
import sys

import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
out = os.path.join(here, 'out_toleranceTable_v002')
os.makedirs(out, exist_ok=True)
T = pd.read_csv(os.path.join(here, 'out_toleranceTable_v001', 'tolerance_table.csv'))
a = T.sdObserved / T.sdRobinson
b = T.sdObserved / T.sdBand
R = pd.DataFrame({'median': [a.median(), b.median(), (b ** 2).median()], 'min': [a.min(), b.min(), (b ** 2).min()], 'max': [a.max(), b.max(), (b ** 2).max()]},
                 index=['observed SD / (1 / sqrt(m))', 'observed SD / sdBand', 'implied variance factor over sdBand^2'])
text = f'Cells: {len(T)} (fs x N x alpha range, exact references)\n' + R.round(3).to_string() + '\n'
with open(os.path.join(out, 'report.txt'), 'w') as f:
    f.write(text)
print(text)
