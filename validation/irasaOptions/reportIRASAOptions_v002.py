# reportIRASAOptions_v002.py  estimateIRASA bias under the fixed default (Resampling="irasa") against legacy, from checkIRASAOptions_v002.
# Run:  python3 -I reportIRASAOptions_v002.py   (from any folder). Reads runs/v002/checkIRASAOptions_v002.csv; writes runs/v002/report.txt and prints it.
# Bias = median(alphaHat) - alpha per cell (10 series), as in v001 and the estimateIRASA help RANGE. A failure is a thrown error (alphaHat NaN).
# The help RANGE figures (legacy, measured before v2.0.0) are listed beside the legacy rows here so the rerun can be checked against them.
import os
import sys

import pandas as pd

here = os.path.dirname(os.path.abspath(sys.argv[0]))
run = os.path.join(here, 'runs', 'v002')
pd.set_option('display.width', 230, 'display.max_columns', 40, 'display.max_rows', 300)
T = pd.read_csv(os.path.join(run, 'checkIRASAOptions_v002.csv'), keep_default_na=False, na_values=['NaN', ''])
T['alphaHat'] = pd.to_numeric(T.alphaHat, errors='coerce')
lines = []


def show(title, frame):
    lines.append('\n=== ' + title + ' ===\n' + frame.to_string())


g = T.groupby(['generator', 'fs', 'alpha', 'method', 'resampling'])
S = pd.DataFrame({'bias': g.alphaHat.median() - g.alpha.first(), 'nFail': g.alphaHat.apply(lambda v: int(v.isna().sum())), 'n': g.size()}).reset_index()
S['col'] = S.method + ('/' + S.resampling).where(S.resampling != 'none', '')
P = S.pivot_table(index=['generator', 'fs', 'alpha'], columns='col', values='bias').round(2)
order = [c for m in ['irasa_pmtm', 'irasa_welch', 'irasa_pmtm_pw', 'irasa_welch_pw'] for c in (m + '/irasa', m + '/legacy')] + ['whittle']
show('Bias (median alphaHat - alpha), 10 series per cell; method/resampling', P[order])
F = S[S.nFail > 0][['generator', 'fs', 'alpha', 'col', 'nFail', 'n']]
show('Cells with failures (thrown errors)', F if len(F) else pd.DataFrame({'none': ['no failures']}))
W = S.groupby(['col', 'generator']).bias.agg(lambda v: v.abs().max()).unstack().round(2)
show('Worst absolute bias over alpha 1.5-5.5 and fs 100/240, per method/resampling and generator', W.loc[[c for c in order if c in W.index]])
Dd = S[S.resampling != 'none'].pivot_table(index=['generator', 'fs', 'alpha', 'method'], columns='resampling', values='bias')
Dd['absDiff'] = (Dd['irasa'] - Dd['legacy']).abs()
show('Largest |bias(irasa) - bias(legacy)| per method, over all generators, fs and alpha', Dd.groupby('method').absDiff.max().round(2).to_frame())
show('Worst absolute bias with Prewhiten, both PSDs, per resampling', S[S.method.str.endswith('_pw')].groupby('resampling').bias.agg(lambda v: v.abs().max()).round(2).to_frame())
claims = [('exact', 240, 3.5, 'irasa_pmtm', '+0.5'), ('exact', 240, 4.5, 'irasa_pmtm', '+1.7 to +2.1'), ('xu1', 240, 4.5, 'irasa_pmtm', '+1.7 to +2.1'),
          ('exact', 100, 5.5, 'irasa_pmtm', '-2.1'), ('xu099', 100, 5.5, 'irasa_pmtm', 'within +/-0.1 to 5.5'), ('xu099', 240, 5.5, 'irasa_pmtm', 'within +/-0.1 to 5.5')]
rows = []
for gen, fs, a, m, said in claims:
    def b(rs):
        v = S[(S.generator == gen) & (S.fs == fs) & (S.alpha == a) & (S.method == m) & (S.resampling == rs)].bias
        return round(float(v.iloc[0]), 2) if len(v) else float('nan')
    rows.append({'generator': gen, 'fs': fs, 'alpha': a, 'method': m, 'help says (legacy)': said, 'legacy now': b('legacy'), 'irasa (default)': b('irasa')})
show('The estimateIRASA help RANGE figures: legacy rerun and the new default, same series', pd.DataFrame(rows))
text = '\n'.join(lines) + '\n'
with open(os.path.join(run, 'report.txt'), 'w') as f:
    f.write(text)
print(text)
