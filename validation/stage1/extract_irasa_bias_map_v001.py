"""extract_irasa_bias_map_v001.py  Python twin of extractIRASABiasMap_v001.m (knife item K-IRASA1).

Same table as the MATLAB function, for use while MATLAB is busy and for non-MATLAB users.
Usage: python3 extract_irasa_bias_map_v001.py <runDir> [estimator] [alo ahi] [N1,N2,...]
Writes <runDir>/biasMap_matlab_py_v001.csv. Median is pandas' (same as MATLAB median); IQR uses
numpy linear quantiles, which can differ slightly from MATLAB iqr().
CAVEAT: stage-1 runs at fs = 100 Hz; the pre-reg IRASA runs at 240 Hz. Indicative only.
"""
import sys, os
import numpy as np, pandas as pd

run = sys.argv[1]
est = sys.argv[2] if len(sys.argv) > 2 else "matlab_HB_irasa"
alo, ahi = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else (2.5, 6.75)
Ns = [int(x) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else [2048, 4096]

E = pd.read_csv(os.path.join(run, "estimates_matlab.csv"), na_values=["NA", "nan", "NaN"])
M = pd.read_csv(os.path.join(run, "manifest.csv"))
miss = ~E.cell_id.isin(M.cell_id)
if miss.any():
    raise RuntimeError(f"{miss.sum()} estimate rows have no manifest cell")
E = E.merge(M, on="cell_id", how="left")
E["alphaEff"] = np.nan
for b in E.band.unique():
    col = f"alpha_eff_{b}"
    if col not in M.columns:
        raise RuntimeError(f"manifest lacks {col}")
    E.loc[E.band == b, "alphaEff"] = E.loc[E.band == b, col]
E["ok"] = (E.status == "ok") & np.isfinite(E.alpha_hat)

def agg(g):
    a = g.alpha_hat[g.ok]
    return pd.Series(dict(n=int(g.ok.sum()), nFail=int((~g.ok).sum()), median=a.median(),
                          IQR=a.quantile(.75) - a.quantile(.25), alphaEff=g.alphaEff.iloc[0]))

T = E.groupby(["estimator", "band", "N", "alpha"]).apply(agg, include_groups=False).reset_index()
T["bias"] = T["median"] - T.alpha
T["biasEff"] = T["median"] - T.alphaEff
T["dMed"] = T.groupby(["estimator", "band", "N"])["median"].diff()
T.to_csv(os.path.join(run, "biasMap_matlab_py_v001.csv"), index=False)

S = T[(T.estimator == est) & T.alpha.between(alo, ahi) & T.N.isin(Ns)]
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 500)
print(f"{est}, alpha {alo}-{ahi}, N {Ns} (fs 100 Hz: indicative only)")
print(S[["band", "N", "alpha", "n", "nFail", "median", "IQR", "bias", "biasEff", "dMed"]].round(3).to_string(index=False))
nm = S[S.dMed <= 0]
if len(nm):
    print(f"WARNING: {len(nm)} non-monotone steps in the slice")
