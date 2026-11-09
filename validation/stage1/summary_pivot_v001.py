"""summary_pivot_v001.py  Pivot of stage-1 summary.csv medians at one N for selected alphas.
Usage: python3 -I summary_pivot_v001.py <runDir> [N] [alphas comma list]
Prints spectral estimators (fixed band) and band-free estimators (band 'none'), plus nFail."""
import sys, os, pandas as pd
run = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 4096
al = [float(a) for a in sys.argv[3].split(",")] if len(sys.argv) > 3 else [2.5, 4.0, 4.5, 5.5, 6.5]
T = pd.read_csv(os.path.join(run, "summary.csv"))
T = T[(T.N == N) & T.alpha.isin(al)]
pd.set_option("display.width", 200)
for band in ["fixed", "none"]:
    s = T[T.band == band]
    print(f"== band {band}, N {N}: median alpha-hat"); print(s.pivot_table(index="estimator", columns="alpha", values="median").round(2).to_string())
    print(f"== band {band}, N {N}: nFail"); print(s.pivot_table(index="estimator", columns="alpha", values="nFail").to_string())
