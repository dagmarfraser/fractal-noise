"""xcheck_pivot_v001.py  Pivot crossCheckIRASAGenerator_v001 summary: median alpha-hat (iraAlphaSigma_v001) by generator x fs x N vs alpha.
Usage: python3 -I xcheck_pivot_v001.py <runDir> [estimator]"""
import sys, os, pandas as pd
run = sys.argv[1]; est = sys.argv[2] if len(sys.argv) > 2 else "iraAlphaSigma_v001"
T = pd.read_csv(os.path.join(run, "xcheck_summary_v001.csv")); T = T[T.estimator == est]
pd.set_option("display.width", 200)
print(f"== {est}: median alpha-hat"); print(T.pivot_table(index=["gen", "fs", "N"], columns="alpha", values="median").round(2).to_string())
print("== band-effective alpha (target)"); print(T.pivot_table(index=["gen", "fs", "N"], columns="alpha", values="alphaEff").round(2).to_string())
print("== sd"); print(T.pivot_table(index=["gen", "fs", "N"], columns="alpha", values="sd").round(2).to_string())
