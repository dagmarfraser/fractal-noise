"""summarise_dfa45_v001.py  DFA-4/5 supplementary pass: median alpha-hat and failures by estimator x alpha at selected N.
Usage: python3 -I summarise_dfa45_v001.py <runDir> [N comma list]"""
import sys, os, pandas as pd
run = sys.argv[1]; Ns = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [1024, 4096]
M = pd.read_csv(os.path.join(run, "manifest.csv"))[["cell_id", "alpha", "N"]]
E = pd.concat([pd.read_csv(os.path.join(run, f), na_values=["NA", "nan", "NaN"]) for f in ["estimates_python_dfa45.csv", "estimates_R_dfa45.csv"]])
E = E.merge(M, on="cell_id", how="left")
if E.alpha.isna().any(): raise RuntimeError("estimate rows without a manifest cell")
ok = (E.status == "ok") & E.alpha_hat.notna()
print(f"rows {len(E)}, FAIL {int((~ok).sum())}")
pd.set_option("display.width", 220)
for N in Ns:
    s = E[ok & (E.N == N) & (E.alpha >= 2.5)]
    print(f"== N {N}: median alpha-hat"); print(s.pivot_table(index="estimator", columns="alpha", values="alpha_hat", aggfunc="median").round(2).to_string())
