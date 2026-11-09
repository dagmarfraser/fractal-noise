"""compareSeeds_v001.py  Cross-seed reproducibility of stage-1 summaries.
Usage: python3 -I compareSeeds_v001.py <runsDir> <outDir> run1 run2 run3 ...
Reads <runsDir>/<run>/summary.csv for each run, matches rows on estimator x band x cell_id, and writes
  <outDir>/seedSpread_v001.csv   per row: median alpha-hat per seed, spread (max-min), mean bias/biasEff, total nFail
  <outDir>/seedSpread_v001.txt   per estimator x band: rows, median and max spread, rows with spread > 0.1,
                                 total nFail, mean |biasEff| per seed (to see whether any estimator's verdict moves).
Missing rows in any seed are reported, never dropped silently."""
import sys, os, pandas as pd
runs_dir, out_dir, runs = sys.argv[1], sys.argv[2], sys.argv[3:]
os.makedirs(out_dir, exist_ok=True)
key = ["estimator", "band", "cell_id", "alpha", "N"]
frames = []
for i, r in enumerate(runs):
    T = pd.read_csv(os.path.join(runs_dir, r, "summary.csv"))[key + ["median", "bias", "biasEff", "nFail", "nOK"]]
    frames.append(T.rename(columns={c: f"{c}_s{i}" for c in ["median", "bias", "biasEff", "nFail", "nOK"]}))
M = frames[0]
for F in frames[1:]:
    M = M.merge(F, on=key, how="outer", indicator=False)
med = M[[f"median_s{i}" for i in range(len(runs))]]
missing = med.isna().any(axis=1)
M["spread"] = med.max(axis=1) - med.min(axis=1)
M["meanBias"] = M[[f"bias_s{i}" for i in range(len(runs))]].mean(axis=1)
M["meanBiasEff"] = M[[f"biasEff_s{i}" for i in range(len(runs))]].mean(axis=1)
M["nFailTotal"] = M[[f"nFail_s{i}" for i in range(len(runs))]].sum(axis=1)
M.to_csv(os.path.join(out_dir, "seedSpread_v001.csv"), index=False)
g = M.groupby(["estimator", "band"])
S = pd.DataFrame({
    "rows": g.size(),
    "spreadMedian": g.spread.median(), "spreadMax": g.spread.max(),
    "nSpreadGt0.1": g.spread.apply(lambda s: int((s > 0.1).sum())),
    "nFailTotal": g.nFailTotal.sum(),
    **{f"mean|biasEff|_s{i}": g[f"biasEff_s{i}"].apply(lambda s: s.abs().mean()) for i in range(len(runs))},
}).round(3)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
txt = (f"runs: {runs}\nrows: {len(M)}; rows missing in >=1 seed: {int(missing.sum())}\n\n" + S.to_string() +
       "\n\nTop 15 spreads:\n" + M.nlargest(15, "spread")[key + [f'median_s{i}' for i in range(len(runs))] + ['spread', 'nFailTotal']].round(3).to_string())
open(os.path.join(out_dir, "seedSpread_v001.txt"), "w").write(txt); print(txt)
