"""Build the markdown tables in V1L_prototype_results_v001.md from a v1l_leakage_v001.py output CSV.
Usage: python v1l_tables_v001.py <v1l.csv>  (prints markdown; median alpha_hat minus band-effective alpha per cell)"""
import sys, csv, collections, statistics as st
d = collections.defaultdict(list); eff = {}
for r in csv.DictReader(open(sys.argv[1])):
    key = (float(r["alpha"]), int(r["N"]), r["band"]); d[key + (r["method"],)].append(float(r["alpha_hat"])); eff[key] = float(r["alpha_eff"])
meths = ["M0_raw", "M0a_adaptive", "M1_prewhiten_oraclek", "M1e_prewhiten_estk", "M2_more_tapers", "M3_lin_detrend", "M4_end_match", "M5_hann"]
short = ["raw", "adaptive", "prewhiten (oracle k)", "prewhiten (est. k)", "more tapers", "lin. detrend", "end-match", "Hann"]
for b in ["record", "fixed"]:
    print(f"\n**{b} band** (median α̂ − band-effective α; {len(d[next(iter(d))])} reps per cell)\n")
    print("| α | N | α_eff | " + " | ".join(short) + " |"); print("|" + "---|" * (3 + len(short)))
    for key in sorted(eff):
        if key[2] != b: continue
        print(f"| {key[0]} | {key[1]} | {eff[key]:.2f} | " + " | ".join(f"{st.median(d[key + (m,)]) - eff[key]:+.2f}" for m in meths) + " |")
