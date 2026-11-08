"""Check of the Kiyono (2015) DFA ceiling on exact V0 series: median DFA-m alpha-hat (m = 1..5), 20 reps, seed 1.
Numbers quoted in thread 2026-11-08 ~11:00 (e.g. DFA4/5 track 6.5 at N = 4096). Usage: python checkDFAceiling_v001.py"""
import numpy as np, statistics as st
from exact_ref_v001 import exact_ref_series
from v1l_leakage_v001 import dfa_alpha
rng = np.random.default_rng(1)
for a in [4.0, 5.5, 6.5]:
    for N in [256, 4096]:
        X, _ = exact_ref_series(N, a, 20, rng)
        print(a, N, [round(st.median([dfa_alpha(X[:, r], m) for r in range(20)]), 2) for m in (1, 2, 3, 4, 5)])
