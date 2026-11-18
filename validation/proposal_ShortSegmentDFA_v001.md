# Proposal: the DFA pre-estimate on short segments in estimateWhittleAlpha (2026-11-18)

Status: proposal only. No toolbox code has changed for this.

## The problem

`estimateWhittleAlpha` chooses its differencing order k from a DFA pre-estimate on 15 log-spaced box sizes between 10 and floor(n/10), rounded and de-duplicated.

| Series length n | Distinct box sizes | Behaviour |
|---|---|---|
| Below 100 | none used | Error `estimateWhittleAlpha:short` |
| 100 to 109 | 1 | The slope is a line through one point: MATLAB warns `polyfit:PolyNotUnique` and the pre-estimate is arbitrary |
| 110 to 119 | 2 | A slope from two points, erratic |
| 120 and above | 3 or more | An ordinary regression, on few points at first |

## Evidence

- `checkEstimatorShortDFA_v004.m` (md5 65bc0709c80d9080e724a5feff682c0e), fs 120 Hz, one random walk per length (true alpha 2). The pre-estimate was -1.68, -1.75 and -1.95 at n = 100, 105 and 109, which sets k = 0. At n = 110, 111 and 119 it was -2.95, 3.56 and 4.11. At n = 120 it was 2.36. One series per length shows the behaviour but does not measure its size.
- The warning appeared 917 times in the v004 layer run (`full_run_stdout_v004.txt`), because segments this short reach the estimator through `shapeNoise`.
- The 14 Nov end-mismatch report recorded 1238 of the same warnings in one dataset's section, with the source untraced. This is the likely source.
- Real residuals, each axis alone (`reportRangeEvidence_v001.py`, table "chosen segment length"):

| Stripping | Axes with a segment | L below 100 | L 100 to 109 | L 110 to 119 | L below 200 | Shortest L |
|---|---|---|---|---|---|---|
| K = 1 | 7945 | 0 | 1 | 1 | 9 | 101 |
| K = 4 | 7944 | 0 | 0 | 1 | 8 | 112 |

The problem is rare on real residuals: 2 of 7945 axes at K = 1. It is more common on steep break-free synthetic series, whose segments are short. The saved outputs do not hold segment lengths per series, so that rate needs a rerun.

## Candidate remedies

| Remedy | What it does | What it would change in the validated results |
|---|---|---|
| A. Raise the minimum length to 120 samples, using the existing error `estimateWhittleAlpha:short` | Segments of 100 to 119 samples stop with a loud error, and `shapeNoise` reports `segmentAlpha` | Real residuals at K = 1: 2 axes move from scored to `segmentAlpha`. Exact references: failure rates can only rise, most likely above alpha 2.5 where segments are short. A rerun of v002 and v003 is needed to state by how much |
| B. Document the minimum only | Help text states that fewer than 120 samples give an unreliable pre-estimate | No numbers change, and the arbitrary k remains |
| C. Take k from the whole series when the segment is short | `shapeNoise` passes `K` from an estimate on the full residual to the segment's estimate | The target alpha and the draws change for short segments only. It adds a code path and needs its own tests, and the whole-series k is itself an estimate, not a true value |
| D. Warn with a new identifier, without stopping | `estimateWhittleAlpha:shortDFA` warns for n below 120 | No numbers change. The user can see it, but the estimate is still used |

## Recommendation

Adopt A, with 120 samples as the minimum, after a rerun of the v002 and v003 grids confirms the change in failure rates. It is the smallest change and it fails loudly, as the toolbox does elsewhere. First measure the pre-estimate's error against series length on exact references, because 120 rests on the box-size count and on one series per length.
