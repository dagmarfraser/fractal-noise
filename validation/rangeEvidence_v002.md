# Evidence for the measured range of shapeNoise, MinLength 120 (2026-11-21)

Every number in the RANGE section of `toolbox/shapeNoise.m`, in `CHANGELOG.md` and in `README.md` comes from the tables written by
`validation/reportRangeEvidence_v002.py`. These are the runs of `rangeEvidence_v001.md` (kept, unedited) repeated with the new `shapeNoise` default
`MinLength = 120` (it was 32). Toolbox: commit 8aa06fe plus the uncommitted MinLength change recorded in `provenance_minL120default_v001.txt`.
Run outputs stay on the machine that made them; their md5s are listed here.

## Runs

| Run | Script (md5) | Configuration | Output (md5) |
|---|---|---|---|
| Low alpha | `checkShapeNoiseEndMatched_v005.m` (e4f9c8c9674059228839013cf8b4693c), runner `runCheckShapeNoiseMinLength_v001.m` (dceea5627525ff4c3d48be745bf3879e) | The v002 grid, `ShapeArgs = {"MinLength", 120}`; master seed 1 | `out_lowAlpha_minL120_v001/cells.csv` (19eec22dbb28a468bcbaea9540bb7029) |
| High alpha | `checkShapeNoiseEndMatched_v006.m` (b2af6cbb1ed096f0de05281533a22503), same runner | The v003 grid, `ShapeArgs = {"MinLength", 120}`; master seed 3 | `out_highAlpha_minL120_v001/cells.csv` (2f16fca3e3a2169612441bf3a0fa5f1e) |
| Layers | `checkShapeNoiseLayers_v004.m` and helpers (unchanged), runner `runCheckShapeNoiseLayersMinL120_v001.m` (6b97d0b1ac2ed19b26756e092805def2) | The v004 grid under the MinLength 120 default; master seed 2 | `out_layers_minL120_v001/layers_cells.csv` (25fa24571e7cb481e86e704b1efab5af) |
| Real residuals | `runCheckShapeNoiseRealResidualsMinL120_v001.m` (ddd7ecbd07b2a3f34c84ac17e8b63975, kept outside this repository: it loads the project's private datasets; held for v2.1.0) | As v001, under the MinLength 120 default | `out_realResiduals_minL120_v001/axes.csv` (524f3772bbf44979a28186c54f260a83; private, trial level) |
| Tables | `reportRangeEvidence_v002.py` (7ad3faed541cc9b7a43d9c7299b96bf2) | Reads the four outputs above | `out_rangeEvidence_v002/report_tables.txt` (b7c391250ace4185c54ba03dde93f53a) |
| Comparison, grids | `compareMinLength120_v001.py` (d6b647e4879982be9aa5fd28d8ba22ee) | MinLength 32 against 120, cell by cell | `out_compareMinLength120_v001/report.txt` (1e432f6423cf690b51fbb8d4f598db3d) |
| Comparison, real | `compareRealResidualsMinL120_v001.py` (d3ccca211a5cf3bac018dc2f12e1929d) | MinLength 32 against 120, counts only | `out_compareRealResidualsMinL120_v001/report.txt` (f79107ca911a25092eef46139b29adae) |

Checks before the reruns: v005 and v006 with `ShapeArgs = {}` gave exactly what v002 and v003 give, and both reproduced the saved 18 Nov smoke outputs
(log of `runCheckShapeNoiseMinLength_v001.m`). The reruns under the new default start with a guard that `shapeNoise` refuses a 119-sample record.
The reference-estimator columns are identical between the MinLength 32 and 120 runs, as they should be, since they do not involve `shapeNoise`.

## Exact references: share of series that fail (error or no accepted draw)

Pooled over the eight fs x N configurations; MinLength 32 (v001) in brackets.

| Alpha | Failed | Worst configuration |
|---|---|---|
| -0.75 to 2.25 (14 values) | 0 to 0.19% (unchanged) | 1.0% |
| 2.50 | 1.31% (1.06%) | 3.5% |
| 2.75 | 8.06% (7.44%) | 21.0% |
| 2.95 to 4 | 25.25 to 50.50% (24.12 to 48.00%) | 63.5% |
| 4.25 to 4.5 | 69.44 to 81.19% (68.06 to 79.81%) | 89.0% |
| 4.75 | 92.25% (91.69%) | 97.5% |
| 4.95 to 6.75 | 95.81 to 99.69% (95.75 to 99.69%) | 100% |

The largest change is +2.6 percentage points (alpha 4.0, N 1000). What changes is the kind of failure: on the high grid, segmentAlpha failures caused by a
segment too short for `estimateWhittleAlpha` fall from 2674 to 0, and noSegment rises from 18389 to 21503 (default call). Surrogate scores on the series
that succeed are unchanged at low alpha (largest change 0.026 in surrogate bias); at high alpha, where few series succeed, one cell's median bias moved by 0.20.

## Other references

| Reference | Result (MinLength 32 in brackets) | Source |
|---|---|---|
| noiseXu Phi 0.99, accepted by in-band slope | 100% at 2-3.5, 99.9% at 3.5-4, 98.1% (98.2%) at 4-4.5, 92.3% (93.0%) at 4.5-5, 81.2% (82.5%) at 5-5.5, 69.9% (70.2%) at 5.5-6, 43.4% (43.8%) at 6-7 | Layers run |
| Real residuals, stripping K = 1 harmonic | 7948 axes, 97.0% accepted (unchanged); 93.0%, 98.7%, 98.5% and 93.7% for residual alpha below 3, 3-4, 4-5 and above 5 (unchanged) | Real-residual run |
| Real residuals, what changed | No axis moved between accepted and failed. 10 failures changed kind over K = 1 and 4: 6 now noSegment, 4 records of 101 samples now `shapeNoise:short` | Comparison, real |
| estimateWhittleAlpha on exact references | Unchanged (it does not involve `shapeNoise`): worst absolute bias 0.079 for alpha -1.5 to 9.5 | Layers run |

## The batch warning on real data

Unchanged from v001: with MinLength 120 every dataset has the same accepted count, so the batch-warning outcome is the same.
Per-dataset acceptance on the real datasets is held for the BRM paper and v2.1.0 (2026-11-21) and kept outside this repository.
