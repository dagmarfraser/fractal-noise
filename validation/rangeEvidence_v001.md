# Evidence for the measured range of shapeNoise (2026-11-18)

Every number in the RANGE section of `toolbox/shapeNoise.m`, in `CHANGELOG.md` and in `README.md` comes from the tables written by
`validation/reportRangeEvidence_v001.py`. All runs used the toolbox at commit e1d132b. Run outputs stay on the machine that made them; their md5s are listed here.

## Runs

| Run | Script (md5) | Configuration | Output (md5) |
|---|---|---|---|
| Low alpha | `checkShapeNoiseEndMatched_v002.m` (538e0f0a7c6f17f554cb83fc6af19d8f), runner `runCheckShapeNoiseEndMatched_v002.m` (9404399825d5f7ba2a7bd52c166a7a1f) | Exact references and noiseXu; nominal alpha -1 to 3.05; fs 120, 133, 240, 300 Hz; N 1000, 3000; 200 series per cell; band 1-20 Hz; master seed 1 | `out_lowAlpha_v002/cells.csv` (648a15b3b9e6a2acae0329cd5ec6d41a) |
| High alpha | `checkShapeNoiseEndMatched_v003.m` (8bf0e9036ee39fa77ad5475e69ad0ad4), runner `runCheckShapeNoiseEndMatched_v003.m` (651fc8066bffd45f5f80e3a057455bf5) | Exact references; nominal alpha 2.95 to 6.75; same fs, N and band; 200 series per cell; master seed 3 | `out_highAlpha_v003/cells.csv` (9dd08e83060dded4d17d957d4458b7f5) |
| Layers | `checkShapeNoiseLayers_v004.m` (a9fb73312b4424874946514a0697b4af), `exactRefExtended_v004.m` (764b4fda3e3d3e748b2e8796a94c881b), `classifyShapeNoise_v004.m` (ad5de7b5e30173a410bca17c30b56fad), runner `runCheckShapeNoise_v004.m` (2a5743c91c244aa01655f82c10ca3209) | Exact, noiseXu Phi 0.99 and Phi 1; nominal alpha -4 to 10; same fs, N and band; 100 series per cell; master seed 2 | `out_v004/layers_cells.csv` (d8168bce6f7974a24d5316528e35c9f5) |
| Real residuals | `runCheckShapeNoiseRealResiduals_v001.m` (5d826baa5d8e9c1f501737bd99ff875c, kept outside this repository: it loads the project's private datasets; held for v2.1.0) | Nine movement datasets, stripping K = 1 and K = 4 harmonics, each axis alone and both axes together; band 1-20 Hz | `out_realResiduals_v001/axes.csv` (91644b32aaa6e094479a0d7b441df080) |
| Tables | `reportRangeEvidence_v001.py` (dde84647f8280687631fe296f9238a5e) | Reads the four outputs above | `out_rangeEvidence_v001/report_tables.txt` (339422731cfc47b4e7e7b6a480ecda0a) |

## Exact references: share of series that fail (error or no accepted draw)

Pooled over the eight fs x N configurations. Alpha is the generator's exponent. The in-band slope is the band-effective alpha over 1-20 Hz, shown as its range across configurations.

| Alpha | In-band slope | Failed | Worst configuration |
|---|---|---|---|
| -0.75 to 2.25 (14 values) | -0.98 to 2.25 | 0 to 0.19% | 1.0% |
| 2.50 | 2.49 to 2.50 | 1.06% | 3.5% |
| 2.75 | 2.74 to 2.75 | 7.44% | 21.0% |
| 2.95 | 2.95 | 45.62% | 62.0% |
| 3.05 | 2.35 to 2.61 | 45.94% | 61.0% |
| 3.25 | 2.95 to 3.14 | 26.50% | 45.0% |
| 3.50 | 3.35 to 3.46 | 24.12% | 33.0% |
| 3.75 | 3.66 to 3.73 | 31.56% | 47.0% |
| 4.00 | 3.93 to 3.99 | 48.00% | 58.0% |
| 4.25 | 4.20 to 4.24 | 68.06% | 81.0% |
| 4.50 | 4.46 to 4.49 | 79.81% | 89.0% |
| 4.75 | 4.71 to 4.74 | 91.69% | 97.5% |
| 4.95 to 6.75 (9 values) | 4.31 to 6.74 | 95.75 to 99.69% | 100% |

Failure follows the generator's exponent and not the in-band slope. Alpha 3.05 has an in-band slope below 2.61 and fails as often as alpha 2.95, so a range cannot be stated in the in-band slope alone.
No exact series exists at odd-integer alpha (H would be 0 or 1), so those values have no known-truth check.

## Other references

| Reference | Result | Source |
|---|---|---|
| noiseXu Phi 0.99, accepted by in-band slope | 100% at 2-3.5, 99.9% at 3.5-4, 98.2% at 4-4.5, 93.0% at 4.5-5, 82.5% at 5-5.5, 70.2% at 5.5-6, 43.8% at 6-7 | Layers run |
| Real residuals, stripping K = 1 harmonic | 7948 axes, 97.0% accepted; 93.0%, 98.7%, 98.5% and 93.7% for residual alpha below 3, 3-4, 4-5 and above 5; 0 noSegment, 3 segmentAlpha, 235 no accepted draw | Real-residual run |
| estimateWhittleAlpha on exact references | Worst absolute bias 0.079 against the in-band slope for alpha -1.5 to 9.5, at most 1.0% of series with an estimator error; 44.6% errors at alpha -2, 32.8% at alpha 10, at least 99.9% at alpha -2.5 and below (AlphaBounds [-2 10]) | Layers run |

## The batch warning on real data

`summariseShapeNoiseBatch` warns only when acceptance is significantly below a 95% margin (exact one-sided binomial, level 0.01, at least 90 channels). The 95% is a conventional margin chosen after the real-data acceptance results were seen; it is not an independent derivation. The minimum of 90 channels is the smallest batch from which the test keeps at least 80% power to detect a true acceptance of 85% (exact binomial, checked for batches up to 419 channels). The 99% design standard comes from exact references and is reported for information only.

Per-dataset acceptance on the real datasets is held for the BRM paper and v2.1.0 (2026-11-21) and kept outside this repository.
