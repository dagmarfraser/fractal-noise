# fractalnoise

**Make, match and measure 1/f^alpha noise in MATLAB, including the steep spectra (alpha 3-7) that movement and other behavioural signals have and that most tools were never checked on.**

Use it when you need to:

1. **Simulate noise with a chosen spectral exponent**, for example to sweep alpha in a simulation grid or to test an analysis pipeline against known noise.
   `x = noiseXu(3000, 4.5, 1, 240);`
2. **Build a null model of a signal you measured**: remove its periodic part, then make random-phase surrogates of the residual that keep both its spectral exponent and its derivative noise at steep alpha, where ordinary Fourier phase randomisation of the whole record does not (see the [worked story](#a-worked-story-from-a-recorded-movement-to-a-null-model)).
3. **Measure the spectral exponent of a signal**, with an estimator that stays accurate at steep alpha, and separate a scale-free background from rhythmic peaks.
   `alpha = estimateWhittleAlpha(x, 240, 1, 20);`
4. **Check whether a spectrum is really a power law** or bends at a knee.
   `r = checkSpectralKnee(x, 240);`

Typical users: movement science and motor control, psychophysics, physiology, and anyone who needs coloured noise outside the white-to-brown range (alpha 0-2) that built-in MATLAB generators cover.

MIT licence. Base MATLAB for the generator, the surrogate and the Whittle estimator; Signal Processing Toolbox for the other two. Current release: see [CHANGELOG.md](CHANGELOG.md). Also on MATLAB File Exchange (184251).

> **v2.0.0.** Two fixes: `estimateIRASA` did not remove oscillatory peaks up to v1.2.0, so its exponent was a plain spectral slope; and `shapeNoise` gave a shallow surrogate about one time in five at alpha 3.5-4.5, so it now builds its surrogate from an end-matched segment of the residual (shorter than the input). Both old behaviours stay available by name. New: `harmonicTemplate` and `harmonicContent` remove and diagnose the periodic part of a record, `shapeNoiseWrapper` and `cetNoiseWrapper` do the whole path in one call, and `summariseShapeNoiseBatch` reports acceptance across many surrogates. Also changed: `shapeNoise` needs at least 120 samples by default and warns when no `AlphaBand` is given, and the wrappers require K. See [Upgrading from 1.x](#upgrading-from-1x).

## Quick start

Every code block in this README runs as written, top to bottom, in one MATLAB session (`tests/testReadmeAndExamples.m` checks this).

```matlab
rng(1);
x = noiseXu(3000, 4.5, 1.0, 240);                  % alpha 4.5, unit SD, 240 Hz
[a, k] = estimateWhittleAlpha(x, 240, 1, 20)       % exponent over 1-20 Hz, and the differencing order used
aI = estimateIRASA(x, 240, 1, 20, Prewhiten=k)     % IRASA on the same band, safe at steep alpha
[s, info] = shapeNoise(x, 240, [], AlphaBand="kinematic");   % a random-phase surrogate with the same exponent (length info.L)
r = checkSpectralKnee(x, 240);  r.verdict           % knee or no knee? (noiseXu's default has a real break near 0.0016*fs)
```

## A worked story: from a recorded movement to a null model

The question this toolbox was built for: in a recorded movement, is a pattern you see more than the noise would produce? To answer it you need noise that looks like the movement's own noise. The steps below use a synthetic stand-in for one axis of a traced ellipse, so they run anywhere; with your data, start from your own position trace.

```matlab
% 1. The record: a 1.2 Hz movement with a second harmonic, plus steep noise (240 Hz, 12.5 s).
fs = 240;  t = (0:1/fs:12.5)';
rng(1);
pos = 2*sin(2*pi*1.2*t) + 0.5*sin(2*pi*2.4*t + 0.7) + noiseXu(numel(t), 3.5, 0.05, fs, Phi=1);

% 2. How many harmonics stand out above the noise? Look at the ratios, not only the count.
[K, hc] = harmonicContent(pos, fs, 1.2);
hc.ratio(1:4)                                    % harmonics 1 and 2 stand far above the noise beside them; 3 and 4 do not, so K = 2

% 3. Strip the periodic part. What is left is the noise you want to model.
res = harmonicTemplate(pos, fs, 1.2, K);

% 4. Measure it: the exponent over 1-20 Hz, and whether the spectrum bends.
alphaRes = estimateWhittleAlpha(res, fs, 1, 20)
knee = checkSpectralKnee(res, fs);
[knee.verdict, " knee at " + knee.kneeHz + " Hz, in the 1-20 Hz band: " + knee.kneeInBand]

% 5. Many surrogates of the residual: random phases, same exponent and derivative noise.
infos = cell(1, 100);
for i = 1:100
    [sur, infos{i}] = shapeNoise(res, fs, [], AlphaBand="kinematic");
end
S = summariseShapeNoiseBatch(infos);             % warns only if acceptance is significantly low
fprintf("%d of %d surrogates accepted\n", S.nAccepted, S.nChannels)
```

What the output tells you, step by step:

- **Step 2.** `harmonicContent` counts the harmonics that stand out (K = 2 here). See [Choosing K for your data](#choosing-k-for-your-data) before you rely on one count.
- **Step 4.** The knee check finds a knee near 0.3 Hz, below the 1-20 Hz band (`kneeInBand` false). The template removes a constant and a trend in each window of four cycles, which takes out the slowest drift. A knee below the band does not bias the exponent you measured over 1-20 Hz; a knee inside the band would.
- **Step 5.** If `shapeNoise` warns `estimateWhittleAlpha:peakInBand`, a narrow peak is left in the residual (often a harmonic: raise K); `info.peakInBand` records it per channel. Each surrogate is a draw of noise like the residual's. To test a pattern, compute your statistic on the movement and on the surrogates (each added back to the template if your statistic needs the movement), and compare. `summariseShapeNoiseBatch` tells you whether the surrogates were built reliably across the batch.

Steps 3 and 5 in one call: `[sur, info] = shapeNoiseWrapper(pos, fs, F0=1.2, K=K, AlphaBand="kinematic")`. The measured range of `shapeNoise` (how steep a record it handles) is in its help, under RANGE.

## What's in it

| Function | What it does | Needs |
|---|---|---|
| `noiseXu` | Generates 1/f^alpha noise by fractional differencing (Xu 2019), alpha from -2 to 7 | Base MATLAB |
| `shapeNoise` | Surrogate of a measured signal: random phases on the amplitude spectrum of its longest end-matched segment, kept only if its exponent matches | Base MATLAB (`Method="legacy"` carriers `"pinknoise"`/`"dsp"` need Signal Processing / DSP System Toolbox) |
| `harmonicTemplate` | Windowed harmonic-regression template: removes the periodic part of a record (f0 given or estimated, K harmonics, window in cycles) and returns the residual `shapeNoise` takes; K and the window are yours to choose | Base MATLAB |
| `harmonicContent` | Diagnostic: how many harmonics of f0 stand above the noise between them (content rule); a count to look at, not a selector of record | Base MATLAB |
| `summariseShapeNoiseBatch` | Acceptance share and tries across many `shapeNoise` calls; warns when acceptance is significantly below 95% | Base MATLAB |
| `shapeNoiseWrapper`, `cetNoiseWrapper` | One call from a record to a surrogate (optional harmonic template, then `shapeNoise`); K, the number of harmonics stripped, is required and has no default (see [Choosing K for your data](#choosing-k-for-your-data)); `cetNoiseWrapper` is the worked example on the monthly Central England Temperature series, not hand movement | Base MATLAB |
| `estimateWhittleAlpha` | Spectral exponent over a band by prewhitened local Whittle; the one to use at steep alpha | Base MATLAB |
| `estimateIRASA` | IRASA (Wen & Liu 2016): spectral exponent and noise magnitude of the scale-free part, with oscillatory peaks removed | Signal Processing Toolbox; optional FieldTrip cross-check |
| `checkSpectralKnee` | Lorentzian knee against a fixed power law (Preston, Smith & Voytek 2026) | Signal Processing Toolbox |

Optional features (FieldTrip, `dsp.ColoredNoise`, `pinknoise` carriers) degrade visibly: a missing dependency gives a flagged "unavailable" result and a warning, never a silent substitute.

## Choosing K for your data

`shapeNoiseWrapper` (whenever it strips a template) and `cetNoiseWrapper` need K, the number of harmonics of f0 to remove. There is no default: a missing K is an error (`shapeNoiseWrapper:KRequired`, `cetNoiseWrapper:KRequired`). Too few harmonics leave movement in the "noise"; too many carve notches into the noise. `harmonicContent` is the per-record diagnostic for this choice. It is not a selector of record; an information-criterion selector is planned for v2.1.0.

1. **Know f0.** Give it if you know it, or let `harmonicContent` estimate it from a band (`F0Band=[fLow fHigh]`).
2. **Look at the harmonics.** `[Kc, info] = harmonicContent(x, fs, f0)`. `info.ratio(k)` is harmonic k's amplitude over the noise just beside it; Kc is the last harmonic, counting up from 1, at or above `Cut` (2). If `info.noHarmonic` is true, not even the fundamental stands out: Kc is 1 only by the floor, and you may not need a template at all.
3. **Choose K and pass it.** Kc is a starting point. With many records, look at the spread of Kc rather than one record.
4. **Check that your result does not depend on it.** Rerun with neighbouring K (for example Kc - 1, Kc, Kc + 1, using `harmonicTemplate` or the wrapper). If anything you report moves with K, report the ladder, not one value.

Worked example, the bundled Central England Temperature series (tested in `tests/testWrappers.m`):

```matlab
cet = cetNoiseWrapper(LoadOnly=true, FirstYear=1900, LastYear=2025);   % monthly, fs = 12 per year
[Kc, info] = harmonicContent(cet.value, 12, 1);   % f0 = 1 per year: ratios 73.8, 6.0, then 0.56 at the third harmonic, so Kc = 2
n = cetNoiseWrapper(FirstYear=1900, LastYear=2025, K=Kc);
```

The template is fitted in windows of `Window` cycles (4 years here), so it also removes variability slower than about 4 years. That suits a null model on shorter timescales, but it is the wrong residual for asking whether temperature has long memory: over 2 to 50 year periods the template residual looks strongly anti-persistent (alpha about -1), which is the template's doing, not the climate's. For that question, remove only each calendar month's mean and fit a low band:

```matlab
cet = cetNoiseWrapper(LoadOnly=true);   % 1659 onwards
anom = cet.value;
for m = 1:12, s = cet.month == m;  anom(s) = anom(s) - mean(anom(s)); end
alphaCET = estimateWhittleAlpha(anom, 12, 0.02, 1);   % periods of 1 to 50 years: alpha 0.39
H = (alphaCET + 1) / 2;   % 0.70
```

For a stationary series the spectral exponent, the Hurst exponent H, the long-memory parameter d, the DFA exponent a and the fractal (Hausdorff) dimension D of the series' cumulative sum are one quantity in different units: alpha = 2H - 1 = 2d = 2a - 1, and D = 2 - H (for a random-walk-like series, alpha = 2H + 1 and D = 2 - H of the series itself). This holds for self-affine Gaussian series with one scaling regime; in general D measures roughness at small scales and H memory at large scales, and they can differ (Gneiting & Schlather 2004, SIAM Review). Here H 0.70 (D 1.30) sits between published values for the Central England series and for land temperature: H 0.72-0.80 for annual CET 1659-2017 (Gonzalez Hidalgo et al. 2020, Cuadernos de Investigacion Geografica) and H about 0.65 for local and regional land records (Rypdal, Østvand & Rypdal 2013, J. Geophys. Res. Atmos.). Keep the upper band edge well below Nyquist: local Whittle assumes a power law near zero frequency, and with this record's length an upper edge at 5 cycles per year (83% of Nyquist) biased exact references by up to 0.45 (`validation/checkWhittleBandEdge_v001.m`; the band 0.02-1 kept every bias within 0.03).

## Which estimator to use

| Your signal | Use |
|---|---|
| Steep (alpha above about 3), or unknown | `estimateWhittleAlpha` |
| Has rhythmic peaks you want ignored (or may have: a harmonic left after the template, tremor, mains hum) | `estimateIRASA`, with `Prewhiten=k` if it may be steep. An in-band peak biases `estimateWhittleAlpha`; it warns (`estimateWhittleAlpha:peakInBand`) when it finds a narrow one |
| Both | Run both and compare: `[~, k] = estimateWhittleAlpha(x, 240, 1, 20); estimateIRASA(x, 240, 1, 20, Prewhiten=k)` |

- **`estimateWhittleAlpha`** differences the series k times (k from a DFA-3 pre-estimate), so low-frequency power cannot leak into the fitting band, then fits a local Whittle likelihood over your band. On exact reference series it is within 0.08 of the in-band slope for alpha -1.5 to 9.5, and it errors rather than return a number at its search bounds (`validation/rangeEvidence_v002.md`).
- **`estimateIRASA`** resamples the signal by a set of factors so that rhythmic peaks move while the scale-free background stays put, then takes the median. Its spectra are whole-record multitaper estimates, which leak at steep alpha: without `Prewhiten`, noise that has no low-frequency break can read several units wrong. `Prewhiten=k` removes that. The fitting band must satisfy `fHigh*max(Hset) < fs/2` (with the default `Hset`, `fHigh < fs/3.8`).

## What `noiseXu` generates (read this before using alpha > 2)

`noiseXu` implements Xu's generalized Gauss-Markov recursion, x = (1 - Phi*B)^(-alpha/2) w, with Phi = 0.99 by default.

- For alpha <= 2 this is a 1/f^alpha process down to a low-frequency break at f_b = (1 - Phi)*fs/(2*pi), about 0.0016*fs. Below the break the spectrum flattens.
- For alpha > 2 it is Xu's range-limited approximation of a steeper power law. The generator is numerically stable to alpha = 7, but the series is a power law with a break.
- `Phi = 1` removes the break (plain fractional integration). Use it to test an estimator on break-free data.
- The output is rescaled to the standard deviation you ask for. At steep alpha most of that variance is slow drift, so the power inside any fixed band depends on record length and alpha as well as on the SD.

The break makes spectral estimators look better on default `noiseXu` output than on break-free data. Keep "the generator is stable" and "an estimator can recover alpha" apart.

## Where `noiseXu` sits among MATLAB options

| | Method | Alpha range | Requires |
|---|---|---|---|
| `randn` | White noise | 0 | Base MATLAB |
| `pinknoise` | Pink noise | 1 | Signal Processing Toolbox |
| `dsp.ColoredNoise` | Kasdin (1995) filter recursion | [-2, 2] | DSP System Toolbox |
| `fractalcoef` | Kasdin (1995), raw filter coefficients | (0, 2) | Sensor Fusion and Tracking Toolbox |
| `wfbm` | Fractional Brownian motion | (1, 3) | Wavelet Toolbox |
| `noiseXu` | Xu (2019) fractional differencing | [-2, 7] stable; break-free with `Phi = 1` | Base MATLAB |

## Upgrading from 1.x

- `estimateIRASA` now performs IRASA. Versions up to 1.2.0 labelled the resampled copies with their true sampling rates, so peaks never moved and nothing was removed. On signals without peaks the exponent is unchanged; on signals with peaks it now differs, and is correct.
- `sigmaEst` is smaller than before, because the valid range now stops below `fs/(2*max(Hset))`.
- Bands that reach above `fs/(2*max(Hset))` now raise an error.
- To reproduce a 1.x result exactly, pass `Resampling="legacy"`. It warns on every call.
- `shapeNoise` builds a different surrogate. It finds the longest length L at which each channel has a segment whose ends meet (the wrap-around step and slopes are no larger than the largest step expected among the segment's own steps; the 1% energy budget and the fixed tolerance of the v2.0.0 drafts are options), puts uniformly random phases on that segment's amplitude spectrum, and redraws until the exponent is within 0.3 of the target. Before, 1-27% of surrogates at alpha 3-6 came out more than 0.5 shallower than the signal they were built from, because a Fourier surrogate of a record whose ends do not meet spreads the end jump across all frequencies.
- The output has L rows, not `numel(resid)`. Add it to something of length L.
- A second output, `info`, gives the segment (`i0`, `i1`, `L`, `Lfrac`), the target and achieved exponent, and the draw count. A channel with no accepted draw comes back as NaN with a warning.
- Steep records may have no segment whose ends meet; `shapeNoise` then errors rather than return something else. On exact references with no low-frequency break, at most 1.3% of series fail for alpha -0.75 to 2.5, 25-51% from 2.95 to 4, and 96-100% from 4.95 (see RANGE in the `shapeNoise` help and `validation/rangeEvidence_v002.md`).
- `alpha` is now the target for the selection step. Pass `[]` to use the segment's own exponent.
- `shapeNoise` searches segments of at least `MinLength` = 120 samples, so a record shorter than 120 samples is an error (`shapeNoise:short`). Pass `MinLength` lower only if you check the exponent yourself.
- `shapeNoise` warns (`shapeNoise:defaultAlphaBand`) when `AlphaBand` is not given: the default band [1 20] assumes movement data in Hz. Pass `AlphaBand="kinematic"` for that band without the warning, or `[fLow fHigh]` in the units of `fs`.
- A residual with a run of exact zeros at either end (for example an uncovered template tail) is an error. Pass only the part the template covered.
- `PhaseSource` now applies only to `Method="legacy"`, which reproduces 1.x bit for bit and warns on every call.

## Installation

**Toolbox file (recommended):** install "Fractal Noise Toolbox" from MATLAB File Exchange (184251) through the Add-On Explorer, or download the `.mltbx` from the repository's GitHub Releases and double-click it. This also puts `external/fftnoise/` on the path.

**Source folder:** from the folder you cloned the repository into:
```matlab
addpath(fullfile(pwd, "toolbox"));
addpath(fullfile(pwd, "toolbox", "external", "fftnoise"));   % only if you need fftnoise
```

Requires MATLAB R2025a or later for the documentation (plain-text Live Code format). The functions may run on older releases; untested.

## Examples and documentation

- `toolbox/doc/GettingStarted.m`: a short tour of every function.
- `toolbox/doc/WhyHomebrewIRASA.m`: what `estimateIRASA` does, the v2.0.0 fix, and what has been checked.
- `toolbox/examples/CompareFftnoiseVsShapedXu.m`: why Fourier phase randomisation of a whole record (`fftnoise`) fails to keep a steep spectrum, and how `shapeNoise` avoids it.
- `toolbox/examples/CompareNoiseGenerators.m`: `noiseXu` beside the built-in generators.
- `toolbox/examples/KneeCheckWorkedExample.m`: `checkSpectralKnee` against synthetic spectra with a known knee.

## Testing

From the folder you cloned the repository into:
```matlab
results = runtests(fullfile(pwd, "tests"));
```

Tests cover known-alpha recovery, peak removal, surrogate fidelity at steep alpha, fail-loud error paths, optional-dependency handling, a check that the toolbox works with only its own folder on the path, and a run of every README code block, example and documentation page. `validation/runMutationChecks.m` checks that the tests catch deliberate one-line bugs. Source-equivalence tests against the research code the functions came from need that code's folder in the environment variable `FRACTALNOISE_LEGACY_SRC`; without it they are reported as filtered, and with `FRACTALNOISE_REQUIRE_LEGACY=1` a missing folder is an error. CI runs the suite on every push to `main` and on pull requests.

## Validation

`validation/` holds the scripts behind the accuracy statements above: exact reference series (k-fold integrated Davies-Harte fGn, alpha up to 7, analytic spectra) and an estimator calibration in MATLAB, Python and R. Large run outputs are not in git.

## Credits

The bundled example series is the Hadley Centre Central England Temperature (HadCET) dataset. Contains public sector information licensed under the Open Government Licence v3.0. Source: Met Office Hadley Centre, www.metoffice.gov.uk/hadobs. Cite Parker, Legg & Folland (1992), Int. J. Climatol. 12, 317-342.

- `noiseXu`: Xu, C. (2019). An easy algorithm to generate colored noise sequences. *The Astronomical Journal*, 157(3), 127.
- `estimateWhittleAlpha`: Robinson, P. M. (1995). Gaussian semiparametric estimation of long range dependence. *Annals of Statistics*, 23(5), 1630-1661; Hurvich, C. M., & Chen, W. W. (2000). An efficient taper for potentially overdifferenced long-memory time series. *Journal of Time Series Analysis*, 21(2), 155-180; Velasco, C. (1999). Gaussian semiparametric estimation of non-stationary time series. *Journal of Time Series Analysis*, 20(1), 87-127; Hurvich, C., Lang, G., & Soulier, P. (2005). Estimation of long memory in the presence of a smooth nonparametric trend. *Journal of the American Statistical Association*, 100(471), 853-871; DFA order limit from Kiyono, K. (2015), *Physical Review E*, 92, 042925.
- `estimateIRASA`: Wen, H., & Liu, Z. (2016). Separating fractal and oscillatory components in the power spectrum of neurophysiological signal. *Brain Topography*, 29(1), 13-26.
- `checkSpectralKnee`: Preston, M., Smith, S., & Voytek, B. (2026). *Nature Human Behaviour*.
- `external/fftnoise/fftnoise.m`: (c) Aslak Grinsted 2011, BSD 2-clause, bundled unmodified.

The functions were first written for a kinematics project ([velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation)); each docstring names its source file. The toolbox has no dependency on it.

## Getting help and contributing

Questions and bug reports: open an issue at <https://github.com/dagmarfraser/fractal-noise/issues> (include your MATLAB release, `ver` output and a minimal seeded example). Pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md). To cite the toolbox, use [CITATION.cff](CITATION.cff) (GitHub's "Cite this repository" button).

## License

MIT (see `license.txt`), except `external/fftnoise/`, which keeps its BSD 2-clause license.
