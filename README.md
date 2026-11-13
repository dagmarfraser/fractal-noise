# fractalnoise

**Make, match and measure 1/f^alpha noise in MATLAB, including the steep spectra (alpha 3-7) that movement and other behavioural signals have and that most tools were never checked on.**

Use it when you need to:

1. **Simulate noise with a chosen spectral exponent**, for example to sweep alpha in a simulation grid or to test an analysis pipeline against known noise.
   `x = noiseXu(3000, 4.5, 1, 240);`
2. **Build a null model of a signal you measured**: a surrogate with the same spectrum as your residual but random phase, which still holds its shape at steep alpha where ordinary Fourier phase randomisation does not.
   `s = shapeNoise(myResidual, 240, alpha);`
3. **Measure the spectral exponent of a signal**, with an estimator that stays accurate at steep alpha, and separate a scale-free background from rhythmic peaks.
   `alpha = estimateWhittleAlpha(x, 240, 1, 20);`
4. **Check whether a spectrum is really a power law** or bends at a knee.
   `r = checkSpectralKnee(myResidual, 240);`

Typical users: movement science and motor control, psychophysics, physiology, and anyone who needs coloured noise outside the white-to-brown range (alpha 0-2) that built-in MATLAB generators cover.

MIT licence. Base MATLAB for the generator, the surrogate and the Whittle estimator; Signal Processing Toolbox for the other two. Current release: see [CHANGELOG.md](CHANGELOG.md). Also on MATLAB File Exchange (184251).

> **v2.0.0 fixes two things.** `estimateIRASA` did not remove oscillatory peaks up to v1.2.0, so its exponent was a plain spectral slope. `shapeNoise` gave a shallow surrogate about one time in five at alpha 3.5-4.5. Both defaults are now correct, and both old behaviours stay available by name. See [Upgrading from 1.x](#upgrading-from-1x).

## Quick start

```matlab
x = noiseXu(3000, 4.5, 1.0, 240);                  % alpha 4.5, unit SD, 240 Hz
[a, k] = estimateWhittleAlpha(x, 240, 1, 20)       % exponent over 1-20 Hz, and the differencing order used
aI = estimateIRASA(x, 240, 1, 20, Prewhiten=k);    % IRASA on the same band, safe at steep alpha
s = shapeNoise(myResidual, 240, a);                % surrogate of a measured residual
r = checkSpectralKnee(myResidual, 240);            % knee or no knee?
```

## What's in it

| Function | What it does | Needs |
|---|---|---|
| `noiseXu` | Generates 1/f^alpha noise by fractional differencing (Xu 2019), alpha from -2 to 7 | Base MATLAB |
| `shapeNoise` | Surrogate of a measured signal: its amplitude spectrum on a chosen phase carrier (default `noiseXu`), swapped on the differenced series so steep spectra survive | Base MATLAB (two optional carriers need Signal Processing / DSP System Toolbox) |
| `estimateWhittleAlpha` | Spectral exponent over a band by prewhitened local Whittle; the one to use at steep alpha | Base MATLAB |
| `estimateIRASA` | IRASA (Wen & Liu 2016): spectral exponent and noise magnitude of the scale-free part, with oscillatory peaks removed | Signal Processing Toolbox; optional FieldTrip cross-check |
| `checkSpectralKnee` | Lorentzian knee against a fixed power law (Preston, Smith & Voytek 2026) | Signal Processing Toolbox |

Optional features (FieldTrip, `dsp.ColoredNoise`, `pinknoise` carriers) degrade visibly: a missing dependency gives a flagged "unavailable" result and a warning, never a silent substitute.

## Which estimator to use

| Your signal | Use |
|---|---|
| Steep (alpha above about 3), or unknown | `estimateWhittleAlpha` |
| Has rhythmic peaks you want ignored | `estimateIRASA`, with `Prewhiten=k` if it may be steep |
| Both | Run both and compare: `[~, k] = estimateWhittleAlpha(x, fs, 1, 20); estimateIRASA(x, fs, 1, 20, Prewhiten=k)` |

- **`estimateWhittleAlpha`** differences the series k times (k from a DFA-3 pre-estimate), so low-frequency power cannot leak into the fitting band, then fits a local Whittle likelihood over your band. On exact reference series and on `noiseXu` output it is within about +/-0.2 for alpha up to 6.5.
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
- `shapeNoise` now differences the signal and the carrier before swapping magnitude and phase (`Prewhiten="auto"`), then integrates back. Before, 10-25% of surrogates at alpha 3.5-5 came out more than 0.5 shallower than the signal they were built from; now 2% or fewer do. Output is unchanged below alpha 1. With `Prewhiten` above 0 the surrogate matches the magnitude spectrum of the differenced signal, not `abs(fft(x))` itself. `Prewhiten=0` reproduces 1.x exactly.

## Installation

**Toolbox file (recommended):** double-click `release/Fractal Noise Toolbox.mltbx`, or `matlab.addons.toolbox.installToolbox("release/Fractal Noise Toolbox.mltbx")`. This also puts `external/fftnoise/` on the path.

**Source folder:**
```matlab
addpath('/path/to/fractalnoise/toolbox');
addpath('/path/to/fractalnoise/toolbox/external/fftnoise');   % only if you need fftnoise
```

Requires MATLAB R2025a or later for the documentation (plain-text Live Code format). The functions may run on older releases; untested.

## Examples and documentation

- `examples/CompareFftnoiseVsShapedXu.m`: why Fourier phase randomisation (`fftnoise`) fails to keep a steep spectrum, and how `shapeNoise` avoids it.
- `examples/CompareNoiseGenerators.m`: `noiseXu` beside the built-in generators.
- `doc/`: getting started and function pages.

## Testing

```matlab
results = runtests('/path/to/fractalnoise/tests');
```

Tests cover known-alpha recovery, peak removal, surrogate fidelity at steep alpha, fail-loud error paths, optional-dependency handling, and a check that the toolbox works with only its own folder on the path. Source-equivalence tests against the research code the functions came from need that code's folder in the environment variable `FRACTALNOISE_LEGACY_SRC`; without it they are reported as filtered, and with `FRACTALNOISE_REQUIRE_LEGACY=1` a missing folder is an error. CI runs the suite on every push to `main` and on pull requests.

## Validation

`validation/` holds the scripts behind the accuracy statements above: exact reference series (k-fold integrated Davies-Harte fGn, alpha up to 7, analytic spectra) and an estimator calibration in MATLAB, Python and R. Large run outputs are not in git.

## Credits

- `noiseXu`: Xu, C. (2019). An easy algorithm to generate colored noise sequences. *The Astronomical Journal*, 157(3), 127.
- `estimateWhittleAlpha`: Robinson, P. M. (1995). Gaussian semiparametric estimation of long range dependence. *Annals of Statistics*, 23(5), 1630-1661; DFA order limit from Kiyono, K. (2015), *Physical Review E*, 92, 042925.
- `estimateIRASA`: Wen, H., & Liu, Z. (2016). Separating fractal and oscillatory components in the power spectrum of neurophysiological signal. *Brain Topography*, 29(1), 13-26.
- `checkSpectralKnee`: Preston, M., Smith, S., & Voytek, B. (2026). *Nature Human Behaviour*.
- `external/fftnoise/fftnoise.m`: (c) Aslak Grinsted 2011, BSD 2-clause, bundled unmodified.

The functions were first written for a kinematics project ([velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation)); each docstring names its source file. The toolbox has no dependency on it.

## Getting help and contributing

Questions and bug reports: open an issue at <https://github.com/dagmarfraser/fractal-noise/issues> (include your MATLAB release, `ver` output and a minimal seeded example). Pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md). To cite the toolbox, use [CITATION.cff](CITATION.cff) (GitHub's "Cite this repository" button).

## License

MIT (see `license.txt`), except `external/fftnoise/`, which keeps its BSD 2-clause license.
