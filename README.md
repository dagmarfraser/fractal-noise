# fractalnoise

A small, dependency-light MATLAB toolbox for generating, shaping and estimating
1/f^alpha ("fractal" or "colored") noise, with an emphasis on **steep spectra**
(alpha 2-7), where most generators and estimators were never checked.

It is useful wherever you need noise with a controlled spectral exponent:
simulation grids over alpha, null models and surrogate data for a measured
signal, or estimating the aperiodic exponent of a residual.

**Status:** development branch for the next release (in preparation). Last
release v1.1.0 (2026-07-15, MATLAB File Exchange 184251). See
[CHANGELOG.md](CHANGELOG.md). Source:
[github.com/dagmarfraser/fractal-noise](https://github.com/dagmarfraser/fractal-noise).

## What's in it

| Function | Purpose | Needs |
|---|---|---|
| `noiseXu` | Generate 1/f^alpha noise by fractional differencing (Xu 2019) | Base MATLAB |
| `shapeNoise` | Surrogate of a measured signal: its amplitude spectrum on a chosen phase carrier (default `noiseXu`) | Base MATLAB (two optional carriers need Signal Processing / DSP System Toolbox) |
| `estimateWhittleAlpha` | Spectral exponent by prewhitened, band-limited local Whittle; reliable at steep alpha | Base MATLAB |
| `estimateIRASA` | IRASA spectral exponent and noise magnitude (Wen & Liu 2016) | Signal Processing Toolbox (`pmtm`); optional FieldTrip cross-check |
| `checkSpectralKnee` | Lorentzian knee vs fixed power law (Preston, Smith & Voytek 2026) | Signal Processing Toolbox |

Optional features (FieldTrip, `dsp.ColoredNoise`, `pinknoise` carriers) degrade
visibly: if the dependency is missing you get a flagged "unavailable" result and
a warning, never a silent substitute.

## Which estimator to use

- **Steep spectra (alpha > ~3), or you don't know:** `estimateWhittleAlpha`. It
  differences the series k times (k chosen from a DFA-3 pre-estimate) so that
  low-frequency power cannot leak into the fitting band, then fits a local
  Whittle likelihood over your band. On exact reference series and on
  `noiseXu` output it recovers the in-band exponent within about +/-0.1 for
  alpha 2.5-6.5 at N = 1000-3000.
- **`estimateIRASA`** is kept for continuity and for signals with oscillatory
  peaks. Its spectrum step is a whole-record multitaper (`pmtm`), and at steep
  alpha that is vulnerable to red-noise leakage: on noise *without* a
  low-frequency break it can be off by several units, in a direction that
  depends on sampling rate and record length. On noise *with* a break (which
  includes `noiseXu` output and, in our tests, template-subtracted movement
  residuals) the bias is about +0.1 to +0.3. `Prewhiten=k` (difference k
  times, then recolour) removes the break-free bias (within about +/-0.2 to
  alpha 5.5); `PSD="welch"` alone only halves it. If in doubt, run both
  estimators, e.g. `estimateIRASA(x, fs, 1, 20, Prewhiten=k)` with `k` taken
  from `[~, k] = estimateWhittleAlpha(x, fs, 1, 20)`.

## What `noiseXu` generates (read this before using alpha > 2)

`noiseXu` implements Xu's generalized Gauss-Markov (GGM) recursion,
x = (1 - Phi*B)^(-alpha/2) w, with Phi = 0.99 by default. Two consequences:

- For alpha <= 2 this is a 1/f^alpha process down to a low-frequency break at
  f_b = (1 - Phi)*fs/(2*pi) (about 0.0016*fs). Below the break the spectrum
  flattens.
- For alpha > 2 it is Xu's *range-limited* approximation of a steeper power law
  (Xu 2019 discusses alpha up to about 4). The generator is numerically stable
  to alpha = 7, but the series is a power law with a break, not a pure
  1/f^alpha process. Score estimators against the GGM spectrum, not nominal
  alpha, and use exact references if you need a pure power law.
- `Phi = 1` removes the break: plain fractional integration, which in our
  tests reads the same as exact integrated fGn. Use `noiseXu(N, alpha, sd,
  fs, Phi=1)` to test an estimator on break-free data.

The break is deliberate (Xu chose it to reduce red-noise leakage), and it makes
spectral estimators look better on `noiseXu` output than on break-free data.
Keep "the generator is stable" and "an estimator can recover alpha" apart.

## Where `noiseXu` sits among MATLAB options

| | Method | Alpha range | Requires |
|---|---|---|---|
| `randn` | White noise | 0 | Base MATLAB |
| `pinknoise` | Pink noise | 1 | Signal Processing Toolbox |
| `dsp.ColoredNoise` | Kasdin (1995) filter recursion | [-2, 2] | DSP System Toolbox |
| `fractalcoef` | Kasdin (1995), raw filter coefficients | (0, 2) | Sensor Fusion and Tracking Toolbox |
| `noiseXu` | Xu (2019) GGM fractional differencing | [-2, 7] stable; pure power law to 2 (see above) | Base MATLAB |

## Installation

**Toolbox file (recommended):** double-click `release/Fractal Noise Toolbox.mltbx`,
or `matlab.addons.toolbox.installToolbox("release/Fractal Noise Toolbox.mltbx")`.
This also puts `external/fftnoise/` on the path.

**Source folder:**
```matlab
addpath('/path/to/fractalnoise/toolbox');
addpath('/path/to/fractalnoise/toolbox/external/fftnoise');   % only if you need fftnoise
```

Requires MATLAB R2025a or later for the documentation (plain-text Live Code
format). The functions themselves may run on older releases if copied
directly; untested.

## Quick example

```matlab
x = noiseXu(3000, 4.5, 1.0, 240);                 % alpha 4.5, unit SD, 240 Hz
[a, k] = estimateWhittleAlpha(x, 240, 1, 20)      % exponent over 1-20 Hz, and differencing order used
[aI, sI] = estimateIRASA(x, 240, 1, 20);          % IRASA for comparison
s = shapeNoise(myResidual, 240, a);               % surrogate of a measured residual
r = checkSpectralKnee(myResidual, 240);           % knee or no knee?
```

`examples/CompareFftnoiseVsShapedXu.m` shows why plain Fourier phase
randomisation (`fftnoise`) fails to keep a steep spectrum and how `shapeNoise`
avoids it. `doc/` has the function pages.

## Testing

```matlab
results = runtests('/path/to/fractalnoise/tests');
```

Tests cover known-alpha recovery, fail-loud error paths, optional-dependency
handling, and a check that the toolbox works with only its own folder on the
path. Source-equivalence tests against the research code the functions came
from skip themselves when that code is absent.

## Validation

`validation/` holds the scripts behind the claims above: exact reference series
(k-fold integrated Davies-Harte fGn, alpha up to 7, analytic spectra) and an
estimator calibration in MATLAB, Python and R. Large run outputs are not in git.

## Credits

- `noiseXu`: Xu, C. (2019). An easy algorithm to generate colored noise
  sequences. *The Astronomical Journal*, 157(3), 127.
- `estimateWhittleAlpha`: Robinson, P. M. (1995). Gaussian semiparametric
  estimation of long range dependence. *Annals of Statistics*, 23(5),
  1630-1661; DFA order limit from Kiyono, K. (2015), *Physical Review E*, 92,
  042925.
- `estimateIRASA`: Wen, H., & Liu, Z. (2016). Separating fractal and
  oscillatory components in the power spectrum of neurophysiological signal.
  *Brain Topography*, 29(1), 13-26.
- `checkSpectralKnee`: Preston, M., Smith, S., & Voytek, B. (2026). *Nature
  Human Behaviour*.
- `external/fftnoise/fftnoise.m`: (c) Aslak Grinsted 2011, BSD 2-clause,
  bundled unmodified.

The functions were first written for a kinematics project
([velocity-curvature-power-law-simulation](https://github.com/dagmarfraser/velocity-curvature-power-law-simulation));
each docstring names its source file. The toolbox has no dependency on it.

## License

MIT (see `license.txt`), except `external/fftnoise/`, which keeps its BSD
2-clause license.
