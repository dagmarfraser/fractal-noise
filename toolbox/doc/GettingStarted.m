%[text] # Getting Started with `fractalnoise`
%[text] A short tour: generate, estimate, shape, remove a periodic part, and diagnose 1/f^alpha noise. Run the sections in order.
%[text] ## Installation
%[text] Add the toolbox folder to your path. That's it -- no compiled components, no install script.
%[text] ```matlabCodeExample
%[text] addpath('/path/to/fractalnoise/toolbox');
%[text] ```
%[text] Everything below runs with nothing else on path. Two of the functions (`shapeNoise`, `estimateIRASA`) have **optional** extra features that need extra MathWorks toolboxes or FieldTrip -- covered at the end of this document -- but none of that is required for normal use.
%[text] ## 1. Generate: `noiseXu`
%[text] Generates 1/f^alpha noise via Xu's (2019) fractional-differencing method. Four arguments: length, spectral exponent alpha, target amplitude (standard deviation), sampling rate.
rng(1);
n = noiseXu(2000, 2.0, 1.5, 60);
fprintf("Generated %d samples, std=%.3f (target 1.5)\n", numel(n), std(n));
%%
%[text] alpha=0 is white noise, alpha=1 is pink, alpha=2 is Brownian/red, alpha=3 is black. Non-integer and negative values are also valid -- see `noiseXu`'s docstring for the full verified range.
%[text] ## 2. Estimate: `estimateIRASA`
%[text] Recovers alpha (and a noise-magnitude estimate) from a time series via IRASA, which sets rhythmic peaks aside and fits the scale-free background. With the default settings the band must satisfy `fHigh < fs/3.8`. Takes the signal, sampling rate, and a frequency band to fit over.
[alphaEst, sigmaEst] = estimateIRASA(n, 60, 1.0, 15.0);
fprintf("Recovered alpha=%.3f (target 2.0)\n", alphaEst);
fprintf("Noise magnitude sigma=%.3f\n", sigmaEst);
%%
%[text] Note `sigmaEst` integrates the fractal PSD from the lowest frequency up to `fs/(2*max(Hset))` (not restricted to the `[fLow, fHigh]` band used for the alpha fit), so it won't generally equal the generator's amplitude parameter exactly -- it's a spectral-magnitude estimate, not a reconstruction of the original time-domain standard deviation. Use `std()` directly on your signal if you want the latter.
%[text] Ask for a second opinion from FieldTrip (if installed) with `CompareFieldTrip=true`:
[~, ~, ~, ~, ftInfo] = estimateIRASA(n, 60, 1.0, 15.0, CompareFieldTrip=true);
if ftInfo.available
    fprintf("FieldTrip: alpha=%.3f (gap from estimateIRASA: %.3f)\n", ftInfo.alphaFT, ftInfo.gapHBFT);
else
    fprintf("FieldTrip comparison unavailable: %s\n", ftInfo.message);
end
%%
%[text] See `doc/WhyHomebrewIRASA.m` for why this function has its own IRASA implementation, and for the bug fixed in v2.0.0 (earlier versions did not remove peaks).
%[text] ## Steep spectra: `estimateWhittleAlpha`
%[text] Above about alpha=3, whole-record spectra leak and IRASA can read wrong. `estimateWhittleAlpha` differences the series first (it reports the order it used), then fits a local Whittle likelihood over the band.
[alphaW, kUsed] = estimateWhittleAlpha(n, 60, 1.0, 15.0);
fprintf("Whittle alpha=%.3f (differencing order %d)\n", alphaW, kUsed);
%%
%[text] To keep IRASA's peak removal at steep alpha, pass that order on: `estimateIRASA(n, 60, 1.0, 15.0, Prewhiten=kUsed)`.
%[text] ## 3. Shape: `shapeNoise`
%[text] Builds a random-phase surrogate of an **empirical** residual -- data with the residual's spectrum and derivative noise, but not its realisation. A Fourier surrogate of a whole record whose ends do not meet spreads the end jump over every frequency, so `shapeNoise` first finds the longest segment whose ends meet, randomises the phases of that segment, and keeps a draw only if its exponent is within 0.3 of the target (`[]`: the segment's own). The surrogate therefore has `info.L` samples, not `numel(myResidual)`.
myResidual = cumsum(randn(500, 1));
myResidual = myResidual - mean(myResidual);

[surrogate, info] = shapeNoise(myResidual, 100, [], AlphaBand=[1 20]);
seg = myResidual(info.i0:info.i1);
fprintf("Segment %d-%d (L=%d, %.0f%% of the record); std %.3f vs %.3f; alpha %.2f vs target %.2f\n", ...
    info.i0, info.i1, info.L, 100 * info.Lfrac, std(surrogate), std(seg), info.alphaTarget + info.dAlpha, info.alphaTarget);
%%
%[text] `resid` can have one column per channel (e.g. x and y of a movement); all channels share L. For template-subtracted data, pass only the part the template covered: a run of exact zeros at either end is an error. `Method="legacy"` gives the 1.x construction (with its `PhaseSource` carriers) for reproducing earlier results, and warns.
%[text] ## Periodic records: `harmonicContent`, `harmonicTemplate` and the wrappers
%[text] A movement is mostly its periodic part; the noise is what is left after removing it. `harmonicContent` shows how many harmonics of f0 stand above the noise beside them, `harmonicTemplate` removes that many, and `shapeNoise` then works on the residual. `shapeNoiseWrapper` does both steps in one call and requires K (there is no default).
t = (0:1/100:30)';
rng(4);
pos = 2*sin(2*pi*1*t) + 0.5*sin(2*pi*2*t + 0.7) + noiseXu(numel(t), 3, 0.05, 100, Phi=1);
[K, hc] = harmonicContent(pos, 100, 1);
res = harmonicTemplate(pos, 100, 1, K);
fprintf("K = %d (harmonic ratios %s); residual %d samples\n", K, mat2str(round(hc.ratio(1:3), 1)), numel(res));
%%
%[text] See the README section "Choosing K for your data" before relying on one count.
%[text] ## Many surrogates: `summariseShapeNoiseBatch`
%[text] One call gives one surrogate per channel, so acceptance only means something across a batch. Collect the `info` outputs and summarise them; the function warns only when acceptance is significantly below 95%.
infos = cell(1, 50);
for i = 1:50
    [~, infos{i}] = shapeNoise(res, 100, [], AlphaBand=[1 20]);
end
S = summariseShapeNoiseBatch(infos);
fprintf("%d of %d accepted; median tries %g; decision: %s\n", S.nAccepted, S.nChannels, S.triesMedian, S.decision);
%%
%[text] ## 4. Diagnose: `checkSpectralKnee`
%[text] Tests whether a residual's spectrum is better explained by a Lorentzian (kneed) model than a plain power law -- relevant because IRASA-based alpha estimates (like `estimateIRASA` above) assume a scale-free, non-kneed spectrum.
result = checkSpectralKnee(myResidual, 100);
fprintf("Verdict: %s (deltaAIC=%.2f, knee at %.3f Hz, in fit band: %d)\n", ...
    result.verdict, result.deltaAIC, result.kneeHz, result.kneeInBand);
%%
%[text] See `examples/KneeCheckWorkedExample.m` for this diagnostic against synthetic spectra with a known knee: below the fitting band, inside it, and absent.
%[text] ## Optional toolbox dependencies
%[text] Two functions have optional features that degrade **visibly**, not silently, when the underlying toolbox isn't installed: requesting the feature without the dependency returns a clearly-flagged unavailable result and a warning, never a silent fallback or an outright crash.
%[text] - `shapeNoise(..., Method="legacy", PhaseSource="pinknoise")` needs MATLAB's built-in  `pinknoise` function (Signal Processing Toolbox).
%[text] - `shapeNoise(..., Method="legacy", PhaseSource="dsp")` needs `dsp.ColoredNoise` (DSP System Toolbox).
%[text] - `estimateIRASA(..., CompareFieldTrip=true)` needs FieldTrip ([https://www.fieldtriptoolbox.org](https://www.fieldtriptoolbox.org)) on the path, specifically  `ft_defaults` having been run. \
%[text] None of the core functions' default behaviour needs any of these.
%[text] ## Where to go next
%[text] - `doc/WhyHomebrewIRASA.m` -- what `estimateIRASA` does, the v2.0.0 fix, and what has been checked since
%[text] - `examples/KneeCheckWorkedExample.m` -- `checkSpectralKnee` against a known ground truth
%[text] - `examples/CompareNoiseGenerators.m` -- a side-by-side look at the different noise-generation and phase-shaping options
%[text] - `tests/` -- the test suite doubles as executable documentation of expected behaviour for every function \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
