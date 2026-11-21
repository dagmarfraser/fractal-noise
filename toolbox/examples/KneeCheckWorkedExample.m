%[text] # `checkSpectralKnee` Worked Example: Synthetic Ground Truth
%[text] `checkSpectralKnee` applies the Lorentzian knee diagnostic of Preston, Smith & Voytek (2026) to one residual. This example tests it against **synthetic data with an exactly known ground-truth knee**, so its verdicts and recovered numbers can be checked, in three cases: a knee below the fitting band, a knee inside it, and no knee.
%[text] ## Background: why check for a knee at all
%[text] IRASA (used by `estimateIRASA`) assumes its input is scale-free -- a pure power-law spectrum with no preferred timescale. Real electrophysiological and, potentially, kinematic signals can instead follow a **Lorentzian** form: flat below a "knee frequency" `f_knee`, decaying as a power law above it. If a knee falls inside the fitting band, IRASA's alpha estimate is biased by the flat-to-decaying transition it's straddling.
%[text] `checkSpectralKnee` answers two questions for a single residual: (1) is there statistical support for a knee at all (via AIC comparison against a plain fixed power-law model), and (2), if so, does that knee fall inside the band you're fitting alpha over?
%[text] ## Building synthetic data with an EXACT known knee
%[text] The real data can tell you a diagnostic's **qualitative** verdict was right, but not whether its **recovered numbers** (kneeHz, alpha) were close to a true underlying value -- the true value isn't observable in real data. This section builds surrogate residuals with an analytically exact Lorentzian spectrum, so the recovered numbers can be checked against a genuine ground truth.
%[text] The construction: compute the target Lorentzian magnitude spectrum `(kneeHz^2 + f^2)^(-alpha/4)` (square root of the Lorentzian PSD used by `checkSpectralKnee`'s own model) directly at every FFT bin frequency for a signal of length N at sampling rate fs, then hand that magnitude array to the toolbox's bundled `fftnoise` (Grinsted, 2011) to randomise its phase while preserving the exact magnitude -- precisely what `fftnoise` is for, just fed an analytic target spectrum instead of an empirical one. This example calls fftnoise directly, so it needs `external/fftnoise/` on the path in addition to `toolbox/` -- a plain `addpath('.../toolbox')` does not reach subfolders, so that's added explicitly below rather than assumed.
exampleDir = fileparts(mfilename("fullpath"));
addpath(fullfile(exampleDir, "..", "external", "fftnoise"));

fs = 133;
N = 4000;

lorentzianTemplate = @(alpha, kneeHz, sigma) buildLorentzianTemplate(N, fs, alpha, kneeHz, sigma);
%%
%[text] ## Case 1: sub-floor knee
%[text] True alpha=4.5, true knee=0.28 Hz (below the standard 1 Hz fit floor)
rng(1);
subFloorSignal = lorentzianTemplate(4.5, 0.28, 3.0);
resultSubFloor = checkSpectralKnee(subFloorSignal, fs, WideBandLow=0.05);
fprintf("Sub-floor knee: verdict=%-9s deltaAIC=%8.2f  kneeHz=%.4f (true 0.28)  alphaWide=%.3f (true 4.50)  kneeInBand=%d\n", ...
    resultSubFloor.verdict, resultSubFloor.deltaAIC, resultSubFloor.kneeHz, resultSubFloor.alphaWide, resultSubFloor.kneeInBand);
%%
%[text] Strong knee support, but `kneeInBand=false`: a knee exists below the fitting band and does not bias alpha. Because the truth is known, the recovered `kneeHz` can be checked directly against the value used to generate the signal. `alphaWide` is the fixed power-law fit over the wide band (the model the knee model is compared against), not the knee model's own exponent; with the knee this far below the band it stays near the true alpha.
%[text] ## Case 2: in-band knee
%[text] True alpha=4.5, true knee=3 Hz (genuinely inside the default \[1,fHigh\] fit band). This is the actual bias-risk scenario the diagnostic exists to catch.
rng(2);
inBandSignal = lorentzianTemplate(4.5, 3.0, 3.0);
resultInBand = checkSpectralKnee(inBandSignal, fs);
fprintf("In-band knee:   verdict=%-9s deltaAIC=%8.2f  kneeHz=%.4f (true 3.00)  alphaWide=%.3f (true 4.50)  kneeInBand=%d\n", ...
    resultInBand.verdict, resultInBand.deltaAIC, resultInBand.kneeHz, resultInBand.alphaWide, resultInBand.kneeInBand);
fprintf("  alphaStd (naive fixed-band fit, straddling the knee): %.3f -- biased away from true 4.50\n", resultInBand.alphaStd);
%%
%[text] Here `kneeInBand=true`, correctly flagging the risk. And critically, `alphaStd` -- the naive power-law fit over the standard band, which straddles the knee -- comes out well below the true alpha=4.5, and so does `alphaWide`, the fixed fit over the wide band: neither fixed power law can follow a knee inside its range. This is the case `checkSpectralKnee` exists to catch: not just "a knee exists" but "a knee is biasing your alpha estimate", shown here against a known ground truth. (The knee model's own exponent is not returned; use `kneeHz` and `kneeInBand` to decide whether a fixed-band alpha can be trusted.)
%[text] ## Case 3: no knee at all, for contrast
%[text] Break-free noise (`noiseXu` with `Phi=1`). The default `noiseXu` (`Phi=0.99`) is not knee-free: it has a real low-frequency break at about 0.0016\*fs (0.21 Hz here), and `checkSpectralKnee` detects it.
rng(3);
flatSignal = noiseXu(N, 3.8, 3.0, fs, Phi=1);
resultFlat = checkSpectralKnee(flatSignal, fs);
fprintf("No knee:        verdict=%-9s deltaAIC=%8.2f  kneeHz=%.4f\n", ...
    resultFlat.verdict, resultFlat.deltaAIC, resultFlat.kneeHz);
%[text] deltaAIC near +2 (the pure parameter-penalty signature), fitted knee collapsed to (near) zero.
%[text] ## Takeaway for using this function
%[text] - A "none" verdict means treat the fixed power-law alpha estimate as reliable for that residual.
%[text] - A "marginal" or "strong" verdict with `kneeInBand=false` means a knee exists but isn't biasing your alpha fit -- informative, not alarming (Case 1 above).
%[text] - A "marginal" or "strong" verdict with `kneeInBand=true` is the case that actually warrants caution -- Case 2 above quantifies exactly how much bias that can mean.
%[text] - Run this per trial, not on a single pooled PSD across many trials: pooling across trials with different structure can produce spurious or boundary-artefact knee fits. \
%[text] ## Cross-references
%[text] - `checkSpectralKnee.m` docstring
%[text] - Preston, M., Smith, S. & Voytek, B. (2026). Potential mechanisms and functional significance of aperiodic neural activity. Nature Human Behaviour.
%[text] - `external/fftnoise/fftnoise.m`, Grinsted (2011) -- the phase randomiser used to build the synthetic templates above \
function x = buildLorentzianTemplate(N, fs, alpha, kneeHz, sigma)
%BUILDLORENTZIANTEMPLATE Synthetic residual with an exact known Lorentzian spectrum.
%   Example-local helper, not a toolbox function -- scoped to this demo.
%   Builds the target magnitude spectrum (kneeHz^2+f^2)^(-alpha/4)
%   directly at each FFT bin frequency (MATLAB's standard bin ordering:
%   DC, ascending positive frequencies, then negative), then hands it to
%   the bundled fftnoise (Grinsted, 2011) to randomise phase while
%   preserving that exact magnitude.
freqIdx = (0:N - 1)';
freqHz = freqIdx * fs / N;
freqHz(freqHz > fs / 2) = freqHz(freqHz > fs / 2) - fs;
magSpectrum = (kneeHz ^ 2 + freqHz .^ 2) .^ (-alpha / 4);
x = fftnoise(magSpectrum);
x = x - mean(x);
x = x * (sigma / std(x));
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
