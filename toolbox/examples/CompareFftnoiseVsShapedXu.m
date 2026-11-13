%[text] # Comparing fftnoise and shaped_xu Surrogates
%[text] Phase-randomised surrogates ([fftnoise](https://www.mathworks.com/matlabcentral/fileexchange/32104-fftnoise-generate-noise-with-a-given-power-spectrum), Grinsted, 2011) preserve a signal's amplitude spectrum exactly, so they inherit its autocorrelation structure -- in theory. In practice, phase randomisation destroys the long-range temporal correlations that sustain a steep 1/f^alpha slope, and the surrogate's *recovered* alpha collapses toward a much lower value once the target exceeds roughly alpha=3.
%[text] `shapeNoise` with `PhaseSource="xu"` (the *shaped_xu* construction) fixes this by imposing the empirical magnitude spectrum onto a genuinely fractional phase carrier (Xu, 2019) instead of a randomised one, so the phase structure -- not just the magnitude -- carries the intended alpha.
%[text] ## Why check with more than one estimator
%[text] A single spectral-exponent estimator agreeing with itself is not independent confirmation. This project checks with three deliberately different mechanisms:
%[text] - **HB** -- `estimateIRASA`'s own homebrew IRASA (resampling-median on pmtm PSDs, Wen \& Liu, 2016)
%[text] - **FT** -- FieldTrip's independent IRASA implementation (`CompareFieldTrip=true`), a different PSD estimator (Hanning-windowed FFT) but the *same* resampling-median mechanism as HB
%[text] - **PM** -- a direct pmtm log-log slope fit with *no* resampling step at all (`ComparePMTM=true`), the one genuinely independent check, since it shares neither PSD estimator nor resampling mechanism with the other two where it matters
%[text] Brookshire (2022, *Nature Human Behaviour*) makes the general case that naive single-method aperiodic/spectral characterisation is exactly where false confidence creeps in. That risk is highest here: this comparison runs well outside the range MATLAB's own `dsp.ColoredNoise` supports (InverseFrequencyPower in [-2, 2]) -- the whole point of `noiseXu` existing in the first place.
%[text] Before v2.0.0 the first of these was not a separate check: `estimateIRASA` did not relocate peaks, so HB was in effect PM (see `doc/WhyHomebrewIRASA.m`). From v2.0.0 the three are distinct.
%[text] All three are whole-record spectral fits and leak at steep alpha. A fourth column, **W**, is `estimateWhittleAlpha`, which differences the series before fitting and is the one to read for the exponent itself.
%[text] ## Setup
addpath(fullfile(fileparts(mfilename("fullpath")), "..", "external", "fftnoise"));

fs = 133;
N = 3000;
fLow = 1.0;
fHigh = min(20, fs/2 - 1);
targetAlphas = [2.0 3.0 4.0 5.0];
nReps = 10;
%[text] ## Generate ground truth, then both surrogate types
%[text] For each target alpha: generate a genuine Xu-fractional reference signal via `noiseXu`, then build two surrogates from it -- a `fftnoise` phase-randomised copy, and a `shapeNoise(...,PhaseSource="xu")` shaped_xu copy -- and estimate the recovered alpha of *each surrogate* with all four methods.
alphaFftHB = nan(numel(targetAlphas), nReps);
alphaFftFT = nan(numel(targetAlphas), nReps);
alphaFftPM = nan(numel(targetAlphas), nReps);
alphaXuHB  = nan(numel(targetAlphas), nReps);
alphaXuFT  = nan(numel(targetAlphas), nReps);
alphaXuPM  = nan(numel(targetAlphas), nReps);
alphaFftW  = nan(numel(targetAlphas), nReps);
alphaXuW   = nan(numel(targetAlphas), nReps);
alphaRefW  = nan(numel(targetAlphas), nReps);
nWhittleFail = 0;  lastWhittleMsg = "";   % estimateWhittleAlpha errors loudly; failures are counted and reported below

for ai = 1:numel(targetAlphas)
    alphaTarget = targetAlphas(ai);
    for r = 1:nReps
        rng(1000 * ai + r);
        reference = noiseXu(N, alphaTarget, 1.0, fs);
        try, alphaRefW(ai, r) = estimateWhittleAlpha(reference, fs, fLow, fHigh); catch err, nWhittleFail = nWhittleFail + 1; lastWhittleMsg = string(err.message); end

        fftnoiseSurrogate = fftnoise(fft(reference), 1);
        shapedXuSurrogate = shapeNoise(reference, fs, alphaTarget, PhaseSource="xu");

        [aHB, ~, ~, ~, ftInfo, pmInfo] = estimateIRASA(fftnoiseSurrogate, fs, fLow, fHigh, ...
            CompareFieldTrip=true, ComparePMTM=true);
        alphaFftHB(ai, r) = aHB;
        alphaFftFT(ai, r) = ftInfo.alphaFT;
        alphaFftPM(ai, r) = pmInfo.alphaPM;
        try, alphaFftW(ai, r) = estimateWhittleAlpha(fftnoiseSurrogate, fs, fLow, fHigh); catch err, nWhittleFail = nWhittleFail + 1; lastWhittleMsg = string(err.message); end

        [aHB, ~, ~, ~, ftInfo, pmInfo] = estimateIRASA(shapedXuSurrogate, fs, fLow, fHigh, ...
            CompareFieldTrip=true, ComparePMTM=true);
        alphaXuHB(ai, r) = aHB;
        alphaXuFT(ai, r) = ftInfo.alphaFT;
        alphaXuPM(ai, r) = pmInfo.alphaPM;
        try, alphaXuW(ai, r) = estimateWhittleAlpha(shapedXuSurrogate, fs, fLow, fHigh); catch err, nWhittleFail = nWhittleFail + 1; lastWhittleMsg = string(err.message); end
    end
end
%[text] ## Result
%[text] Means over the repetitions. The first column is the reference series itself, read by W. Judge each surrogate against that column, not against the nominal target. Up to about alpha=3 both surrogates follow the reference. Above it the `fftnoise` surrogate stops following and stays between 2 and 3 on all four estimators. The shaped_xu surrogate follows the reference on W at every target. HB, FT and PM are whole-record spectral fits and scatter around it by several tenths once alpha passes 3, so read W there.
fprintf("%-6s | %-6s | %-33s | %-33s\n", "target", "ref W", "fftnoise (HB / FT / PM / W)", "shaped_xu (HB / FT / PM / W)");
fprintf("%s\n", repmat('-', 1, 87));
mu = @(A, ai) mean(A(ai, :), 'omitnan');
for ai = 1:numel(targetAlphas)
    fprintf("%-6.1f | %6.2f | %6.2f / %6.2f / %6.2f / %6.2f | %6.2f / %6.2f / %6.2f / %6.2f\n", targetAlphas(ai), mu(alphaRefW, ai), ...
        mu(alphaFftHB, ai), mu(alphaFftFT, ai), mu(alphaFftPM, ai), mu(alphaFftW, ai), ...
        mu(alphaXuHB, ai), mu(alphaXuFT, ai), mu(alphaXuPM, ai), mu(alphaXuW, ai));
end
if nWhittleFail > 0
    fprintf("estimateWhittleAlpha failed on %d of %d series; those are left out of the W means. Last message: %s\n", ...
        nWhittleFail, 3 * numel(alphaFftW), lastWhittleMsg);
end
%[text] ## Visual comparison
figure;
hold on;
plot(targetAlphas, mean(alphaFftHB, 2, 'omitnan'), 'o-', 'DisplayName', 'fftnoise (HB)', 'LineWidth', 1.5);
plot(targetAlphas, mean(alphaFftFT, 2, 'omitnan'), 's--', 'DisplayName', 'fftnoise (FT)', 'LineWidth', 1.5);
plot(targetAlphas, mean(alphaFftPM, 2, 'omitnan'), '^:', 'DisplayName', 'fftnoise (PM)', 'LineWidth', 1.5);
plot(targetAlphas, mean(alphaXuHB, 2, 'omitnan'), 'o-', 'DisplayName', 'shaped\_xu (HB)', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
plot(targetAlphas, mean(alphaFftW, 2, 'omitnan'), 'd-.', 'DisplayName', 'fftnoise (W)', 'LineWidth', 1.5);
plot(targetAlphas, mean(alphaXuW, 2, 'omitnan'), 'd-', 'DisplayName', 'shaped\_xu (W)', 'LineWidth', 1.5);
plot(targetAlphas, targetAlphas, 'k--', 'DisplayName', 'target = recovered', 'LineWidth', 1);
xlabel('Target alpha (generator)');
ylabel('Recovered alpha (estimator)');
legend('Location', 'northwest');
title('fftnoise ceiling vs shaped\_xu fidelity, four estimators');
grid on;
hold off;
%[text] ## Takeaways
%[text] - Above alpha~3, fftnoise is not a reliable surrogate for red/black biological noise regardless of which estimator checks it.
%[text] - `shapeNoise` with `PhaseSource="xu"` (shaped_xu) is the substitute in that regime. Its default, `Prewhiten="auto"`, swaps magnitude and phase on the differenced series. Without it (`Prewhiten=0`, the behaviour before v2.0.0) about one surrogate in five came out more than 0.5 shallow between alpha 3.5 and 4.5; see the `shapeNoise` help.
%[text] - Checking with `CompareFieldTrip` alone would not have been enough here, since FT shares HB's resampling-median mechanism; `ComparePMTM` is the estimator that is actually independent of that shared step.
%[text] - The cause is the same in every case: an FFT treats the record as periodic, and a steep series whose ends do not meet leaks as f^-2. Differencing first (W, `Prewhiten`) removes it.
%[text] ## References
%[text] Xu, C. (2019). An Easy Algorithm to Generate Colored Noise Sequences. *The Astronomical Journal*, 157:127.
%[text] Wen, H., \& Liu, Z. (2016). Separating Fractal and Oscillatory Components in the Power Spectrum of Neurophysiological Signal. *Brain Topography*, 29(1), 13-26.
%[text] Grinsted, A. (2011). fftnoise: Generate noise with a given power spectrum. MATLAB Central File Exchange.
%[text] Brookshire, G. (2022). Putative rhythms in attentional switching can be explained by aperiodic temporal structure. *Nature Human Behaviour*. https://doi.org/10.1038/s41562-022-01364-0

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
