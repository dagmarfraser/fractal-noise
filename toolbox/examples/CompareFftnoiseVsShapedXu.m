%[text] # Comparing fftnoise and End-Matched Surrogates
%[text] Phase-randomised surrogates ([fftnoise](https://www.mathworks.com/matlabcentral/fileexchange/32104-fftnoise-generate-noise-with-a-given-power-spectrum), Grinsted, 2011) preserve a signal's amplitude spectrum exactly -- in theory. In practice the FFT treats the record as periodic: when a steep series' ends do not meet, the end jump is a step, and its f^-2 leakage dominates the amplitude spectrum above the lowest bins. Random phases turn that leakage into real f^-2 noise, and the surrogate's *recovered* alpha collapses toward 2-3 once the target exceeds roughly alpha=3 (Schreiber \& Schmitz, 2000).
%[text] `shapeNoise` (default, end-matched) removes the cause rather than the symptom. It randomises the phases of the longest segment whose ends meet (wrap-around step and slope no larger than the segment's own average steps), then keeps a draw only if its exponent is within 0.3 of the segment's. The phases are just as random as fftnoise's; the difference is the segment.
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
%[text] For each target alpha: generate a genuine Xu-fractional reference signal via `noiseXu`, then build two surrogates from it -- a `fftnoise` phase-randomised copy of the whole record, and an end-matched `shapeNoise` copy -- and estimate the recovered alpha of *each surrogate* with all four methods. The end-matched surrogate is shorter than the reference (`info.L`); a draw that is never accepted comes back NaN and is counted, not scored.
alphaFftHB = nan(numel(targetAlphas), nReps);
alphaFftFT = nan(numel(targetAlphas), nReps);
alphaFftPM = nan(numel(targetAlphas), nReps);
alphaXuHB  = nan(numel(targetAlphas), nReps);
alphaXuFT  = nan(numel(targetAlphas), nReps);
alphaXuPM  = nan(numel(targetAlphas), nReps);
alphaFftW  = nan(numel(targetAlphas), nReps);
alphaXuW   = nan(numel(targetAlphas), nReps);
alphaRefW  = nan(numel(targetAlphas), nReps);
Lfrac      = nan(numel(targetAlphas), nReps);  nNoAccept = 0;
nWhittleFail = 0;  lastWhittleMsg = "";   % estimateWhittleAlpha errors loudly; failures are counted and reported below

for ai = 1:numel(targetAlphas)
    alphaTarget = targetAlphas(ai);
    for r = 1:nReps
        rng(1000 * ai + r);
        reference = noiseXu(N, alphaTarget, 1.0, fs);
        try, alphaRefW(ai, r) = estimateWhittleAlpha(reference, fs, fLow, fHigh); catch err, nWhittleFail = nWhittleFail + 1; lastWhittleMsg = string(err.message); end

        fftnoiseSurrogate = fftnoise(fft(reference), 1);
        [shapedXuSurrogate, info] = shapeNoise(reference, fs, []);   % end-matched (variable name kept from 1.x)
        Lfrac(ai, r) = info.Lfrac;

        [aHB, ~, ~, ~, ftInfo, pmInfo] = estimateIRASA(fftnoiseSurrogate, fs, fLow, fHigh, ...
            CompareFieldTrip=true, ComparePMTM=true);
        alphaFftHB(ai, r) = aHB;
        alphaFftFT(ai, r) = ftInfo.alphaFT;
        alphaFftPM(ai, r) = pmInfo.alphaPM;
        try, alphaFftW(ai, r) = estimateWhittleAlpha(fftnoiseSurrogate, fs, fLow, fHigh); catch err, nWhittleFail = nWhittleFail + 1; lastWhittleMsg = string(err.message); end

        if ~info.accepted, nNoAccept = nNoAccept + 1;  continue, end
        [aHB, ~, ~, ~, ftInfo, pmInfo] = estimateIRASA(shapedXuSurrogate, fs, fLow, fHigh, ...
            CompareFieldTrip=true, ComparePMTM=true);
        alphaXuHB(ai, r) = aHB;
        alphaXuFT(ai, r) = ftInfo.alphaFT;
        alphaXuPM(ai, r) = pmInfo.alphaPM;
        try, alphaXuW(ai, r) = estimateWhittleAlpha(shapedXuSurrogate, fs, fLow, fHigh); catch err, nWhittleFail = nWhittleFail + 1; lastWhittleMsg = string(err.message); end
    end
end
%[text] ## Result
%[text] Means over the repetitions. The first column is the reference series itself, read by W. Judge each surrogate against that column, not against the nominal target. Up to about alpha=3 both surrogates follow the reference. Above it the `fftnoise` surrogate stops following and stays between 2 and 3 on all four estimators. The end-matched surrogate follows the reference on W at every target (its own target is its segment's exponent, which can differ from the whole reference's by a little). HB, FT and PM are whole-record spectral fits; on the end-matched surrogate they also stay close (on the M5, R2026a, within 0.2 at alpha 5), because its ends meet, but on other whole records they can scatter by several tenths once alpha passes 3, so read W there.
fprintf("%-6s | %-6s | %-33s | %-33s | %s\n", "target", "ref W", "fftnoise (HB / FT / PM / W)", "end-matched (HB / FT / PM / W)", "L/N");
fprintf("%s\n", repmat('-', 1, 94));
mu = @(A, ai) mean(A(ai, :), 'omitnan');
for ai = 1:numel(targetAlphas)
    fprintf("%-6.1f | %6.2f | %6.2f / %6.2f / %6.2f / %6.2f | %6.2f / %6.2f / %6.2f / %6.2f | %4.2f\n", targetAlphas(ai), mu(alphaRefW, ai), ...
        mu(alphaFftHB, ai), mu(alphaFftFT, ai), mu(alphaFftPM, ai), mu(alphaFftW, ai), ...
        mu(alphaXuHB, ai), mu(alphaXuFT, ai), mu(alphaXuPM, ai), mu(alphaXuW, ai), mu(Lfrac, ai));
end
if nNoAccept > 0
    fprintf("End-matched: %d of %d references gave no accepted draw (NaN); left out of its columns.\n", nNoAccept, numel(Lfrac));
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
plot(targetAlphas, mean(alphaXuHB, 2, 'omitnan'), 'o-', 'DisplayName', 'end-matched (HB)', 'LineWidth', 1.5, 'Color', [0.85 0.33 0.10]);
plot(targetAlphas, mean(alphaFftW, 2, 'omitnan'), 'd-.', 'DisplayName', 'fftnoise (W)', 'LineWidth', 1.5);
plot(targetAlphas, mean(alphaXuW, 2, 'omitnan'), 'd-', 'DisplayName', 'end-matched (W)', 'LineWidth', 1.5);
plot(targetAlphas, targetAlphas, 'k--', 'DisplayName', 'target = recovered', 'LineWidth', 1);
xlabel('Target alpha (generator)');
ylabel('Recovered alpha (estimator)');
legend('Location', 'northwest');
title('fftnoise ceiling vs end-matched fidelity, four estimators');
grid on;
hold off;
%[text] ## Takeaways
%[text] - Above alpha~3, fftnoise is not a reliable surrogate for red/black biological noise regardless of which estimator checks it.
%[text] - `shapeNoise` (end-matched) is the substitute in that regime. The 1.x construction (`Method="legacy"`, *shaped_xu*: the whole record's magnitude on a `noiseXu` phase carrier) held the median but left about one surrogate in five more than 0.5 shallow between alpha 3.5 and 4.5; see the `shapeNoise` help.
%[text] - Checking with `CompareFieldTrip` alone would not have been enough here, since FT shares HB's resampling-median mechanism; `ComparePMTM` is the estimator that is actually independent of that shared step.
%[text] - The cause is the same in every case: an FFT treats the record as periodic, and a steep series whose ends do not meet leaks as f^-2. Differencing before the fit (W, `estimateIRASA`'s `Prewhiten`) removes it from an estimate; choosing a segment whose ends meet removes it from a surrogate.
%[text] ## References
%[text] Schreiber, T., \& Schmitz, A. (2000). Surrogate time series. *Physica D*, 142(3-4), 346-382.
%[text] Nakamura, T., Small, M., \& Hirata, Y. (2006). Testing for nonlinearity in irregular fluctuations with long-term trends. *Physical Review E*, 74, 026205.
%[text] Xu, C. (2019). An Easy Algorithm to Generate Colored Noise Sequences. *The Astronomical Journal*, 157:127.
%[text] Wen, H., \& Liu, Z. (2016). Separating Fractal and Oscillatory Components in the Power Spectrum of Neurophysiological Signal. *Brain Topography*, 29(1), 13-26.
%[text] Grinsted, A. (2011). fftnoise: Generate noise with a given power spectrum. MATLAB Central File Exchange.
%[text] Brookshire, G. (2022). Putative rhythms in attentional switching can be explained by aperiodic temporal structure. *Nature Human Behaviour*. https://doi.org/10.1038/s41562-022-01364-0

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
