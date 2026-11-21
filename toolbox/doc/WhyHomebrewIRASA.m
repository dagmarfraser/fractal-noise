%[text] # `estimateIRASA`: Own Implementation, and What Changed in v2.0.0
%[text] `estimateIRASA` implements IRASA (Wen and Liu 2016) itself instead of wrapping FieldTrip's `ft_freqanalysis(method="irasa")`. This page says why, what was wrong with it up to v1.2.0, and what has and has not been checked since.
%[text] ## Why not just call FieldTrip
%[text] - **No dependency.** The function needs only Signal Processing Toolbox. FieldTrip stays an optional cross-check (`CompareFieldTrip=true`).
%[text] - **Choice of spectrum.** Multitaper (`PSD="pmtm"`, default) or Welch (`PSD="welch"`).
%[text] - **Steep spectra.** `Prewhiten=k` differences the series before IRASA and recolours afterwards, which removes the leakage that biases any whole-record spectrum once alpha passes about 3. FieldTrip has no equivalent. \
%[text] ## The bug in v1.2.0 and earlier
%[text] IRASA resamples the signal by a set of factors h and takes the spectrum of every copy **as if it were still sampled at fs**. A scale-free background keeps its shape under that relabelling, while a rhythmic peak at f0 moves to f0\*h and f0/h. The median over h then rejects the peak.
%[text] Up to v1.2.0 this function gave each copy its true rate (fs\*h, fs/h) and interpolated back to the original frequency axis. Every copy then estimated the same spectrum at the same frequencies. Peaks never moved, the median had nothing to reject, and the returned exponent was a plain spectral slope. The research code it was extracted from had the same form from its first version.
%[text] It went unnoticed because every test used peak-free noise, where correct IRASA, the old code and a direct slope fit all agree.
%[text] ## What that does to the earlier three-way validation
%[text] Earlier versions of this page reported a three-way check: this function (HB), FieldTrip (FT) and a direct `pmtm` slope with no resampling (PM), on `noiseXu` output for alpha from -2 to 7. HB and PM agreed within 0.005 up to alpha 4, and that was read as independent confirmation.
%[text] It was not independent. With the old resampling, HB **was** PM apart from interpolation, so the agreement was close to an identity. Two further cautions apply to those numbers:
%[text] - They were measured on default `noiseXu` output, which has a low-frequency break (see the `noiseXu` help). The break hides leakage that appears on break-free spectra.
%[text] - The "verified for alpha in \[-2, 6\]" statement that rested on them is withdrawn. \
%[text] The table is kept in the pre-reg project (`docs/README_IRASAThreeWayAgreement_v001.md`) as a record of the legacy code, and can be reproduced with `Resampling="legacy"`.
%[text] ## What has been checked for v2.0.0
%[text] All on `noiseXu` output, N = 3000, fs = 240 Hz, band 1-20 Hz, in the toolbox test suite (`tests/testEstimateIRASA.m`):
checked = table( ...
    ["10 Hz sine (0.5 SD) on alpha 1.5 noise"; "Peak-free noise, alpha 2"; "Break-free noise (Phi = 1), alpha 4.5, Prewhiten = 2"; "Band reaching above fs/(2*max(Hset))"], ...
    ["Peak left in pFractal: about 1.5 dB (legacy: 14.7 dB; raw spectrum: 14.9 dB)"; "Median within 0.2 of 2, and within 0.15 of legacy"; "Median within 0.3 of 4.5, pmtm and Welch"; "Error, not a silent answer"], ...
    VariableNames=["Signal", "Result"]);
disp(checked)
%%
%[text] **Not yet re-measured under the new default:** bias across the full alpha range, and the alpha 7 behaviour. Until the validation grid is rerun (`validation/`), treat the list above as the whole of the evidence. At steep alpha use `estimateWhittleAlpha`, or pass its `kUsed` as `Prewhiten`.
%[text] ## See it remove a peak
%[text] The same signal through the default and through the legacy path. Peak excess is the largest value of `pFractal` in 9.5-10.5 Hz over its median in the neighbouring bands, in dB.
fs = 240;  N = 3000;  t = (0:N-1)' / fs;
rng(20261109);
x = noiseXu(N, 1.5, 1, fs) + 0.5 * sin(2*pi*10*t + 2*pi*rand);
peakDb = @(f, p) 10 * log10(max(p(f >= 9.5 & f <= 10.5)) / median(p((f >= 7 & f <= 8.5) | (f >= 12 & f <= 14)), 'omitnan'));

[aNew, ~, pNew, f] = estimateIRASA(x, fs, 1, 20);
warnState = warning('off', 'estimateIRASA:legacyResampling');   % demonstration only: the warning is expected here
[aOld, ~, pOld] = estimateIRASA(x, fs, 1, 20, Resampling="legacy");
warning(warnState);
fprintf("default: alpha %.2f, 10 Hz peak %.1f dB\n", aNew, peakDb(f, pNew));
fprintf("legacy : alpha %.2f, 10 Hz peak %.1f dB\n", aOld, peakDb(f, pOld));

figure;  loglog(f, pOld, f, pNew, 'LineWidth', 1.2);  xlim([1 40]);  grid on;
xlabel('Frequency (Hz)');  ylabel('pFractal');  legend('legacy (v1.x)', 'default (v2.0.0)');
title('A 10 Hz rhythm on 1/f^{1.5} noise');
%%
%[text] ## A second opinion from FieldTrip
%[text] FieldTrip's IRASA uses a Hanning-windowed FFT in place of multitaper but the same resampling-median step, so it is now a like-for-like comparison. A missing FieldTrip gives a flagged result and a warning, never a substitute.
[aHB, ~, ~, ~, ftInfo] = estimateIRASA(x, fs, 1, 20, CompareFieldTrip=true);
if ftInfo.available
    fprintf("estimateIRASA %.3f, FieldTrip %.3f, gap %.3f\n", aHB, ftInfo.alphaFT, ftInfo.gapHBFT);
else
    fprintf("estimateIRASA %.3f (FieldTrip not available: %s)\n", aHB, ftInfo.message);
end
%%
%[text] ## Using it
%[text] - The band must satisfy `fHigh*max(Hset) < fs/2` (with the default `Hset`, `fHigh < fs/3.8`). `pFractal` is NaN above `fs/(2*max(Hset))`.
%[text] - `ComparePMTM=true` returns the direct slope with no resampling. A gap between it and `alphaEst` now means peaks were present and removed.
%[text] - `Resampling="legacy"` exists only to reproduce results from v1.2.0 and earlier. It warns on every call. \
%[text] ## References
%[text] Wen, H., and Liu, Z. (2016). Separating fractal and oscillatory components in the power spectrum of neurophysiological signal. *Brain Topography*, 29(1), 13-26.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
