function [alphaEst, sigmaEst, pFractal, fVec, ftInfo, pmInfo] = estimateIRASA(x, fs, fLow, fHigh, options)
%ESTIMATEIRASA IRASA spectral-exponent and noise-magnitude estimator.
%
% Toolbox wrap of PowerLawSimulationPreReg's iraAlphaSigma_v003.m (not
% v001 -- v003 adds the pFractal/fVec outputs needed for Lorentzian-knee
% work at no downside; same core estimator otherwise). Core algorithm
% (resampling-median IRASA on pmtm PSDs, Wen & Liu 2016 mechanism) is
% unchanged from source.
%
% SYNTAX:
%   [alphaEst, sigmaEst] = estimateIRASA(x, fs, fLow, fHigh)
%   [alphaEst, sigmaEst, pFractal, fVec] = estimateIRASA(x, fs, fLow, fHigh)
%   [___] = estimateIRASA(___, Hset=customHset)
%   [___, ftInfo] = estimateIRASA(___, CompareFieldTrip=true)
%   [___, ftInfo, pmInfo] = estimateIRASA(___, ComparePMTM=true)
%   [___] = estimateIRASA(___, Prewhiten=k)          % recommended for steep spectra
%   [___] = estimateIRASA(___, PSD="welch")
%
% INPUTS:
%   x     - Time series, real column vector
%   fs    - Sampling rate in Hz
%   fLow  - Lower bound of the log-log fitting band, Hz
%   fHigh - Upper bound of the log-log fitting band, Hz
%   Hset  - Name-value. Resampling ratio set for the IRASA median.
%           Default 1.1:0.05:1.9 (matches source default).
%   CompareFieldTrip - Name-value logical, default false. If true, also
%           runs FieldTrip's ft_freqanalysis(method='irasa') on the same
%           signal and band for comparison, returned in ftInfo.
%   ComparePMTM - Name-value logical, default false. If true, also fits
%           a direct log-log slope to the un-resampled PSD (pmtm, or
%           Welch with PSD="welch"; recoloured if Prewhiten > 0) (no
%           IRASA resampling-median step at all) over the same [fLow,
%           fHigh] band, returned in pmInfo. This is the comparison that
%           actually tests the resampling-median step in isolation:
%           CompareFieldTrip varies the PSD estimator (pmtm vs Hanning
%           FFT) while keeping resampling-median fixed, so HB-vs-FT
%           agreement cannot rule out both sharing a resampling-median
%           blind spot. ComparePMTM holds the PSD estimator fixed
%           (pmtm, same as HB) and removes resampling-median entirely,
%           isolating that step as the variable under test -- see
%           README_IRASAThreeWayAgreement_v001.md Section 2.
%   PSD   - Name-value, "pmtm" (default, unchanged from source) or
%           "welch" (Hann, 50% overlap, segment min(1024, floor(N/2))
%           samples at fs; resampled copies keep the same segment
%           duration). Applies to every PSD inside IRASA and to pRaw.
%   Prewhiten - Name-value, integer k in 0..3, default 0. Differences x
%           k times before IRASA, then recolours pFractal and pRaw by
%           1/(2 - 2cos(2*pi*f/fs))^k, so alphaEst is still the exponent
%           of x. Removes the leakage that biases IRASA on steep,
%           break-free spectra (see RANGE). Pick k = floor((alpha + 1)/2)
%           from a rough alpha, or take the kUsed output of
%           estimateWhittleAlpha. Errors if combined with CompareFieldTrip.
%
% OUTPUTS:
%   alphaEst - Spectral exponent (homebrew IRASA estimate)
%   sigmaEst - Noise magnitude, sqrt(integrated fractal PSD * 2), where
%              the integral is taken over the full available frequency
%              range of the fractal PSD estimate (all valid bins), NOT
%              restricted to [fLow, fHigh] -- that band is used only for
%              the alpha log-log fit. Do not expect sigmaEst to equal
%              std(x); it is a spectral-magnitude estimate, not a
%              reconstruction of the time-domain standard deviation.
%   pFractal - Aperiodic PSD, linear units, geometric median across
%              h-pairs (same length as fVec)
%   fVec     - Frequency axis, Hz
%   ftInfo   - Struct, populated only when CompareFieldTrip=true:
%       .available - logical, whether FieldTrip was found on the path
%       .alphaFT   - FieldTrip's recovered alpha (NaN if unavailable or
%                    the fit failed)
%       .gapHBFT   - abs(alphaEst - alphaFT) (NaN if either is NaN)
%       .message   - string, empty on success; explains why
%                    available=false or alphaFT=NaN otherwise
%       When CompareFieldTrip=false (default), ftInfo is returned as a
%       struct with available=false and message="CompareFieldTrip not requested".
%   pmInfo   - Struct, populated only when ComparePMTM=true:
%       .available - logical, always true when requested (no external
%                    dependency; pmtm ships with Signal Processing
%                    Toolbox, already required by the core estimator)
%       .alphaPM   - direct pmtm log-log slope estimate, no resampling
%                    (NaN if fewer than 3 valid bins in [fLow, fHigh])
%       .gapHBPM   - abs(alphaEst - alphaPM) (NaN if either is NaN)
%       .message   - string, empty on success; explains why alphaPM=NaN
%                    otherwise
%       When ComparePMTM=false (default), pmInfo is returned as a
%       struct with available=false and message="ComparePMTM not requested".
%
% RANGE (validation/irasaOptions/checkIRASAOptions_v001, N 3000, band
% 1-20 Hz, fs 100 and 240, alpha 1.5-5.5, 10 series per cell):
%   The default (pmtm, no prewhitening) depends on the noise, not just
%   on alpha. On noise with a low-frequency break (noiseXu default Phi)
%   it is within +/-0.1 to alpha 5.5. On break-free noise (exact
%   integrated fGn, or noiseXu Phi = 1) it overshoots from alpha 3.5 at
%   fs 240 (+0.5 at 3.5, +1.7 to +2.1 at 4.5) and collapses at alpha 5.5,
%   fs 100 (-2.1). PSD="welch" alone halves but does not remove this
%   (+0.4 to +0.6 at alpha >= 4.5). Prewhiten=k, with either PSD, keeps
%   the bias within +/-0.23 on all three generators, the same as
%   estimateWhittleAlpha. For real data, where you do not know whether
%   there is a break, use Prewhiten for alpha above about 3, or use
%   estimateWhittleAlpha.
%   ComparePMTM is NOT an independent check of the default: it fits the
%   same un-resampled pmtm PSD, so it shares the leakage, and agreement
%   with alphaEst does not validate it on steep spectra. CompareFieldTrip
%   changes only the taper (one whole-record Hann) and keeps IRASA.

% FIELDTRIP DEPENDENCY POLICY:
%   CompareFieldTrip defaults to false, so estimateIRASA has no hard
%   FieldTrip dependency for its core (homebrew) estimate. Requesting
%   CompareFieldTrip=true without FieldTrip on the path does NOT error --
%   it returns ftInfo.available=false with an explanatory message and a
%   visible warning banner (Fail Loud, Never Fake: the degraded mode is
%   disclosed, not silently skipped). The `which('ft_defaults')` gate
%   (not `which('ft_freqanalysis')`) is used to detect FieldTrip
%   availability, following the corrected pattern in
%   benchmarkIRASA_v004.m: ft_freqanalysis.m sits in the FieldTrip root
%   and `which` can find it before ft_defaults has run, but FieldTrip's
%   internal private/ utility calls will still fail until ft_defaults has
%   added them to the path. Checking for ft_defaults itself is the
%   correct availability test.
%
% See also: noiseXu, checkSpectralKnee, estimateWhittleAlpha
%
% Toolbox: fractalnoise
% Extracted from: src/functions/iraAlphaSigma_v003.m (algorithm inlined
% below as a local function -- the toolbox does not depend on the
% original file existing outside its own tree)
% (PowerLawSimulationPreReg, session 52)
% Citation: Wen, H., & Liu, Z. (2016). Separating Fractal and
%   Oscillatory Components in the Power Spectrum of Neurophysiological
%   Signal. Brain Topography, 29(1), 13-26.
%   https://doi.org/10.1007/s10548-015-0448-0

arguments
    x (:, 1) double
    fs (1, 1) double {mustBePositive}
    fLow (1, 1) double {mustBePositive}
    fHigh (1, 1) double {mustBePositive}
    options.Hset (1, :) double = 1.1:0.05:1.9
    options.CompareFieldTrip (1, 1) logical = false
    options.ComparePMTM (1, 1) logical = false
    options.PSD (1, 1) string {mustBeMember(options.PSD, ["pmtm", "welch"])} = "pmtm"
    options.Prewhiten (1, 1) double {mustBeInteger, mustBeNonnegative, mustBeLessThanOrEqual(options.Prewhiten, 3)} = 0
end

if options.CompareFieldTrip && options.Prewhiten > 0
    error("estimateIRASA:prewhitenFieldTrip", "%s", ...
        "CompareFieldTrip runs FieldTrip on the raw series, so it cannot be " + ...
        "compared with a prewhitened estimate. Use Prewhiten=0 with CompareFieldTrip.");
end
if numel(x) - options.Prewhiten < 8
    error("estimateIRASA:short", "Series too short (N = %d) for Prewhiten = %d.", numel(x), options.Prewhiten);
end

[alphaEst, sigmaEst, pFractal, fVec, pRaw] = estimateIRASACore(x, fs, fLow, fHigh, options.Hset, options.PSD, options.Prewhiten);

if options.CompareFieldTrip
    ftInfo = estimateIRASAFieldTripCompare(x, fs, fLow, fHigh, options.Hset, alphaEst);
else
    ftInfo = struct("available", false, "alphaFT", NaN, "gapHBFT", NaN, ...
        "message", "CompareFieldTrip not requested");
end

if options.ComparePMTM
    pmInfo = estimateIRASAPmtmCompare(fVec, pRaw, fLow, fHigh, alphaEst);
else
    pmInfo = struct("available", false, "alphaPM", NaN, "gapHBPM", NaN, ...
        "message", "ComparePMTM not requested");
end

end

function pmInfo = estimateIRASAPmtmCompare(fVec, pRaw, fLow, fHigh, alphaEst)
%ESTIMATEIRASAPMTMCOMPARE Direct pmtm log-log slope, no resampling-median.
%   Uses the same un-resampled pmtm PSD (pRaw, fVec) that
%   estimateIRASACore already computes once at h=1 before entering the
%   resampling-median loop -- no extra pmtm call. Isolates the
%   resampling-median step as the variable under test against alphaEst
%   (see README_IRASAThreeWayAgreement_v001.md Section 2, "PM").

fitMask = fVec >= fLow & fVec <= fHigh & ~isnan(pRaw) & pRaw > 0;
if sum(fitMask) < 3
    pmInfo = struct("available", true, "alphaPM", NaN, "gapHBPM", NaN, ...
        "message", "PMTM fit failed: fewer than 3 valid bins in fitting band");
    return
end

pFit = polyfit(log10(fVec(fitMask)), log10(pRaw(fitMask)), 1);
alphaPM = -pFit(1);

if isfinite(alphaEst) && isfinite(alphaPM)
    gapHBPM = abs(alphaEst - alphaPM);
else
    gapHBPM = NaN;
end

pmInfo = struct("available", true, "alphaPM", alphaPM, "gapHBPM", gapHBPM, ...
    "message", "");

end

function ftInfo = estimateIRASAFieldTripCompare(x, fs, fLow, fHigh, hset, alphaEst)
%ESTIMATEIRASAFIELDTRIPCOMPARE Optional FieldTrip IRASA comparison panel.
%   Availability gate follows benchmarkIRASA_v004.m: check ft_defaults,
%   not ft_freqanalysis (see FIELDTRIP DEPENDENCY POLICY above).

if isempty(which("ft_defaults"))
    warning("estimateIRASA:NoFieldTrip", "%s", ...
        "CompareFieldTrip=true requested but FieldTrip not found on path " + ...
        "(checked for ft_defaults). Returning ftInfo.available=false. " + ...
        "Install from https://www.fieldtriptoolbox.org and addpath the " + ...
        "FieldTrip root folder to enable this comparison.");
    ftInfo = struct("available", false, "alphaFT", NaN, "gapHBFT", NaN, ...
        "message", "FieldTrip not found on path (ft_defaults absent)");
    return
end

ft_defaults;
if isempty(which("ft_freqanalysis"))
    warning("estimateIRASA:FieldTripIncomplete", "%s", ...
        "ft_defaults ran but ft_freqanalysis still not found -- check your " + ...
        "FieldTrip installation. Returning ftInfo.available=false.");
    ftInfo = struct("available", false, "alphaFT", NaN, "gapHBFT", NaN, ...
        "message", "ft_defaults ran but ft_freqanalysis not found");
    return
end

N = numel(x);
ftData.trial = {x(:)'};
ftData.time = {(0:N - 1) / fs};
ftData.fsample = fs;
ftData.label = {'ch1'};
ftData.sampleinfo = [1, N]; %#ok<STRNU> ftData used inside evalc() string below; analyzer cannot trace it. Also avoids FieldTrip's fixsampleinfo reconstruction warning.
% FieldTrip's internal validation requires char vectors here, not string
% objects -- {"ch1"} throws "Cell array input must be a cell array of
% character vectors" inside ft_freqanalysis. Single-quoted char literal
% is required, not a style choice.

ftCfg = [];
ftCfg.method = 'irasa';
ftCfg.output = 'fractal';
ftCfg.foilim = [fLow fHigh];
ftCfg.taper = 'hanning';
ftCfg.pad = 'nextpow2';
ftCfg.hset = hset;
ftCfg.keeptrials = 'no';

alphaFT = NaN;
message = "";
try
    [~, ftFreq] = evalc("ft_freqanalysis(ftCfg, ftData)");
    fFT = ftFreq.freq;
    pFT = squeeze(ftFreq.powspctrm);
    fitMask = fFT >= fLow & fFT <= fHigh & ~isnan(pFT) & pFT > 0;
    if sum(fitMask) >= 3
        pFit = polyfit(log10(fFT(fitMask)), log10(pFT(fitMask)), 1);
        alphaFT = -pFit(1);
    else
        message = "FieldTrip fit failed: fewer than 3 valid bins in fitting band";
    end
catch ME
    message = "FieldTrip call failed: " + string(ME.message);
end

if isfinite(alphaEst) && isfinite(alphaFT)
    gapHBFT = abs(alphaEst - alphaFT);
else
    gapHBFT = NaN;
end

ftInfo = struct("available", true, "alphaFT", alphaFT, "gapHBFT", gapHBFT, ...
    "message", message);

end

function [alphaEst, sigmaEst, pFractal, fVec, pRaw] = estimateIRASACore(x, fs, fLow, fHigh, hset, psdMethod, k)
%ESTIMATEIRASACORE Inlined copy of iraAlphaSigma_v003's IRASA algorithm.
%   Bundled here, not called externally, so the toolbox has no dependency
%   on PowerLawSimulationPreReg's src/functions/ directory existing on
%   the path. With psdMethod="pmtm" and k=0 (the defaults) the algorithm
%   is unchanged from source: resampling-median IRASA on pmtm PSDs (Wen &
%   Liu 2016 mechanism), and the output is bit-identical to it.
%
%   psdMethod="welch": every PSD (original, up- and down-sampled) is a
%   Hann-window Welch estimate with 50% overlap and the same segment
%   DURATION, so all three share one frequency resolution. Segment length
%   at fs is min(1024, floor(N/2)) samples (the stage-1 validation grid).
%
%   k > 0: x is differenced k times before IRASA, and pFractal and pRaw
%   are recoloured by 1/(2 - 2cos(2*pi*f/fs))^k before the fit, so
%   alphaEst is the exponent of x itself. The DC bin is set to NaN
%   (division by zero there), so sigmaEst excludes it.
%
%   pRaw (fifth output) is the un-resampled h=1 PSD, computed here
%   regardless (needed for fVec). Kept and returned so ComparePMTM can
%   reuse it without a second PSD call.

for j = 1:k
    x = diff(x);
end

[pRaw, fVec] = psdLocal(x, fs, psdMethod, numel(x));
nF = numel(fVec);
geoMeans = nan(nF, numel(hset));

for hi = 1:numel(hset)
    h = hset(hi);

    xUp = resample(x, round(h * 1000), 1000);
    [pU, fU] = psdLocal(xUp(:), fs * h, psdMethod, numel(x));
    pUI = interp1(fU, pU, fVec, 'linear', NaN);

    xDn = resample(x, 1000, round(h * 1000));
    [pD, fD] = psdLocal(xDn(:), fs / h, psdMethod, numel(x));
    pDI = interp1(fD, pD, fVec, 'linear', NaN);

    geoMeans(:, hi) = sqrt(pUI .* pDI);
end

pFractal = median(geoMeans, 2, 'omitnan');

if k > 0
    gain = (2 - 2 * cos(2 * pi * fVec / fs)) .^ k;
    gain(gain == 0) = NaN;
    pFractal = pFractal ./ gain;
    pRaw = pRaw ./ gain;
end

fitMask = fVec >= fLow & fVec <= fHigh & ~isnan(pFractal) & pFractal > 0;
if sum(fitMask) < 3
    alphaEst = NaN;
    sigmaEst = NaN;
    return
end
pFit = polyfit(log10(fVec(fitMask)), log10(pFractal(fitMask)), 1);
alphaEst = -pFit(1);

validMask = ~isnan(pFractal) & pFractal > 0;
noiseVar = trapz(fVec(validMask), pFractal(validMask)) * 2;
sigmaEst = sqrt(max(noiseVar, 0));

end

function [p, f] = psdLocal(x, fsLocal, psdMethod, nOrig)
%PSDLOCAL PSD of one (possibly resampled) series. pmtm path is the source call.
if psdMethod == "pmtm"
    [p, f] = pmtm(x, 4, [], fsLocal);
else
    % Segment length at the original fs is min(1024, floor(nOrig/2)) samples;
    % resampled series scale it by their length ratio, keeping the same
    % segment DURATION and hence the same frequency resolution.
    L = min(numel(x), round(min(1024, floor(nOrig / 2)) * numel(x) / nOrig));
    [p, f] = pwelch(x, hann(L, 'periodic'), floor(L / 2), max(256, 2 ^ nextpow2(L)), fsLocal);
end
end
