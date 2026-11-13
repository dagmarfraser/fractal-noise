function n = shapeNoise(resid, fs, alpha, options)
%SHAPENOISE Impose an empirical amplitude spectrum onto a chosen phase carrier.
%
% Toolbox generalisation of shapeXu_local, a private local function inside
% PowerLawSimulationPreReg's generateLoopClosureNoise_v003.m (~line 129).
% The magnitude-imposition / sigma-restoration / conjugate-symmetry logic
% is unchanged and source-agnostic. What changes is the phase carrier: the
% original hardcoded generateCustomNoise_v004 (Xu-only); this version
% exposes PhaseSource so the carrier composes with dsp.colorednoise,
% MATLAB's built-in pinknoise, plain white noise, or any custom generator,
% not just Xu.
%
% CONSTRUCTION with Prewhiten = 0 (unchanged from source):
%   F_shaped = |FFT(resid)| .* exp(1i*angle(FFT(carrier)))
%   n        = real(ifft(F_shaped)), rescaled so std(n) == std(resid)
% CONSTRUCTION with Prewhiten = k > 0 (new in v2.0.0):
%   the same swap applied to diff(resid, k) and diff(carrier, k), then
%   integrated k times (cumsum), centred on mean(resid) and rescaled so
%   std(n) == std(resid).
% Conjugate symmetry is automatic: |FFT(resid)| is symmetric and
% angle(FFT(carrier)) is antisymmetric for a real carrier, so F_shaped is
% a valid real-signal spectrum. real() guards residual numerical asymmetry.
%
% SYNTAX:
%   n = shapeNoise(resid, fs, alpha)
%   n = shapeNoise(resid, fs, alpha, PhaseSource="xu")
%   n = shapeNoise(resid, fs, alpha, PhaseSource="white")
%   n = shapeNoise(resid, fs, alpha, PhaseSource="pinknoise")
%   n = shapeNoise(resid, fs, alpha, PhaseSource="dsp")
%   n = shapeNoise(resid, fs, alpha, PhaseSource=@myCarrierFcn)
%   n = shapeNoise(___, Prewhiten=k)             % k = "auto" or 0..3
%
% INPUTS:
%   resid - Empirical residual time series, real column vector. Its
%           amplitude spectrum (magnitude only) is what gets imposed.
%   fs    - Sampling rate in Hz, passed through to the phase carrier.
%   alpha - Target spectral exponent for the phase carrier. Required by
%           PhaseSource "xu" and "dsp"; ignored (with a warning if it
%           looks like the caller expected it to matter) by "white" and
%           "pinknoise", both of which are fixed-colour carriers.
%   PhaseSource - Name-value. One of:
%       "xu"        (default) noiseXu(M, alpha, sig, fs) -- genuine
%                   fractional-phase carrier, matches source behaviour
%                   exactly when left at default.
%       "white"     Plain randn(M,1)*sig. Random phase; equivalent
%                   carrier to the project's separate 'fftnoise' model,
%                   offered here for composability/comparison.
%       "pinknoise" MATLAB's built-in pinknoise function (Signal
%                   Processing Toolbox). Fixed alpha~1 carrier; alpha
%                   input is not used to shape it.
%       "dsp"       dsp.ColoredNoise (DSP System Toolbox), using its
%                   InverseFrequencyPower property (continuous, range
%                   [-2, 2]). Requested alpha outside that range is
%                   clamped to the nearest bound, with a warning.
%       function handle @(M, fs, alpha) -> real Mx1 vector. Most general
%                   option; must return a real column vector of length M.
%   Prewhiten - Name-value, "auto" or an integer 0..3. Default "auto".
%       Number of times resid and the carrier are differenced before the
%       magnitude/phase swap. "auto" uses k = min(3, floor((alpha + 1)/2))
%       (0 below alpha = 1), so alpha must then be a rough exponent of
%       resid whatever the carrier.
%       WHY: abs(fft(resid)) of a steep series whose ends do not meet is,
%       above the lowest bins, mostly end-mismatch leakage (f^-2). With
%       Prewhiten = 0 the surrogate inherits that magnitude and stays
%       steep only when the carrier's own end jump supplies matching
%       phases. Measured (fs 133, N 3000, 1-20 Hz, 100 seeds, "xu"
%       carrier): 21, 23, 13 and 10 percent of surrogates at alpha 3.5,
%       4, 4.5 and 5 came out more than 0.5 shallower than resid; with
%       "auto", none did. With Prewhiten = k > 0 the surrogate's
%       k-th difference has the magnitude spectrum of resid's k-th
%       difference; abs(fft(n)) itself no longer equals abs(fft(resid)).
%       Prewhiten = 0 is the construction of v1.2.0 and earlier, and of
%       the source shapeXu_local, bit for bit.
%
% OUTPUT:
%   n - Real surrogate column vector, length numel(resid), with
%       std(n) numerically equal to std(resid).
%
% NOTES:
%   - "dsp" and "pinknoise" are optional-toolbox paths. noiseXu ("xu",
%     the default) and "white" require no toolbox beyond base MATLAB and
%     are the guaranteed-available fallback.
%   - Whether a given carrier's phase actually survives the magnitude
%     swap to still deliver the intended alpha is an empirical question
%     (see PowerLawSimulationPreReg Finding #67 and the
%     alpha-preservation checks referenced there). This function performs
%     the construction; it does not itself validate alpha preservation.
%     For the "xu" carrier see Prewhiten above, and
%     validation/checkShapeNoiseMechanism_v001.m.
%
% See also: noiseXu, estimateIRASA
%
% Toolbox: fractalnoise
% Extracted from: src/functions/generateLoopClosureNoise_v003.m (shapeXu_local)
% (PowerLawSimulationPreReg, session 52)

arguments
    resid (:, 1) double
    fs (1, 1) double {mustBePositive}
    alpha (1, 1) double
    options.PhaseSource = "xu"
    options.Prewhiten = "auto"
end

resid = resid(:);
M = numel(resid);
sig = std(resid, 0, 1);

if sig == 0
    n = zeros(M, 1);
    return
end

k = prewhitenOrder(options.Prewhiten, alpha, M);   % validated before any random numbers are drawn
carrier = shapeNoiseCarrier(options.PhaseSource, M, fs, alpha, sig);

if ~isreal(carrier) || ~iscolumn(carrier) || numel(carrier) ~= M
    error("shapeNoise:InvalidCarrier", ...
        "Phase carrier must return a real column vector of length %d.", M);
end

if k == 0
    empiricalMagnitude = abs(fft(resid));
    carrierSpectrum = fft(carrier);
    shapedSpectrum = empiricalMagnitude .* exp(1i * angle(carrierSpectrum));
    n = real(ifft(shapedSpectrum));

    carrierStd = std(n, 0, 1);
    if carrierStd > 0
        n = n * (sig / carrierStd);
    end
else
    n = real(ifft(abs(fft(diff(resid, k))) .* exp(1i * angle(fft(diff(carrier, k))))));
    for q = 1:k
        n = cumsum([0; n]);
    end
    shapedStd = std(n, 0, 1);
    if shapedStd == 0 || ~isfinite(shapedStd)
        error("shapeNoise:PrewhitenDegenerate", ...
            "Prewhiten=%d gave a surrogate with std %g. Use a smaller Prewhiten.", k, shapedStd);
    end
    n = (n - mean(n)) * (sig / shapedStd) + mean(resid);
end

end

function k = prewhitenOrder(prewhiten, alpha, M)
%PREWHITENORDER Resolve the Prewhiten option to an integer differencing order.
if (isstring(prewhiten) || ischar(prewhiten)) && isscalar(string(prewhiten)) && lower(string(prewhiten)) == "auto"
    if ~isfinite(alpha)
        error("shapeNoise:PrewhitenAutoNeedsAlpha", ...
            "Prewhiten=""auto"" needs a finite alpha; got %g. Pass Prewhiten=k explicitly.", alpha);
    end
    k = min(3, max(0, floor((alpha + 1) / 2)));
elseif isnumeric(prewhiten) && isscalar(prewhiten) && ismember(prewhiten, 0:3)
    k = double(prewhiten);
else
    error("shapeNoise:InvalidPrewhiten", "Prewhiten must be ""auto"" or an integer 0..3.");
end
if k > 0 && M - k < 8
    error("shapeNoise:TooShortForPrewhiten", ...
        "resid has %d samples; Prewhiten=%d needs at least %d.", M, k, k + 8);
end
end

function carrier = shapeNoiseCarrier(phaseSource, M, fs, alpha, sig)
%SHAPENOISECARRIER Dispatch to the requested phase-carrier generator.
%   sig is the residual's std, passed through as the carrier's target
%   amplitude for carriers where that choice matters for exact
%   reproducibility (xu); it has no effect on the final output for any
%   carrier, since only phase survives into the shaped spectrum.

if isa(phaseSource, "function_handle")
    carrier = phaseSource(M, fs, alpha);
    return
end

phaseSource = lower(string(phaseSource));

switch phaseSource
    case "xu"
        carrier = noiseXu(M, alpha, sig, fs);

    case "white"
        carrier = randn(M, 1);

    case "pinknoise"
        if exist("pinknoise", "file") ~= 2 && exist("pinknoise", "builtin") == 0
            error("shapeNoise:PinknoiseUnavailable", ...
                "PhaseSource=""pinknoise"" requires MATLAB's built-in pinknoise " + ...
                "function (Signal Processing Toolbox), which is not on the path.");
        end
        if isfinite(alpha) && abs(alpha - 1) > 0.3
            warning("shapeNoise:PinknoiseAlphaMismatch", ...
                "PhaseSource=""pinknoise"" is a fixed alpha~1 carrier; requested " + ...
                "alpha=%.2f is not used to shape it.", alpha);
        end
        carrier = pinknoise(M);

    case "dsp"
        if exist("dsp.ColoredNoise", "class") ~= 8
            error("shapeNoise:DspUnavailable", ...
                "PhaseSource=""dsp"" requires dsp.ColoredNoise (DSP System " + ...
                "Toolbox), which is not installed.");
        end
        clampedAlpha = min(max(alpha, -2), 2);
        if clampedAlpha ~= alpha
            warning("shapeNoise:DspAlphaClamped", ...
                "PhaseSource=""dsp"" (dsp.ColoredNoise) supports " + ...
                "InverseFrequencyPower in [-2, 2]; requested alpha=%.2f " + ...
                "clamped to %.2f.", alpha, clampedAlpha);
        end
        generator = dsp.ColoredNoise(InverseFrequencyPower=clampedAlpha, ...
            SamplesPerFrame=M, NumChannels=1);
        carrier = generator();

    otherwise
        error("shapeNoise:UnknownPhaseSource", ...
            "Unrecognised PhaseSource ""%s"". Use ""xu"", ""white"", " + ...
            """pinknoise"", ""dsp"", or a function handle.", phaseSource);
end

end
