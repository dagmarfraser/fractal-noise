function [n, info] = shapeNoise(resid, fs, alpha, options)
%SHAPENOISE Surrogate noise with the amplitude spectrum of a measured residual.
%
%   [n, info] = shapeNoise(resid, fs, alpha) returns an end-matched
%   random-phase surrogate of resid (the default since v2.0.0):
%     1. Segment: the longest length L at which every channel of resid has
%        a sub-segment whose ends meet: its wrap-around step, and each of its
%        two wrapped second differences, is no larger than the largest step
%        expected among the segment's own L-1 steps (L-2 second differences).
%        That level is q(n) = chi2inv(1 - 1/n, 1), about 2 ln n, times the
%        mean squared step: the wrap must not stand out from the segment's
%        own steps, and no tolerance is chosen (SegmentRule "detect").
%        Segment measures from Schreiber & Schmitz 2000, Physica D 142,
%        Eqs. 16-17; longest-first search after the stepwise stopping of
%        Nakamura, Small & Hirata 2006, Phys Rev E 74, 026205. Each channel
%        takes its own best segment of that common length.
%     2. Surrogate: |FFT(segment)| on uniformly random phases, mean removed,
%        scaled to std(segment).
%     3. Alpha selection: a draw is accepted when the estimateWhittleAlpha
%        exponent of the surrogate is within AlphaTol of the target; up to
%        MaxDraw draws per channel.
%   The surrogate has length L, not numel(resid): a Fourier surrogate of a
%   steep series whose ends do not meet carries the end jump as broadband
%   leakage (Schreiber & Schmitz 2000), so the whole record cannot be
%   matched. Whatever you add the noise to must have length L too.
%   Under "detect" the wrap adds at most q(L-1)/(L-1) to the first-
%   difference energy (about 4% at L 280, 1.4% at 1000, 0.5% at 3000).
%
%   Validated on real template-subtracted movement residuals (Fraser,
%   Hickman, Cook; PowerLawSimulationPreReg runner v024): filtered velocity
%   and acceleration noise 1.00-1.02 of the residual's (median), with alpha
%   held by the selection step. The 1.x construction failed that check
%   (1-27% of surrogates more than 0.5 too shallow at alpha 3-6; see
%   validation/checkShapeNoiseEndMatched_v001.m and
%   validation/checkShapeNoiseMechanism_v001.m).
%
%   RANGE: steep, short records may have no end-matched segment. In
%   noiseXu references (validation/checkShapeNoiseEndMatched_v001.m) the
%   call errored ("shapeNoise:noSegment", or "shapeNoise:segmentAlpha"
%   when the segment is too short to estimate) for N 3000 at fs 133 in 6%
%   of references at alpha 5.5 and 30% at 6, none below; for N 1000 at
%   fs 240 in 1% at alpha 4, 16% at 4.5, 58% at 5, 79% at 5.5, 82% at 6.
%   Median L/N fell from 0.97 (alpha 3) to 0.53 (alpha 6) at N 3000.
%
%   TEMPLATE-SUBTRACTED DATA: pass only the part of the residual the
%   template actually covered. A stretch left at exactly zero (an uncovered
%   template tail) is an error here, because the segment search would
%   otherwise take it and feed a step into the surrogate.
%
% SYNTAX:
%   [n, info] = shapeNoise(resid, fs, [])          % target alpha = segment's own (recommended)
%   [n, info] = shapeNoise(resid, fs, alpha)       % explicit target, scalar or 1-by-C
%   [n, info] = shapeNoise(___, AlphaTol=0.3, MaxDraw=20, MinLength=32, AlphaBand=[1 20])
%   [n, info] = shapeNoise(___, SegmentRule="budget", SegmentBudget=0.01)   % wrap adds <= 1% derivative energy
%   [n, info] = shapeNoise(___, SegmentRule="tol", FairShareTol=1)          % wrap <= Tol x mean terms
%   n = shapeNoise(resid, fs, alpha, Method="legacy")                      % 1.x construction; warns
%   n = shapeNoise(resid, fs, alpha, Method="legacy", PhaseSource="white") % legacy carriers
%
% INPUTS:
%   resid - Residual, real N-by-C matrix (one column per channel, e.g. the
%           two axes of a movement) or a vector. "legacy" takes a vector.
%   fs    - Sampling rate in Hz.
%   alpha - Target exponent for alpha selection: [] (default use: each
%           channel's segment alpha by estimateWhittleAlpha over AlphaBand),
%           a scalar, or 1-by-C. Under "legacy" the carrier's alpha.
%   Method       - "endMatched" (default) or "legacy".
%   AlphaTol     - Accept a draw when |alpha(draw) - target| <= AlphaTol (0.3).
%   MaxDraw      - Draws per channel before giving up (20). A channel with
%                  no accepted draw is returned as NaN, with a warning
%                  ("shapeNoise:noAcceptedDraw") and info.accepted false.
%   SegmentRule  - How a segment's ends must meet:
%                  "detect" (default) wrap step <= q(L-1) and each wrapped
%                          second difference <= q(L-2) times the segment's
%                          mean squared first / second difference;
%                  "budget" the wrap adds at most SegmentBudget of the
%                          segment's first- and second-difference energy;
%                  "tol"    wrap step and mean wrapped second difference
%                          <= FairShareTol times the mean terms (the rule
%                          of the v2.0.0 drafts; stricter at long L).
%   SegmentBudget - Energy share for "budget" (no default; required there).
%   FairShareTol - Tolerance for "tol" (1).
%   MinLength    - Shortest segment searched (32 samples). No segment at
%                  or above it is an error ("shapeNoise:noSegment").
%   AlphaBand    - Band in Hz for estimateWhittleAlpha
%                  ([1, min(20, fs/2 - 1)]).
%   PhaseSource  - "legacy" only: "xu" (default), "white", "pinknoise",
%                  "dsp" or a function handle @(M, fs, alpha).
%
% OUTPUTS:
%   n    - L-by-C surrogate ("endMatched"), or numel(resid)-by-1 ("legacy").
%   info - "endMatched": SegmentRule, i0, i1 (1-by-C segment start and end
%          in resid), L, Lfrac (L/N), cJump, cSlip, cSlipMax (squared wrap
%          step, mean and larger wrapped second difference, over the
%          segment's mean terms; 1-by-C), qJ, qS (the "detect" levels at L),
%          alphaTarget, nDraw, accepted, dAlpha, nEstFail (draws whose alpha
%          estimate errored, counted as rejected and warned about as
%          "shapeNoise:estimateFailed"; 1-by-C). "legacy": Method only.
%
% LEGACY ("legacy", the 1.x behaviour, bit for bit):
%   F = |FFT(resid)| .* exp(1i*angle(FFT(carrier))), n = real(ifft(F))
%   rescaled to std(resid). The carrier's phases decide whether the
%   surrogate stays steep; at alpha 3-6 a tail of surrogates comes out
%   shallow. Kept only to reproduce earlier results; warns on every call
%   ("shapeNoise:legacyConstruction").
%
% See also: estimateWhittleAlpha, noiseXu, estimateIRASA
%
% Toolbox: fractalnoise
% Extracted from: PowerLawSimulationPreReg src/endToEndSegment_v006.m
% (rules "detect", "budget", "tol"; chi2inv(1 - p, 1) = 2*erfcinv(p)^2, so
% no Statistics Toolbox) and src/functions/endMatchedNoise_v001.m, with
% alphaPWLW_v001 replaced by its toolbox copy estimateWhittleAlpha; legacy
% path from src/functions/generateLoopClosureNoise_v003.m (shapeXu_local).

arguments
    resid double {mustBeReal}
    fs (1, 1) double {mustBePositive}
    alpha double = []
    options.Method (1, 1) string {mustBeMember(options.Method, ["endMatched", "legacy"])} = "endMatched"
    options.AlphaTol (1, 1) double {mustBePositive} = 0.3
    options.MaxDraw (1, 1) double {mustBeInteger, mustBePositive} = 20
    options.SegmentRule (1, 1) string {mustBeMember(options.SegmentRule, ["detect", "budget", "tol"])} = "detect"
    options.SegmentBudget double = []
    options.FairShareTol (1, 1) double {mustBePositive} = 1
    options.MinLength (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.MinLength, 5)} = 32
    options.AlphaBand (1, 2) double = [NaN NaN]
    options.PhaseSource = ""
end

if options.Method == "legacy"
    warning("shapeNoise:legacyConstruction", ...
        "Method=""legacy"" is the 1.x construction: at alpha 3-6 a tail of surrogates comes out more than 0.5 " + ...
        "too shallow. It is kept only to reproduce earlier results.");
    if ~isvector(resid) || isempty(resid)
        error("shapeNoise:legacyVector", "Method=""legacy"" takes a non-empty vector residual.");
    end
    if ~isscalar(alpha), error("shapeNoise:legacyAlpha", "Method=""legacy"" needs a scalar carrier alpha."); end
    ps = options.PhaseSource;
    if (isstring(ps) || ischar(ps)) && string(ps) == "", ps = "xu"; end
    n = legacyShape(resid(:), fs, alpha, ps);
    info = struct("Method", "legacy");
    return
end

if ~(isstring(options.PhaseSource) || ischar(options.PhaseSource)) || string(options.PhaseSource) ~= ""
    error("shapeNoise:phaseSourceLegacyOnly", ...
        "PhaseSource applies only to Method=""legacy""; the end-matched surrogate uses uniformly random phases.");
end
if isvector(resid), resid = resid(:); end
[N, C] = size(resid);
if any(~isfinite(resid), "all"), error("shapeNoise:nonFinite", "resid contains NaN or Inf."); end
band = options.AlphaBand;
if any(isnan(band)), band = [1, min(20, fs/2 - 1)]; end
checkZeroTail(resid);
if any(std(resid, 0, 1) == 0), error("shapeNoise:constant", "a channel of resid is constant."); end

if options.SegmentRule == "budget" && ~(isscalar(options.SegmentBudget) && options.SegmentBudget > 0 && options.SegmentBudget < 1)
    error("shapeNoise:segmentBudget", "SegmentRule=""budget"" needs a SegmentBudget in (0, 1).");
end
[i0, i1, seg] = endMatchedSegments(resid, options.SegmentRule, options.FairShareTol, options.SegmentBudget, options.MinLength);
L = seg.L;
if isempty(alpha)
    target = NaN(1, C);
    for ch = 1:C
        try
            target(ch) = estimateWhittleAlpha(resid(i0(ch):i1(ch), ch), fs, band(1), band(2));
        catch cause
            err = MException("shapeNoise:segmentAlpha", ...
                "channel %d: the end-matched segment (L = %d of %d samples) is too short to estimate its alpha: %s", ...
                ch, L, N, cause.message);
            throw(addCause(err, cause));
        end
    end
elseif isscalar(alpha)
    target = repmat(alpha, 1, C);
elseif numel(alpha) == C
    target = alpha(:)';
else
    error("shapeNoise:alphaSize", "alpha must be [], a scalar, or one value per channel (%d).", C);
end
if any(~isfinite(target)), error("shapeNoise:alphaNotFinite", "target alpha is not finite."); end

n = NaN(L, C);  nDraw = zeros(1, C);  accepted = false(1, C);  dAlpha = NaN(1, C);  nEstFail = zeros(1, C);  estMsg = "";
for ch = 1:C
    [n(:, ch), nDraw(ch), accepted(ch), dAlpha(ch), nEstFail(ch), m] = drawEndMatched(resid(i0(ch):i1(ch), ch), ...
        target(ch), fs, band, options.AlphaTol, options.MaxDraw);
    if m ~= "", estMsg = m; end
end
if ~all(accepted)
    warning("shapeNoise:noAcceptedDraw", ...
        "No draw within AlphaTol %.2f of the target in %d tries for channel(s) %s; returned as NaN.", ...
        options.AlphaTol, options.MaxDraw, mat2str(find(~accepted)));
end
info = struct("Method", "endMatched", "SegmentRule", seg.rule, "i0", i0, "i1", i1, "L", L, "Lfrac", L / N, ...
    "cJump", seg.cJump, "cSlip", seg.cSlip, "cSlipMax", seg.cSlipMax, "qJ", seg.qJ, "qS", seg.qS, "alphaTarget", target, "nDraw", nDraw, "accepted", accepted, "dAlpha", dAlpha, "nEstFail", nEstFail);
if any(nEstFail > 0)
    warning("shapeNoise:estimateFailed", ...
        "estimateWhittleAlpha failed on %d draw(s) (counted as rejected; info.nEstFail). Last message: %s", ...
        sum(nEstFail), estMsg);
end
end

function checkZeroTail(S)
% An exact-zero run at either end is an uncovered template stretch, not noise.
for ch = 1:size(S, 2)
    x = S(:, ch);
    nz = find(x ~= 0);
    if isempty(nz), continue, end
    lead = nz(1) - 1;  trail = numel(x) - nz(end);
    if lead >= 3 || trail >= 3
        error("shapeNoise:zeroTail", ...
            "channel %d has %d leading and %d trailing samples of exactly zero (an uncovered template stretch?). " + ...
            "Pass only the covered residual.", ch, lead, trail);
    end
end
end

function [i0, i1, info] = endMatchedSegments(S, rule, tol, budget, minLength)
% endToEndSegment_v006 (rules "detect", "budget", "tol"); chi2inv(1 - p, 1) written as 2*erfcinv(p)^2 (base MATLAB).
[N, C] = size(S);
if N < minLength
    error("shapeNoise:short", "resid has %d samples, fewer than MinLength %d.", N, minLength);
end
e1 = [zeros(1, C); cumsum(diff(S).^2)];  e2 = [zeros(1, C); cumsum(diff(S, 2).^2)];
q = @(n) 2 * erfcinv(1 / n)^2;
i0 = NaN(1, C);  i1 = NaN(1, C);
for L = N:-1:minLength
    a = (1:N - L + 1)';  b = a + L - 1;
    m1 = (e1(b, :) - e1(a, :)) / (L - 1);  m2 = (e2(b - 1, :) - e2(a, :)) / (L - 2);
    cj = (S(a, :) - S(b, :)).^2 ./ m1;
    w1 = (S(a + 1, :) - 2 * S(a, :) + S(b, :)).^2;  w2 = (S(a, :) - 2 * S(b, :) + S(b - 1, :)).^2;
    switch rule
        case "tol",    c = max(cj / tol, (w1 + w2) / 2 ./ m2 / tol);
        case "budget", c = max(cj / (budget * (L - 1)), (w1 + w2) / 2 ./ m2 / (budget * (L - 2) / 2));
        case "detect", c = max(cj / q(L - 1), max(w1, w2) ./ m2 / q(L - 2));
    end
    c(~isfinite(c)) = Inf;
    [cMin, k] = min(c, [], 1);
    if all(cMin <= 1), i0 = a(k)';  i1 = b(k)';  break, end
end
if any(isnan(i0))
    error("shapeNoise:noSegment", "no length of %d or more samples has an end-matched segment in all %d channels (SegmentRule ""%s"").", ...
        minLength, C, rule);
end
L = i1(1) - i0(1) + 1;
info = struct("L", L, "rule", rule, "qJ", NaN, "qS", NaN, "cJump", NaN(1, C), "cSlip", NaN(1, C), "cSlipMax", NaN(1, C));
if rule == "detect", info.qJ = q(L - 1);  info.qS = q(L - 2); end
for ch = 1:C
    G = S(i0(ch):i1(ch), ch);  m2c = mean(diff(G, 2).^2);
    v1 = (G(2) - 2 * G(1) + G(end))^2;  v2 = (G(1) - 2 * G(end) + G(end - 1))^2;
    info.cJump(ch) = (G(1) - G(end))^2 / mean(diff(G).^2);
    info.cSlip(ch) = (v1 + v2) / 2 / m2c;  info.cSlipMax(ch) = max(v1, v2) / m2c;
end
end

function [u, k, ok, dA, nFail, msg] = drawEndMatched(ref, aRef, fs, band, tol, maxDraw)
% endMatchedNoise_v001 drawAxis, with estimateWhittleAlpha for alphaPWLW_v001.
L = numel(ref);  A = abs(fft(ref - mean(ref)));  h = floor((L - 1) / 2);  sR = std(ref);  ok = false;  dA = NaN;  nFail = 0;  msg = "";
for k = 1:maxDraw
    Z = zeros(L, 1);  Z(2:h + 1) = A(2:h + 1) .* exp(2i * pi * rand(h, 1));  Z(L:-1:L - h + 1) = conj(Z(2:h + 1));
    if mod(L, 2) == 0, Z(L/2 + 1) = A(L/2 + 1) * sign(randn); end
    u = real(ifft(Z));  u = u - mean(u);  u = u / std(u) * sR;
    d = NaN;
    try
        d = estimateWhittleAlpha(u, fs, band(1), band(2)) - aRef;
    catch err
        nFail = nFail + 1;  msg = string(err.message);   % rejected draw, reported by the caller
    end
    if abs(d) <= tol, ok = true;  dA = d;  return, end
end
u = NaN(L, 1);
end

function n = legacyShape(resid, fs, alpha, phaseSource)
% 1.x construction (shapeXu_local), unchanged.
M = numel(resid);
sig = std(resid, 0, 1);
if sig == 0
    n = zeros(M, 1);
    return
end
carrier = shapeNoiseCarrier(phaseSource, M, fs, alpha, sig);
if ~isreal(carrier) || ~iscolumn(carrier) || numel(carrier) ~= M
    error("shapeNoise:InvalidCarrier", ...
        "Phase carrier must return a real column vector of length %d.", M);
end
empiricalMagnitude = abs(fft(resid));
carrierSpectrum = fft(carrier);
shapedSpectrum = empiricalMagnitude .* exp(1i * angle(carrierSpectrum));
n = real(ifft(shapedSpectrum));
carrierStd = std(n, 0, 1);
if carrierStd > 0
    n = n * (sig / carrierStd);
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
