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
%   Validated on real template-subtracted hand-movement residuals from
%   several datasets: filtered velocity
%   and acceleration noise 1.00-1.02 of the residual's (median), with alpha
%   held by the selection step. The 1.x construction failed that check
%   (1-27% of surrogates more than 0.5 too shallow at alpha 3-6; see
%   validation/checkShapeNoiseEndMatched_v001.m and
%   validation/checkShapeNoiseMechanism_v001.m).
%
%   RANGE, measured against known truth on exact references (k-fold integrated fractional Gaussian noise, no low-frequency break;
%   fs 120, 133, 240, 300 Hz; N 1000 and 3000; band 1-20 Hz; 200 series per cell; default MinLength 120). Runs, md5s and tables:
%   validation/rangeEvidence_v002.md (v001 holds the same runs with MinLength 32).
%   A failure is an error or no accepted draw. Alpha is the generator's exponent, because failure follows it and not the in-band slope:
%   alpha 3.05 has an in-band slope of only 2.35-2.61 and still fails like alpha 2.95.
%     alpha -0.75 to 2.5   validated: at most 1.3% of series fail (3.5% in the worst configuration)
%     alpha 2.75           8% fail (21% in the worst configuration)
%     alpha 2.95 to 4      25-51% fail
%     alpha 4.25 to 4.5    69-81% fail
%     alpha 4.75           92% fail
%     alpha 4.95 to 6.75   96-100% fail
%   These failures are loud: "shapeNoise:noSegment" (with MinLength 120 a segment too short to estimate is no longer searched),
%   or "shapeNoise:segmentAlpha" if the segment's exponent still cannot be estimated. Nothing is
%   shortened further or substituted. Odd-integer alpha has no exact reference, so it has no known-truth check.
%   Steeper series work when they have a low-frequency break or leak. noiseXu with Phi 0.99 is accepted in at least 98% of series up to an
%   in-band slope of 4.5, 92% at 4.5-5, 81% at 5-5.5, 70% at 5.5-6 and 43% at 6-7. Real movement residuals after stripping K = 1 harmonic:
%   97.0% of 7948 axes accepted, 93.7% of those with a residual alpha above 5. Against MinLength 32 no axis moved between accepted and
%   failed; 5 failures changed kind (3 now noSegment, 2 records of 101 samples now "shapeNoise:short").
%   The exponent itself is limited by estimateWhittleAlpha: within 0.08 of the in-band slope for alpha -1.5 to 9.5, and an error when the
%   estimate reaches AlphaBounds [-2 10] (45% of series at alpha -2, 33% at alpha 10).
%   Random phases keep the amplitude spectrum only: non-Gaussian structure in the residual is not reproduced.
%   Across many calls, summariseShapeNoiseBatch reports acceptance and tries, which show departures from a clean power law. Mouse data
%   recorded with operating-system pointer smoothing on are an example of what its warning is for.
%
%   TEMPLATE-SUBTRACTED DATA: pass only the part of the residual the
%   template actually covered. A stretch left at exactly zero (an uncovered
%   template tail) is an error here, because the segment search would
%   otherwise take it and feed a step into the surrogate.
%
% SYNTAX:
%   [n, info] = shapeNoise(resid, fs, [])          % target alpha = segment's own (recommended)
%   [n, info] = shapeNoise(resid, fs, alpha)       % explicit target, scalar or 1-by-C
%   [n, info] = shapeNoise(___, AlphaTol=0.3, MaxDraw=20, MinLength=120, AlphaBand=[1 20])   % band in the units of fs; "kinematic" = [1, min(20, fs/2 - 1)] Hz
%   [n, info] = shapeNoise(___, CheckUncovered=false)                      % allow exact zeros at the ends (counts, rainfall)
%   [n, info] = shapeNoise(___, SegmentRule="budget", SegmentBudget=0.01)   % wrap adds <= 1% derivative energy
%   [n, info] = shapeNoise(___, SegmentRule="tol", FairShareTol=1)          % wrap <= Tol x mean terms
%   n = shapeNoise(resid, fs, alpha, Method="legacy")                      % 1.x construction; warns
%   n = shapeNoise(resid, fs, alpha, Method="legacy", PhaseSource="white") % legacy carriers
%
% INPUTS:
%   resid - Residual, real N-by-C matrix (one column per channel, e.g. the
%           two axes of a movement) or a vector. "legacy" takes a vector.
%   fs    - Sampling rate in cycles per unit of time: Hz for movement, per year for monthly data (fs = 12), per day for
%           daily data. Every frequency argument (AlphaBand) is in the same unit.
%   alpha - Target exponent for alpha selection: [] (default use: each
%           channel's segment alpha by estimateWhittleAlpha over AlphaBand),
%           a scalar, or 1-by-C. Under "legacy" the carrier's alpha.
%   Method       - "endMatched" (default) or "legacy".
%   AlphaTol     - Accept a draw when |alpha(draw) - target| <= AlphaTol (0.3). On exact references 0.3 is about 2 standard
%                  deviations of estimateWhittleAlpha over 1-20 Hz (1.1 to 3.5 across fs 120-300 Hz and N 1000-3000). That standard deviation is measured; its size agrees with the Hann-taper variance factor of Velasco (1999), within the limits given in the estimateWhittleAlpha help. If draws
%                  centre on the target, 20 draws then accept at least 99% of channels for any standard deviation up to 1.1, so
%                  repeated failure points to surrogates centred away from the target. Across many calls, see summariseShapeNoiseBatch.
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
%   MinLength    - Shortest segment searched (120 samples; 32 before 2026-11-21). No segment at
%                  or above it is an error ("shapeNoise:noSegment"); a record shorter than it is an error ("shapeNoise:short").
%                  Below about 120 samples a segment is often too short for estimateWhittleAlpha (its DFA pre-estimate needs
%                  100 samples after differencing). On the exact-reference grids, 120 removed every segmentAlpha failure caused
%                  by a too-short segment; those series now end in noSegment, and failure rates rose by at most 2.6 percentage
%                  points (validation/rangeEvidence_v002.md). Lower it only if you check the exponent yourself.
%   AlphaBand    - [fLow fHigh] for estimateWhittleAlpha, in the units of fs (0 < fLow < fHigh < fs/2), or "kinematic"
%                  ([1, min(20, fs/2 - 1)] Hz, the band validated on hand-movement residuals at fs 240; below fs 160 Hz its upper edge
%                  exceeds 25% of Nyquist, so choose a lower one). NOT GIVEN: the kinematic
%                  band is used with the warning "shapeNoise:defaultAlphaBand", because it assumes Hz-scale movement
%                  data; for any other data pass the band you mean. An invalid band is an error ("shapeNoise:alphaBand").
%   CheckUncovered - true (default): exact-zero stretches of 3 or more samples at either end are an error
%                  (an uncovered template stretch, "shapeNoise:zeroTail"). false: skip the check, for data with
%                  legitimate zeros (counts, rainfall); you then own the check that no template tail is left at zero.
%   PeakCheck    - true (default): each channel's segment goes through the estimateWhittleAlpha peak check once, which warns
%                  "estimateWhittleAlpha:peakInBand" when a narrow peak in the band would bias the exponent this function
%                  selects on (a harmonic left after the template: raise K). info.peakInBand and info.peakHz record it either way;
%                  false only silences the warning.
%   PhaseSource  - "legacy" only: "xu" (default), "white", "pinknoise",
%                  "dsp" or a function handle @(M, fs, alpha).
%
% OUTPUTS:
%   n    - L-by-C surrogate ("endMatched"), or numel(resid)-by-1 ("legacy").
%   info - "endMatched": AlphaBand (the [fLow fHigh] actually used), SegmentRule, i0, i1 (1-by-C segment start and end
%          in resid), L, Lfrac (L/N), cJump, cSlip, cSlipMax (squared wrap
%          step, mean and larger wrapped second difference, over the
%          segment's mean terms; 1-by-C), qJ, qS (the "detect" levels at L),
%          alphaTarget, nDraw, accepted, dAlpha, nEstFail (draws whose alpha
%          estimate errored, counted as rejected and warned about as
%          "shapeNoise:estimateFailed"; 1-by-C), peakInBand, peakHz (1-by-C: the
%          estimateWhittleAlpha peak check on each channel's segment, which warns
%          "estimateWhittleAlpha:peakInBand" once per channel; PeakCheck=false turns it
%          off). "legacy": Method only.
%
% NOTES ON SCOPE (read before using it outside movement data):
%   VALIDATED RANGE. See RANGE above. Known-truth checks cover fs 120-300 Hz and N 1000 and 3000 with the 1-20 Hz band. Long records
%     (N of 10^4 or more), other bands and non-Hz units have no known-truth check; look at the surrogate against your data
%     (info.dAlpha, info.Lfrac, a spectrum plot) before relying on it.
%   CHANNELS. Each channel gets its own independent random phases, so the surrogate keeps each channel's amplitude spectrum and
%     length L but NOT any coherence or correlation between channels. That is intended for the x and y noise of a movement; it is
%     not a multivariate surrogate.
%   COST. The segment search tries lengths from N downwards and stops at the first length that works, with vector work of size
%     N - L + 1 per length: fast when the full record nearly matches, up to O(N^2) when only short segments match (steep alpha,
%     large N). Long, steep records can be slow.
%
% LEGACY ("legacy", the 1.x behaviour, bit for bit):
%   F = |FFT(resid)| .* exp(1i*angle(FFT(carrier))), n = real(ifft(F))
%   rescaled to std(resid). The carrier's phases decide whether the
%   surrogate stays steep; at alpha 3-6 a tail of surrogates comes out
%   shallow. Kept only to reproduce earlier results; warns on every call
%   ("shapeNoise:legacyConstruction").
%
% See also: estimateWhittleAlpha, noiseXu, estimateIRASA, summariseShapeNoiseBatch
%
% Toolbox: fractalnoise

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
    options.MinLength (1, 1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.MinLength, 5)} = 120
    options.AlphaBand = [NaN NaN]
    options.CheckUncovered (1, 1) logical = true
    options.PeakCheck (1, 1) logical = true
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
band = resolveAlphaBand(options.AlphaBand, fs);
if options.CheckUncovered, checkZeroTail(resid); end
if any(std(resid, 0, 1) == 0), error("shapeNoise:constant", "a channel of resid is constant."); end

if options.SegmentRule == "budget" && ~(isscalar(options.SegmentBudget) && options.SegmentBudget > 0 && options.SegmentBudget < 1)
    error("shapeNoise:segmentBudget", "SegmentRule=""budget"" needs a SegmentBudget in (0, 1).");
end
[i0, i1, seg] = endMatchedSegments(resid, options.SegmentRule, options.FairShareTol, options.SegmentBudget, options.MinLength);
L = seg.L;
peakInBand = false(1, C);  peakHz = NaN(1, C);   % the estimateWhittleAlpha peak check, run once per channel on its segment
if isempty(alpha)
    target = NaN(1, C);
    for ch = 1:C
        try
            [target(ch), ~, segInfo] = estimateWhittleAlpha(resid(i0(ch):i1(ch), ch), fs, band(1), band(2), PeakCheck=options.PeakCheck);
            [peakInBand(ch), peakHz(ch)] = peakFlag(segInfo);
        catch cause
            err = MException("shapeNoise:segmentAlpha", ...
                "channel %d: the exponent of the end-matched segment (L = %d of %d samples) could not be estimated: %s", ...
                ch, L, N, cause.message);
            throw(addCause(err, cause));
        end
    end
else
    for ch = 1:C   % a given target: the segment is estimated only for the peak check; a failure here is left to the draws
        try
            [~, ~, segInfo] = estimateWhittleAlpha(resid(i0(ch):i1(ch), ch), fs, band(1), band(2), PeakCheck=options.PeakCheck);
            [peakInBand(ch), peakHz(ch)] = peakFlag(segInfo);
        catch
        end
    end
end
if isempty(alpha)
    % target estimated from the segment above
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
        "No draw within AlphaTol = %.2f of the target in MaxDraw = %d tries for channel(s) %s; returned as NaN.", ...
        options.AlphaTol, options.MaxDraw, mat2str(find(~accepted)));
end
info = struct("Method", "endMatched", "SegmentRule", seg.rule, "i0", i0, "i1", i1, "L", L, "Lfrac", L / N, ...
    "cJump", seg.cJump, "cSlip", seg.cSlip, "cSlipMax", seg.cSlipMax, "qJ", seg.qJ, "qS", seg.qS, "alphaTarget", target, "nDraw", nDraw, "accepted", accepted, "dAlpha", dAlpha, "nEstFail", nEstFail, "AlphaBand", band, ...
    "peakInBand", peakInBand, "peakHz", peakHz);
if any(nEstFail > 0)
    warning("shapeNoise:estimateFailed", ...
        "estimateWhittleAlpha failed on %d draw(s) (counted as rejected; info.nEstFail). Last message: %s", ...
        sum(nEstFail), estMsg);
end
end

function band = resolveAlphaBand(ab, fs)
% AlphaBand: [fLow fHigh] in the units of fs, or "kinematic" (= [1, min(20, fs/2 - 1)] Hz, the band validated on hand movement),
% or not supplied: the kinematic band WITH A WARNING, because it assumes Hz-scale movement data.
kin = @() [1, min(20, fs/2 - 1)];
if isstring(ab) || ischar(ab)
    if string(ab) ~= "kinematic"
        error("shapeNoise:alphaBand", "AlphaBand must be a [fLow fHigh] pair or ""kinematic"", got ""%s"".", string(ab));
    end
    band = kin();
elseif isnumeric(ab) && numel(ab) == 2 && all(isnan(ab))
    band = kin();
    warning("shapeNoise:defaultAlphaBand", ...
        "AlphaBand not given: using the kinematic default [%g %g] (the band validated on hand-movement residuals, fs in Hz). " + ...
        "For other data pass AlphaBand=[fLow fHigh] in the units of fs, or AlphaBand=""kinematic"" to silence this warning.", band(1), band(2));
elseif isnumeric(ab) && numel(ab) == 2 && all(isfinite(ab)) && ab(1) > 0 && ab(2) > ab(1)
    band = double(ab(:)');
else
    error("shapeNoise:alphaBand", "AlphaBand must be [fLow fHigh] with 0 < fLow < fHigh, or ""kinematic"".");
end
if band(2) >= fs/2 || band(2) <= band(1)
    error("shapeNoise:alphaBand", ...
        "The alpha band [%g %g] is invalid for fs = %g (needs fLow < fHigh < fs/2). The kinematic default assumes fs of several tens of Hz; pass AlphaBand for other data.", ...
        band(1), band(2), fs);
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
% Rules "detect", "budget", "tol"; chi2inv(1 - p, 1) written as 2*erfcinv(p)^2 (base MATLAB).
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

function [flag, hz] = peakFlag(info)
% The estimateWhittleAlpha peak check's verdict on one segment.
flag = info.peakRatio > info.peakThreshold;  hz = NaN;
if flag, hz = info.peakHz; end
end

function [u, k, ok, dA, nFail, msg] = drawEndMatched(ref, aRef, fs, band, tol, maxDraw)
% One axis: random-phase draws from the segment, kept when estimateWhittleAlpha matches.
L = numel(ref);  A = abs(fft(ref - mean(ref)));  h = floor((L - 1) / 2);  sR = std(ref);  ok = false;  dA = NaN;  nFail = 0;  msg = "";
for k = 1:maxDraw
    Z = zeros(L, 1);  Z(2:h + 1) = A(2:h + 1) .* exp(2i * pi * rand(h, 1));  Z(L:-1:L - h + 1) = conj(Z(2:h + 1));
    if mod(L, 2) == 0, Z(L/2 + 1) = A(L/2 + 1) * sign(randn); end
    u = real(ifft(Z));  u = u - mean(u);  u = u / std(u) * sR;
    d = NaN;
    try
        d = estimateWhittleAlpha(u, fs, band(1), band(2), PeakCheck=false) - aRef;   % checked once, on the segment
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
