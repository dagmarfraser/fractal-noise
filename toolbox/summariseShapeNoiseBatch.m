function S = summariseShapeNoiseBatch(infos, options)
%SUMMARISESHAPENOISEBATCH Acceptance and tries across many shapeNoise calls, with a warning when acceptance is significantly low.
%
%   S = summariseShapeNoiseBatch(infos) takes the info outputs of many shapeNoise calls and reports the share of
%   channels with an accepted draw and the tries used. A single call has one or two channels, so a percentage only means
%   something across a batch.
%
%   THE STANDARD IS SET BY DESIGN, NOT FROM ANY DATASET. On exact references (k-fold integrated fractional Gaussian noise,
%   known exponent, alpha 0.5 to 4.5, fs 133 and 240, N 1000 and 3000) at least 99% of channels that pass the segment rule are
%   accepted within MaxDraw = 20. The batch is compared with that standard by an exact one-sided binomial test: the warning
%   "shapeNoise:lowAcceptance" is raised only when the number accepted is significantly below Standard at Level, and only
%   when there are at least MinChannels channels. Below that number the share is reported without a test. Low acceptance is
%   never an error: the function reports, and the decision stays with you.
%
% SYNTAX:
%   S = summariseShapeNoiseBatch(infos)
%   S = summariseShapeNoiseBatch(infos, Standard=0.99, Level=0.01, MinChannels=42)
%
% INPUTS:
%   infos       - Info outputs of shapeNoise (Method "endMatched"): a cell array of structs, or a struct array. Calls that
%                 ended in an error return no info and are not part of the batch: count those yourself.
%   Standard    - Name-value. Acceptance share the batch is tested against (0.99).
%   Level       - Name-value. One-sided significance level of the binomial test (0.01).
%   MinChannels - Name-value. Fewest channels for which the test is run (42: the smallest batch with at least 80% power to
%                 detect a true acceptance of 90% against the default Standard and Level).
%
% OUTPUT:
%   S - struct: nInfos, nChannels, nAccepted, acceptance (share of channels accepted), perDrawAcceptance (accepted channels
%       over all draws made: the estimated chance that one draw is accepted, allowing for channels that ran out of draws), triesMedian, triesP90 (tries per
%       channel, nearest-rank 90th percentile; a channel with no accepted draw counts its MaxDraw), pValue (probability of
%       nAccepted or fewer if the true acceptance were Standard; NaN when not tested), Standard, Level, MinChannels, and decision:
%       "ok" (tested, not significantly low), "low" (tested, significantly low, warned) or "tooFew" (not tested).
%
% NOTES:
%   Channels are treated as independent. Two axes of one movement share a trial, so with strongly related axes the test is
%   somewhat too ready to warn.
%
% See also: shapeNoise, estimateWhittleAlpha
%
% Toolbox: fractalnoise

arguments
    infos
    options.Standard (1, 1) double {mustBeGreaterThan(options.Standard, 0), mustBeLessThan(options.Standard, 1)} = 0.99
    options.Level (1, 1) double {mustBeGreaterThan(options.Level, 0), mustBeLessThan(options.Level, 1)} = 0.01
    options.MinChannels (1, 1) double {mustBeInteger, mustBePositive} = 42
end

[accepted, tries, nInfos] = collectChannels(infos);
nChannels = numel(accepted);
nAccepted = nnz(accepted);
sorted = sort(tries);
pValue = NaN;
decision = "tooFew";
tested = nChannels >= options.MinChannels;
if tested
    pValue = binomialLowerTail(nAccepted, nChannels, options.Standard);
    decision = "ok";
    isLow = pValue < options.Level;
    if isLow
        decision = "low";
        warning("shapeNoise:lowAcceptance", ...
            "%d of %d channels (%.1f%%) had an accepted draw, significantly below the design standard of %.1f%% (one-sided binomial p = %.2g). " + ...
            "On real data this usually means the residual is not a clean power law, or that AlphaTol is tight for this record length. " + ...
            "Inspect the residual's spectrum (info.dAlpha, info.Lfrac, a spectrum plot), then consider a larger MaxDraw or a looser AlphaTol, and report the acceptance share.", ...
            nAccepted, nChannels, 100 * nAccepted / nChannels, 100 * options.Standard, pValue);
    end
end
S = struct("nInfos", nInfos, "nChannels", nChannels, "nAccepted", nAccepted, "acceptance", nAccepted / nChannels, ...
    "perDrawAcceptance", nAccepted / sum(tries), "triesMedian", median(tries), "triesP90", sorted(ceil(0.9 * nChannels)), "pValue", pValue, ...
    "Standard", options.Standard, "Level", options.Level, "MinChannels", options.MinChannels, "decision", decision);
end

function [accepted, tries, nInfos] = collectChannels(infos)
% One entry per channel from every info struct. Anything that is not a shapeNoise end-matched info is an error.
if isstruct(infos), infos = num2cell(infos); end
if ~iscell(infos), error("summariseShapeNoiseBatch:badInfo", "infos must be a cell array of shapeNoise info structs, or a struct array."); end
nInfos = numel(infos);
if nInfos == 0, error("summariseShapeNoiseBatch:empty", "infos is empty: there is no batch to summarise."); end
accepted = false(0, 1);  tries = zeros(0, 1);
for k = 1:nInfos
    info = infos{k};
    if ~isstruct(info) || ~isscalar(info) || ~isfield(info, "accepted") || ~isfield(info, "nDraw")
        error("summariseShapeNoiseBatch:badInfo", "infos{%d} is not an end-matched shapeNoise info struct (fields accepted and nDraw are needed).", k);
    end
    a = info.accepted;  d = info.nDraw;
    if ~islogical(a) || isempty(a), error("summariseShapeNoiseBatch:badInfo", "infos{%d}.accepted must be a non-empty logical array.", k); end
    if numel(a) ~= numel(d), error("summariseShapeNoiseBatch:badInfo", "infos{%d}: accepted has %d channels but nDraw has %d.", k, numel(a), numel(d)); end
    if ~isnumeric(d) || any(~isfinite(d)) || any(d < 1) || any(d ~= round(d))
        error("summariseShapeNoiseBatch:badInfo", "infos{%d}.nDraw must hold positive whole numbers.", k);
    end
    accepted = [accepted; a(:)]; %#ok<AGROW>
    tries = [tries; double(d(:))]; %#ok<AGROW>
end
end

function p = binomialLowerTail(x, n, q)
% P(X <= x) for X ~ Binomial(n, q), summed on the log scale (base MATLAB only).
i = (0:x)';
logTerm = gammaln(n + 1) - gammaln(i + 1) - gammaln(n - i + 1) + i * log(q) + (n - i) * log1p(-q);
p = min(1, sum(exp(logTerm)));
end
