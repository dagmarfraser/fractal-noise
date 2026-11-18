function [T, E] = checkShapeNoiseEndMatched_v003(opts)
%CHECKSHAPENOISEENDMATCHED_V003 High-alpha validation of shapeNoise on break-free exact references (v001 and v002 kept).
%   [T, E] = checkShapeNoiseEndMatched_v003() runs the grid; T has one row per cell, E lists every error by cell, mode and identifier.
%
%   Grid: alpha 3.25:0.25:6.75 without 5 (no exact series at odd integers), the brackets 4.95 and 5.05 in its place, and the brackets
%   2.95 and 3.05 for continuity with v002; fs 120, 133, 240, 300 Hz crossed with N 1000, 3000; exact k-fold-integrated fGn only
%   (validation/exactReference; truth known; no low-frequency break); 200 replicates per cell; band 1-20 Hz.
%
%   PRIMARY SCORE: recovery of the known exponent by estimateWhittleAlpha on the references, against the band-effective alpha
%   (bandEffectiveAlpha_v001) and against nominal alpha.
%   SECONDARY: shapeNoise on the same reference (default detect rule, alpha = [], AlphaBand = [1 20]): surrogate alpha against truth, gap
%   against the segment's own alpha, L/N, d1 and d2 as in v002 (d1 and d2 are near-automatic, because the surrogate keeps the segment's
%   amplitude spectrum).
%
%   FAILURE BREAKDOWN: each replicate is classed as ok, noSegment, segmentAlpha, noAccepted or other, for two paired calls on the same
%   reference from the same random state: the default (CheckUncovered=true) and CheckUncovered=false. nRescued counts replicates that
%   fail by default and succeed with the check off; nNewFail the reverse; nSameOutcome those with the same class and the same surrogate.
%   Scores use the default call. The segmentAlpha cause identifier (why the segment could not be estimated) is kept in E.
%
%   Seeds: rng(1e9*Master + 1e8*config + 1e7*1 + 1e4*alphaIndex + replicate, 'twister') for every replicate, Master = 3, so results do
%   not depend on the number of workers and do not reuse the v002 streams. Errors are counted, never hidden; a failed replicate is NaN.
arguments
    opts.Alphas (1, :) double = sort([2.95, 3.05, 3.25:0.25:4.75, 4.95, 5.05, 5.25:0.25:6.75])
    opts.Config (:, 2) double = [120 1000; 120 3000; 133 1000; 133 3000; 240 1000; 240 3000; 300 1000; 300 3000]
    opts.NRep (1, 1) double {mustBeInteger, mustBePositive} = 200
    opts.Master (1, 1) double {mustBeInteger, mustBePositive} = 3
    opts.Band (1, 2) double = [1 20]
end
jobs = struct('fs', {}, 'N', {}, 'alpha', {}, 'aEff', {}, 'ci', {}, 'ai', {});
for ci = 1:size(opts.Config, 1)
    fs = opts.Config(ci, 1);  N = opts.Config(ci, 2);
    for ai = 1:numel(opts.Alphas)
        a = opts.Alphas(ai);
        [~, info] = exactRefSeries_v001(8, a, 1);   % errors loudly at an odd-integer alpha
        jobs(end + 1) = struct('fs', fs, 'N', N, 'alpha', a, 'aEff', bandEffectiveAlpha_v001(info.S, fs, opts.Band, N), 'ci', ci, 'ai', ai); %#ok<AGROW>
    end
end
nJob = numel(jobs);  rows = cell(nJob, 1);  errs = cell(nJob, 1);
nRep = opts.NRep;  master = opts.Master;  band = opts.Band;
parfor j = 1:nJob
    [rows{j}, errs{j}] = runCell(jobs(j), nRep, master, band);
end
T = vertcat(rows{:});
E = vertcat(errs{:});
end

function [row, E] = runCell(c, nRep, master, band)
warning('off', 'shapeNoise:noAcceptedDraw');  warning('off', 'shapeNoise:estimateFailed');
aRef = NaN(nRep, 1);  aSur = aRef;  gap = aRef;  d1 = aRef;  d2 = aRef;  LN = aRef;  nDraw = aRef;
clsD = strings(nRep, 1);  clsF = strings(nRep, 1);  same = false(nRep, 1);
ids = strings(0, 1);  estFail = 0;
for r = 1:nRep
    rng(1e9 * master + 1e8 * c.ci + 1e7 + 1e4 * c.ai + r, 'twister');
    try
        ref = exactRefSeries_v001(c.N, c.alpha, 1);  ref = ref(:);
    catch err, ids(end + 1, 1) = "reference:" + err.identifier;  clsD(r) = "noReference";  clsF(r) = "noReference";  continue %#ok<AGROW>
    end
    try, aRef(r) = estimateWhittleAlpha(ref, c.fs, band(1), band(2));
    catch err, ids(end + 1, 1) = "estimateRef:" + err.identifier; %#ok<AGROW>
    end
    state = rng;
    D = callShapeNoise(ref, c.fs, band, true);
    rng(state);
    F = callShapeNoise(ref, c.fs, band, false);
    clsD(r) = D.cls;  clsF(r) = F.cls;  same(r) = D.cls == F.cls && isequaln(D.s, F.s);
    if D.id ~= "", ids(end + 1, 1) = "default:" + D.id; end %#ok<AGROW>
    if F.id ~= "", ids(end + 1, 1) = "checkOff:" + F.id; end %#ok<AGROW>
    if isempty(D.info), continue, end
    estFail = estFail + sum(D.info.nEstFail);  LN(r) = D.info.Lfrac;  nDraw(r) = D.info.nDraw;
    if D.cls ~= "ok", continue, end
    seg = ref(D.info.i0:D.info.i1);
    try
        aSur(r) = estimateWhittleAlpha(D.s, c.fs, band(1), band(2));
        gap(r) = aSur(r) - D.info.alphaTarget;
        d1(r) = std(diff(D.s)) / std(diff(seg));  d2(r) = std(diff(D.s, 2)) / std(diff(seg, 2));
    catch err, ids(end + 1, 1) = "estimateSurrogate:" + err.identifier; %#ok<AGROW>
    end
end
ok = ~isnan(aRef);  n = nnz(ok);  truth = c.aEff;
failD = clsD ~= "ok";  failF = clsF ~= "ok";
row = table(c.fs, c.N, c.alpha, c.aEff, nRep, n, ...
    mean(aRef, 'omitnan') - c.alpha, mean(aRef, 'omitnan') - truth, std(aRef, 'omitnan'), std(aRef, 'omitnan') / sqrt(max(n, 1)), ...
    sqrt(mean((aRef - truth) .^ 2, 'omitnan')), mean(abs(aRef(ok) - truth) <= 0.25), ...
    mean(aSur, 'omitnan') - c.alpha, mean(aSur, 'omitnan') - truth, median(gap, 'omitnan'), prctile(gap, 10), prctile(gap, 90), mean(gap(~isnan(gap)) < -0.5), ...
    median(LN, 'omitnan'), min(LN), median(d1, 'omitnan'), median(d2, 'omitnan'), median(nDraw, 'omitnan'), ...
    nnz(startsWith(ids, "reference:")), nnz(startsWith(ids, "estimateRef:")), nnz(startsWith(ids, "estimateSurrogate:")), estFail, ...
    nnz(clsD == "ok"), nnz(clsD == "noSegment"), nnz(clsD == "segmentAlpha"), nnz(clsD == "noAccepted"), nnz(clsD == "other"), ...
    nnz(clsF == "ok"), nnz(clsF == "noSegment"), nnz(clsF == "segmentAlpha"), nnz(clsF == "noAccepted"), nnz(clsF == "other"), ...
    nnz(failD & ~failF), nnz(~failD & failF), nnz(same), ...
    'VariableNames', {'fs', 'N', 'alpha', 'alphaEffective', 'nRep', 'nEstimated', ...
    'biasNominal', 'biasEffective', 'sdAlpha', 'mcseBias', 'rmseEffective', 'shareWithin025Effective', ...
    'surrBiasNominal', 'surrBiasEffective', 'gapSelfMedian', 'gapSelfP10', 'gapSelfP90', 'gapSelfShareBelowMinus05', ...
    'LfracMedian', 'LfracMin', 'd1Median', 'd2Median', 'nDrawMedian', ...
    'errReference', 'errEstimateRef', 'errEstimateSurrogate', 'nEstimateFailedDraws', ...
    'okDefault', 'noSegmentDefault', 'segmentAlphaDefault', 'noAcceptedDefault', 'otherDefault', ...
    'okCheckOff', 'noSegmentCheckOff', 'segmentAlphaCheckOff', 'noAcceptedCheckOff', 'otherCheckOff', ...
    'nRescued', 'nNewFail', 'nSameOutcome'});
[u, ~, g] = unique(ids);
if isempty(u)
    E = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), 'VariableNames', {'fs', 'N', 'alpha', 'identifier', 'count'});
else
    m = numel(u);
    E = table(repmat(c.fs, m, 1), repmat(c.N, m, 1), repmat(c.alpha, m, 1), u(:), accumarray(g(:), 1), 'VariableNames', {'fs', 'N', 'alpha', 'identifier', 'count'});
end
end

function R = callShapeNoise(ref, fs, band, checkUncovered)
% One shapeNoise call, classed. A failure is returned as data, with its identifier and, for segmentAlpha, the identifier of its cause.
R = struct('cls', "ok", 'id', "", 's', [], 'info', []);
try
    [R.s, R.info] = shapeNoise(ref, fs, [], AlphaBand=band, CheckUncovered=checkUncovered);
    if ~all(R.info.accepted), R.cls = "noAccepted";  R.s = []; end
catch err
    R.id = string(err.identifier);
    switch R.id
        case "shapeNoise:noSegment",    R.cls = "noSegment";
        case "shapeNoise:segmentAlpha", R.cls = "segmentAlpha";
            if ~isempty(err.cause), R.id = R.id + "<" + string(err.cause{1}.identifier); end
        otherwise,                      R.cls = "other";
    end
end
end
