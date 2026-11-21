function [T, E] = checkShapeNoiseEndMatched_v005(opts)
%CHECKSHAPENOISEENDMATCHED_V005 The v002 low-alpha grid, identical except for one added option, ShapeArgs (v002 kept, unedited).
%   ShapeArgs = {} reproduces v002 exactly (same grid, seeds and scores); ShapeArgs = {"MinLength", 120} reruns it with that shapeNoise option.
%   [T, E] = checkShapeNoiseEndMatched_v002() runs the grid; T has one row per cell, E lists every error by cell and identifier.
%
%   Grid: alpha -1:0.25:3 plus the brackets 0.95, 1.05, 2.95, 3.05; fs 120, 133, 240, 300 Hz crossed with N 1000, 3000;
%   three reference generators: exact k-fold-integrated fGn (validation/exactReference; truth known; no series at alpha -1, 1, 3,
%   where no exact representative exists), noiseXu Phi 0.99 (default) and noiseXu Phi 1; 200 replicates per cell; band 1-20 Hz.
%
%   PRIMARY SCORE: recovery of the known exponent by estimateWhittleAlpha on the exact references, against the band-effective
%   alpha (bandEffectiveAlpha_v001: the in-band slope of the analytic spectrum, the estimator's stated target) and against nominal alpha.
%   SECONDARY: shapeNoise on the same reference (default detect rule, alpha = [], AlphaBand = [1 20]): surrogate alpha against truth,
%   gap against the segment's own alpha (guaranteed by the selection step, so a check of the code), L/N, and d1, d2 = std(diff(surrogate, k)) /
%   std(diff(segment, k)). d1 and d2 are near-automatic when L = N, because the surrogate keeps the amplitude spectrum.
%
%   Seeds: rng(1e9*Master + 1e8*config + 1e7*generator + 1e4*alphaIndex + replicate, 'twister') for every replicate, so results do not
%   depend on the number of workers. Errors are counted per cell by stage and identifier, never hidden; a failed replicate is NaN.
arguments
    opts.Alphas (1, :) double = sort([-1:0.25:3, 0.95, 1.05, 2.95, 3.05])
    opts.Config (:, 2) double = [120 1000; 120 3000; 133 1000; 133 3000; 240 1000; 240 3000; 300 1000; 300 3000]
    opts.NRep (1, 1) double {mustBeInteger, mustBePositive} = 200
    opts.Master (1, 1) double {mustBeInteger, mustBePositive} = 1
    opts.Band (1, 2) double = [1 20]
    opts.ShapeArgs (1, :) cell = {}   % added in this version: further shapeNoise name-value pairs, e.g. {"MinLength", 120}
end
gens = ["exact", "noiseXu_Phi0.99", "noiseXu_Phi1"];
jobs = struct('fs', {}, 'N', {}, 'gen', {}, 'alpha', {}, 'aEff', {}, 'ci', {}, 'gi', {}, 'ai', {});
for ci = 1:size(opts.Config, 1)
    fs = opts.Config(ci, 1);  N = opts.Config(ci, 2);
    for gi = 1:numel(gens)
        for ai = 1:numel(opts.Alphas)
            a = opts.Alphas(ai);  aEff = NaN;
            if gens(gi) == "exact"
                if abs(mod(a + 1, 2)) < 1e-9 || abs(mod(a + 1, 2) - 2) < 1e-9, continue, end   % alpha -1, 1, 3: no exact representative
                [~, info] = exactRefSeries_v001(8, a, 1);
                aEff = bandEffectiveAlpha_v001(info.S, fs, opts.Band, N);
            end
            jobs(end + 1) = struct('fs', fs, 'N', N, 'gen', gens(gi), 'alpha', a, 'aEff', aEff, 'ci', ci, 'gi', gi, 'ai', ai); %#ok<AGROW>
        end
    end
end
nJob = numel(jobs);  rows = cell(nJob, 1);  errs = cell(nJob, 1);
nRep = opts.NRep;  master = opts.Master;  band = opts.Band;  shapeArgs = opts.ShapeArgs;
parfor j = 1:nJob
    [rows{j}, errs{j}] = runCell(jobs(j), nRep, master, band, shapeArgs);
end
T = vertcat(rows{:});
E = vertcat(errs{:});
end

function [row, E] = runCell(c, nRep, master, band, shapeArgs)
warning('off', 'shapeNoise:noAcceptedDraw');  warning('off', 'shapeNoise:estimateFailed');
aRef = NaN(nRep, 1);  aSur = aRef;  gap = aRef;  d1 = aRef;  d2 = aRef;  LN = aRef;  nDraw = aRef;
ids = strings(0, 1);  noAcc = 0;  estFail = 0;
for r = 1:nRep
    rng(1e9 * master + 1e8 * c.ci + 1e7 * c.gi + 1e4 * c.ai + r, 'twister');
    try
        switch c.gen
            case "exact",           ref = exactRefSeries_v001(c.N, c.alpha, 1);
            case "noiseXu_Phi0.99", ref = noiseXu(c.N, c.alpha, 1, c.fs);
            case "noiseXu_Phi1",    ref = noiseXu(c.N, c.alpha, 1, c.fs, Phi=1);
        end
        ref = ref(:);
    catch err, ids(end + 1, 1) = "reference:" + err.identifier;  continue %#ok<AGROW>
    end
    try, aRef(r) = estimateWhittleAlpha(ref, c.fs, band(1), band(2));
    catch err, ids(end + 1, 1) = "estimateRef:" + err.identifier; %#ok<AGROW>
    end
    try, [s, info] = shapeNoise(ref, c.fs, [], "AlphaBand", band, shapeArgs{:});
    catch err, ids(end + 1, 1) = "shapeNoise:" + err.identifier;  continue %#ok<AGROW>
    end
    estFail = estFail + sum(info.nEstFail);  LN(r) = info.Lfrac;  nDraw(r) = info.nDraw;
    if ~info.accepted, noAcc = noAcc + 1;  continue, end
    seg = ref(info.i0:info.i1);
    try
        aSur(r) = estimateWhittleAlpha(s, c.fs, band(1), band(2));
        gap(r) = aSur(r) - info.alphaTarget;
        d1(r) = std(diff(s)) / std(diff(seg));  d2(r) = std(diff(s, 2)) / std(diff(seg, 2));
    catch err, ids(end + 1, 1) = "estimateSurrogate:" + err.identifier; %#ok<AGROW>
    end
end
ok = ~isnan(aRef);  n = nnz(ok);  truth = c.aEff;
row = table(c.fs, c.N, c.gen, c.alpha, c.aEff, nRep, n, ...
    mean(aRef, 'omitnan') - c.alpha, mean(aRef, 'omitnan') - truth, std(aRef, 'omitnan'), std(aRef, 'omitnan') / sqrt(max(n, 1)), ...
    sqrt(mean((aRef - truth) .^ 2, 'omitnan')), mean(abs(aRef(ok) - truth) <= 0.25), ...
    mean(aSur, 'omitnan') - c.alpha, mean(aSur, 'omitnan') - truth, median(gap, 'omitnan'), prctile(gap, 10), prctile(gap, 90), mean(gap(~isnan(gap)) < -0.5), ...
    median(LN, 'omitnan'), min(LN), median(d1, 'omitnan'), median(d2, 'omitnan'), median(nDraw, 'omitnan'), ...
    nnz(startsWith(ids, "reference:")), nnz(startsWith(ids, "estimateRef:")), nnz(startsWith(ids, "shapeNoise:")), nnz(startsWith(ids, "estimateSurrogate:")), noAcc, estFail, ...
    'VariableNames', {'fs', 'N', 'generator', 'alpha', 'alphaEffective', 'nRep', 'nEstimated', ...
    'biasNominal', 'biasEffective', 'sdAlpha', 'mcseBias', 'rmseEffective', 'shareWithin025Effective', ...
    'surrBiasNominal', 'surrBiasEffective', 'gapSelfMedian', 'gapSelfP10', 'gapSelfP90', 'gapSelfShareBelowMinus05', ...
    'LfracMedian', 'LfracMin', 'd1Median', 'd2Median', 'nDrawMedian', ...
    'errReference', 'errEstimateRef', 'errShapeNoise', 'errEstimateSurrogate', 'nNoAcceptedDraw', 'nEstimateFailedDraws'});
[u, ~, g] = unique(ids);
if isempty(u)
    E = table(zeros(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), 'VariableNames', {'fs', 'N', 'generator', 'alpha', 'identifier', 'count'});
else
    m = numel(u);
    E = table(repmat(c.fs, m, 1), repmat(c.N, m, 1), repmat(c.gen, m, 1), repmat(c.alpha, m, 1), u(:), accumarray(g(:), 1), 'VariableNames', {'fs', 'N', 'generator', 'alpha', 'identifier', 'count'});
end
end
