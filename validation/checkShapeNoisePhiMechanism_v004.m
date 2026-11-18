function [T, E] = checkShapeNoisePhiMechanism_v004(opts)
%CHECKSHAPENOISEPHIMECHANISM_V004 Why shapeNoise fails near alpha 3 on break-free references, and whether Phi is the reason it does not on noiseXu.
%   [T, E] = checkShapeNoisePhiMechanism_v004() runs the grid; T has one row per cell and manipulation, E lists every failure identifier.
%
%   Question: the end-matched segment search needs two points of the series that meet. A series dominated by a slow trend rarely returns to
%   an earlier level, so long segments with matching ends should not exist. Phi < 1 makes noiseXu mean-reverting below its break.
%   Generators: exact (not at alpha 3) and noiseXu at Phi 1, 0.999, 0.995, 0.99, 0.95. All noiseXu rows of a replicate share one white-noise
%   seed, so Phi is the only difference between them (dose-response, paired).
%   Manipulations of the same series before the default shapeNoise call: raw, and with a least-squares polynomial of order 1 or 2 removed.
%   Diagnostics of the raw series: trendR2 (variance share of the linear trend), endJumpRatio ((x(1) - x(N))^2 over the mean squared step,
%   over the detect level q(N - 1)), and meanCrossings (crossings of the series mean).
%
%   Seeds: rng(1e9*Master + 1e8*config + 1e7*stream + 1e4*alphaIndex + replicate, 'twister'), Master = 2, stream 4 for exact and 5 for noiseXu
%   (checkShapeNoiseLayers_v004 uses streams 1 to 3). Errors are counted, never hidden.
arguments
    opts.Alphas (1, :) double = [2.95 3 3.05]
    opts.Config (:, 2) double = [120 1000; 120 3000; 133 1000; 133 3000; 240 1000; 240 3000; 300 1000; 300 3000]
    opts.Phis (1, :) double = [1 0.999 0.995 0.99 0.95]
    opts.NRep (1, 1) double {mustBeInteger, mustBePositive} = 200
    opts.Master (1, 1) double {mustBeInteger, mustBePositive} = 2
    opts.Band (1, 2) double = [1 20]
end
jobs = struct('fs', {}, 'N', {}, 'phi', {}, 'alpha', {}, 'ci', {}, 'ai', {});
for ci = 1:size(opts.Config, 1)
    for ai = 1:numel(opts.Alphas)
        for phi = [NaN, opts.Phis]   % NaN = exact
            if isnan(phi) && abs(mod(opts.Alphas(ai) + 1, 2)) < 1e-9, continue, end
            jobs(end + 1) = struct('fs', opts.Config(ci, 1), 'N', opts.Config(ci, 2), 'phi', phi, 'alpha', opts.Alphas(ai), 'ci', ci, 'ai', ai); %#ok<AGROW>
        end
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
manip = ["raw", "detrend1", "detrend2"];  orders = [NaN 1 2];  nM = numel(manip);
if isnan(c.phi), gen = "exact";  stream = 4; else, gen = "noiseXu_Phi" + string(c.phi);  stream = 5; end
cls = strings(nRep, nM);  LN = NaN(nRep, nM);  trendR2 = NaN(nRep, 1);  jump = trendR2;  cross = trendR2;  ids = strings(0, 1);
t = linspace(-1, 1, c.N)';
for r = 1:nRep
    rng(1e9 * master + 1e8 * c.ci + 1e7 * stream + 1e4 * c.ai + r, 'twister');
    try
        if isnan(c.phi), x = exactRefExtended_v004(c.N, c.alpha); else, x = noiseXu(c.N, c.alpha, 1, c.fs, Phi=c.phi); end
        x = x(:);
    catch err, ids(end + 1, 1) = "reference:" + err.identifier;  cls(r, :) = "noReference";  continue %#ok<AGROW>
    end
    trendR2(r) = 1 - var(removePoly(x, t, 1)) / var(x);
    jump(r) = (x(1) - x(end)) ^ 2 / mean(diff(x) .^ 2) / (2 * erfcinv(1 / (c.N - 1)) ^ 2);
    cross(r) = nnz(diff(sign(x - mean(x))) ~= 0);
    state = rng;
    for m = 1:nM
        rng(state);
        y = x;  if ~isnan(orders(m)), y = removePoly(x, t, orders(m)); end
        R = classifyShapeNoise_v004(y, c.fs, band);
        cls(r, m) = R.cls;
        if R.id ~= "", ids(end + 1, 1) = manip(m) + ":" + R.id; end %#ok<AGROW>
        if ~isempty(R.info), LN(r, m) = R.info.Lfrac; end
    end
end
row = table();
for m = 1:nM
    row = [row; table(c.fs, c.N, gen, c.alpha, manip(m), nRep, nnz(cls(:, m) == "ok"), nnz(cls(:, m) == "noSegment"), nnz(cls(:, m) == "segmentAlpha"), ...
        nnz(cls(:, m) == "noAccepted"), nnz(cls(:, m) == "other"), nnz(cls(:, m) == "noReference"), median(LN(:, m), 'omitnan'), ...
        median(trendR2, 'omitnan'), median(jump, 'omitnan'), median(cross, 'omitnan'), ...
        'VariableNames', {'fs', 'N', 'generator', 'alpha', 'manipulation', 'nRep', 'ok', 'noSegment', 'segmentAlpha', 'noAccepted', 'other', 'noReference', ...
        'LfracMedian', 'trendR2Median', 'endJumpRatioMedian', 'meanCrossingsMedian'})]; %#ok<AGROW>
end
[u, ~, g] = unique(ids);
if isempty(u)
    E = table(zeros(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), 'VariableNames', {'fs', 'N', 'generator', 'alpha', 'identifier', 'count'});
else
    k = numel(u);
    E = table(repmat(c.fs, k, 1), repmat(c.N, k, 1), repmat(gen, k, 1), repmat(c.alpha, k, 1), u(:), accumarray(g(:), 1), 'VariableNames', {'fs', 'N', 'generator', 'alpha', 'identifier', 'count'});
end
end

function y = removePoly(x, t, order)
A = t .^ (0:order);
y = x - A * (A \ x);
end
