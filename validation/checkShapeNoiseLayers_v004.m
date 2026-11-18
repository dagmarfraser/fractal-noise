function [T, E] = checkShapeNoiseLayers_v004(opts)
%CHECKSHAPENOISELAYERS_V004 Three-layer score (generator, estimator, shapeNoise) over alpha -4 to 10.
%   [T, E] = checkShapeNoiseLayers_v004() runs the grid; T has one row per cell, E lists every error by cell and identifier.
%
%   Grid: alpha -4:0.5:10 plus +/-0.05 brackets at each odd integer from -3 to 9; fs 120, 133, 240, 300 Hz crossed with N 1000, 3000;
%   generators "exact" (exactRefExtended_v004; no row at odd integers, where no exact series exists), "noiseXu_Phi0.99" and "noiseXu_Phi1";
%   100 replicates per cell; band 1-20 Hz.
%
%   LAYER A, generator (no estimateWhittleAlpha): the series are differenced k = max(0, floor((alpha + 1) / 2)) times with k known from the
%   target, Hann-tapered, and their periodograms recoloured by |1 - exp(-i lambda)|^(-2k) and averaged over replicates. The in-band OLS
%   slope of the log mean periodogram is compared with the OLS slope of the analytic spectrum at the same bins: aA - aTarget, its Monte
%   Carlo SE from per-replicate slopes, and the RMS shape error in dB after removing the level.
%   Analytic spectra: exact, info.S; noiseXu, |1 - Phi exp(-i lambda)|^(-alpha) (Xu 2019 eq. 8).
%   LAYER B, estimator: estimateWhittleAlpha on each series against the band-effective alpha (bandEffectiveAlpha_v001). The exact-generator
%   row of a cell is layer B; noiseXu rows carry the same columns for information only.
%   LAYER C, shapeNoise: the default call on each series. Outcome counts by class, and the gap between the surrogate alpha and the alpha of
%   the whole input series. The pass thresholds and the attribution are applied afterwards, in reportLayers_v004.py.
%
%   Seeds: rng(1e9*Master + 1e8*config + 1e7*generator + 1e4*alphaIndex + replicate, 'twister'), Master = 2 (v002 used 1, v003 used 3).
%   Errors are counted per cell by stage and identifier, never hidden; a failed replicate is NaN.
arguments
    opts.Alphas (1, :) double = sort([-4:0.5:10, (-3:2:9) - 0.05, (-3:2:9) + 0.05])
    opts.Config (:, 2) double = [120 1000; 120 3000; 133 1000; 133 3000; 240 1000; 240 3000; 300 1000; 300 3000]
    opts.NRep (1, 1) double {mustBeInteger, mustBePositive} = 100
    opts.Master (1, 1) double {mustBeInteger, mustBePositive} = 2
    opts.Band (1, 2) double = [1 20]
end
gens = ["exact", "noiseXu_Phi0.99", "noiseXu_Phi1"];
jobs = struct('fs', {}, 'N', {}, 'gen', {}, 'alpha', {}, 'ci', {}, 'gi', {}, 'ai', {});
for ci = 1:size(opts.Config, 1)
    for gi = 1:numel(gens)
        for ai = 1:numel(opts.Alphas)
            a = opts.Alphas(ai);
            if gens(gi) == "exact" && abs(mod(a + 1, 2)) < 1e-9, continue, end   % odd integer: not constructible
            jobs(end + 1) = struct('fs', opts.Config(ci, 1), 'N', opts.Config(ci, 2), 'gen', gens(gi), 'alpha', a, 'ci', ci, 'gi', gi, 'ai', ai); %#ok<AGROW>
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
S = targetSpectrum(c);
aEff = bandEffectiveAlpha_v001(S, c.fs, band, c.N);
X = NaN(c.N, nRep);  aIn = NaN(nRep, 1);  aSur = aIn;  gapIn = aIn;  gapSelf = aIn;  LN = aIn;  cls = strings(nRep, 1);
ids = strings(0, 1);  estFail = 0;
for r = 1:nRep
    rng(1e9 * master + 1e8 * c.ci + 1e7 * c.gi + 1e4 * c.ai + r, 'twister');
    try
        x = generate(c);
        if any(~isfinite(x)) || std(x) == 0, error('checkLayers:badSeries', 'Series is non-finite or constant.'); end
    catch err, ids(end + 1, 1) = "reference:" + err.identifier;  cls(r) = "noReference";  continue %#ok<AGROW>
    end
    X(:, r) = x;
    try, aIn(r) = estimateWhittleAlpha(x, c.fs, band(1), band(2));
    catch err, ids(end + 1, 1) = "estimateInput:" + err.identifier; %#ok<AGROW>
    end
    R = classifyShapeNoise_v004(x, c.fs, band);
    cls(r) = R.cls;
    if R.id ~= "", ids(end + 1, 1) = "shapeNoise:" + R.id; end %#ok<AGROW>
    if isempty(R.info), continue, end
    estFail = estFail + sum(R.info.nEstFail);  LN(r) = R.info.Lfrac;
    if R.cls ~= "ok", continue, end
    try
        aSur(r) = estimateWhittleAlpha(R.s, c.fs, band(1), band(2));
        gapIn(r) = aSur(r) - aIn(r);  gapSelf(r) = aSur(r) - R.info.alphaTarget;
    catch err, ids(end + 1, 1) = "estimateSurrogate:" + err.identifier; %#ok<AGROW>
    end
end
good = all(isfinite(X), 1);
A = struct('aA', NaN, 'aTarget', NaN, 'gA', NaN, 'mcseA', NaN, 'rmsDb', NaN, 'k', NaN);
try, A = layerA(X(:, good), c, band, S);
catch err, ids(end + 1, 1) = "layerA:" + err.identifier; %#ok<AGROW>
end
ok = ~isnan(aIn);  n = nnz(ok);  nGap = nnz(~isnan(gapIn));
row = table(c.fs, c.N, c.gen, c.alpha, aEff, nRep, nnz(good), ...
    A.k, A.aA, A.aTarget, A.gA, A.mcseA, A.rmsDb, ...
    n, mean(aIn, 'omitnan') - c.alpha, mean(aIn, 'omitnan') - aEff, std(aIn, 'omitnan'), std(aIn, 'omitnan') / sqrt(max(n, 1)), sqrt(mean((aIn - aEff) .^ 2, 'omitnan')), ...
    nnz(cls == "ok"), nnz(cls == "noSegment"), nnz(cls == "segmentAlpha"), nnz(cls == "noAccepted"), nnz(cls == "other"), estFail, median(LN, 'omitnan'), ...
    nGap, median(gapIn, 'omitnan'), 1.253 * std(gapIn, 'omitnan') / sqrt(max(nGap, 1)), median(gapSelf, 'omitnan'), mean(aSur, 'omitnan') - aEff, ...
    'VariableNames', {'fs', 'N', 'generator', 'alpha', 'alphaEffective', 'nRep', 'nGenerated', ...
    'kPrewhiten', 'aLayerA', 'aTargetA', 'gapA', 'mcseA', 'shapeRmsDb', ...
    'nEstimated', 'biasNominal', 'biasEffective', 'sdAlpha', 'mcseBias', 'rmseEffective', ...
    'okShape', 'noSegment', 'segmentAlpha', 'noAccepted', 'otherShape', 'nEstimateFailedDraws', 'LfracMedian', ...
    'nGap', 'gapInputMedian', 'mcseGapInput', 'gapSelfMedian', 'surrBiasEffective'});
[u, ~, g] = unique(ids);
if isempty(u)
    E = table(zeros(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), strings(0, 1), zeros(0, 1), 'VariableNames', {'fs', 'N', 'generator', 'alpha', 'identifier', 'count'});
else
    m = numel(u);
    E = table(repmat(c.fs, m, 1), repmat(c.N, m, 1), repmat(c.gen, m, 1), repmat(c.alpha, m, 1), u(:), accumarray(g(:), 1), 'VariableNames', {'fs', 'N', 'generator', 'alpha', 'identifier', 'count'});
end
end

function x = generate(c)
switch c.gen
    case "exact",           x = exactRefExtended_v004(c.N, c.alpha);
    case "noiseXu_Phi0.99", x = noiseXu(c.N, c.alpha, 1, c.fs);
    case "noiseXu_Phi1",    x = noiseXu(c.N, c.alpha, 1, c.fs, Phi=1);
end
x = x(:);
end

function S = targetSpectrum(c)
switch c.gen
    case "exact",           [~, info] = exactRefExtended_v004(8, c.alpha);  S = info.S;
    case "noiseXu_Phi0.99", S = ggmSpectrum(0.99, c.alpha);
    case "noiseXu_Phi1",    S = ggmSpectrum(1, c.alpha);
end
end

function S = ggmSpectrum(phi, alpha)
S = @(lam) (1 - 2 * phi * cos(lam(:)) + phi ^ 2) .^ (-alpha / 2);
end

function A = layerA(X, c, band, S)
% Generator check without the estimator: prewhiten with the known order, Hann periodogram, recolour, average, OLS slope in the band.
if size(X, 2) < 2, error('checkLayers:tooFewSeries', 'Fewer than two series generated.'); end
k = max(0, floor((c.alpha + 1) / 2));
Xd = X;  for j = 1:k, Xd = diff(Xd, 1, 1); end
M = size(Xd, 1);  w = 0.5 * (1 - cos(2 * pi * (0:M - 1)' / M));
I = abs(fft((Xd - mean(Xd, 1)) .* w, [], 1)) .^ 2 / sum(w .^ 2);
f = (0:M - 1)' * c.fs / M;  m = f >= band(1) & f <= band(2) & f < c.fs / 2;
lam = 2 * pi * f(m) / c.fs;
Ib = I(m, :) ./ (2 - 2 * cos(lam)) .^ k;
logF = log10(f(m));  D = [ones(size(logF)), logF];
pMean = D \ log10(mean(Ib, 2));
pRep = D \ log10(Ib);
target = log10(S(lam));  pTarget = D \ target;
res = 10 * (log10(mean(Ib, 2)) - target);
A = struct('aA', -pMean(2), 'aTarget', -pTarget(2), 'gA', pTarget(2) - pMean(2), 'mcseA', std(pRep(2, :)) / sqrt(size(X, 2)), 'rmsDb', rms(res - mean(res)), 'k', k);
end
