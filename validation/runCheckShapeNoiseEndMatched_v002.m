% runCheckShapeNoiseEndMatched_v002.m  Runner for checkShapeNoiseEndMatched_v002 (low-alpha extension of the shapeNoise validation).
% Run from a shell:  matlab -batch "SMOKE=false; run('<this file>')"     SMOKE=true: 2 alphas, 2 configurations, 4 replicates, separate output folder.
% Output folder beside this file: out_lowAlpha_v002 (or out_lowAlpha_v002_smoke): cells.csv, summary.csv, errors.csv, results.mat, log.txt.
% Stops rather than overwrites existing output. Toolbox, exact references and the band-effective target come from this clone only.
if ~exist('SMOKE', 'var'), SMOKE = false; end
here = fileparts(mfilename('fullpath'));
outDir = fullfile(here, "out_lowAlpha_v002" + string(repmat('_smoke', 1, SMOKE)));
if isfolder(outDir), error('runCheck:exists', 'Output folder exists; use a new version: %s', outDir); end
mkdir(outDir);
diary(fullfile(outDir, 'log.txt'));  cleanupDiary = onCleanup(@() diary('off'));
fprintf('runCheckShapeNoiseEndMatched_v002 %s  MATLAB %s  SMOKE=%d\n', string(datetime('now')), version, SMOKE);
paths = {fullfile(here, '..', 'toolbox'), here, fullfile(here, 'exactReference'), fullfile(here, 'stage1')};
addpath(paths{:});
fprintf('shapeNoise: %s\nestimateWhittleAlpha: %s\nexactRefSeries_v001: %s\n', which('shapeNoise'), which('estimateWhittleAlpha'), which('exactRefSeries_v001'));
pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', 6); end
wait(parfevalOnAll(pool, @addpath, 0, paths{:}));
fprintf('workers: %d\n', pool.NumWorkers);
tStart = tic;
if SMOKE
    [T, E] = checkShapeNoiseEndMatched_v002(Alphas=[-0.5 1.05], Config=[120 1000; 300 3000], NRep=4);
else
    [T, E] = checkShapeNoiseEndMatched_v002();
end
elapsedMin = toc(tStart) / 60;
writetable(T, fullfile(outDir, 'cells.csv'));
writetable(E, fullfile(outDir, 'errors.csv'));

% Summary per configuration, generator and alpha range (worst cell in each range)
ranges = [-1.01 -0.01; 0 2.01; 2.02 3.06];  names = ["alpha -1 to -0.25", "alpha 0 to 2", "alpha 2.25 to 3.05"];
S = table();
[G, fsG, NG, genG] = findgroups(T.fs, T.N, T.generator);
for g = 1:max(G)
    for k = 1:size(ranges, 1)
        q = T(G == g & T.alpha >= ranges(k, 1) & T.alpha <= ranges(k, 2), :);
        if isempty(q), continue, end
        S = [S; table(fsG(g), NG(g), genG(g), names(k), height(q), ...
            max(abs(q.biasEffective)), max(abs(q.biasNominal)), max(q.rmseEffective), min(q.shareWithin025Effective), ...
            max(abs(q.surrBiasEffective)), max(abs(q.gapSelfMedian)), min(q.LfracMedian), min(q.LfracMin), min(q.d1Median), max(q.d1Median), min(q.d2Median), max(q.d2Median), ...
            sum(q.errReference), sum(q.errEstimateRef), sum(q.errShapeNoise), sum(q.errEstimateSurrogate), sum(q.nNoAcceptedDraw), sum(q.nEstimateFailedDraws), ...
            'VariableNames', {'fs', 'N', 'generator', 'alphaRange', 'nCells', 'worstAbsBiasEffective', 'worstAbsBiasNominal', 'worstRmseEffective', 'worstShareWithin025', ...
            'worstAbsSurrBiasEffective', 'worstAbsGapSelfMedian', 'minLfracMedian', 'minLfracAny', 'd1Low', 'd1High', 'd2Low', 'd2High', ...
            'errReference', 'errEstimateRef', 'errShapeNoise', 'errEstimateSurrogate', 'nNoAcceptedDraw', 'nEstimateFailedDraws'})]; %#ok<AGROW>
    end
end
writetable(S, fullfile(outDir, 'summary.csv'));
save(fullfile(outDir, 'results.mat'), 'T', 'E', 'S', 'elapsedMin', '-v7.3');
fprintf('cells: %d; series: %d; errors listed: %d rows; elapsed %.1f min\n', height(T), sum(T.nRep), height(E), elapsedMin);
disp(S(S.generator == "exact", {'fs', 'N', 'alphaRange', 'worstAbsBiasEffective', 'worstAbsBiasNominal', 'worstRmseEffective', 'worstShareWithin025', 'errEstimateRef', 'errShapeNoise', 'nNoAcceptedDraw'}));
fprintf('DONE\n');
