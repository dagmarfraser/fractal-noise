% runCheckShapeNoiseEndMatched_v003.m  Runner for checkShapeNoiseEndMatched_v003 (high-alpha validation on break-free exact references).
% Run from a shell:  matlab -batch "SMOKE=false; run('<this file>')"     SMOKE=true: 2 alphas, 2 configurations, 4 replicates, separate output folder.
% Output folder beside this file: out_highAlpha_v003 (or out_highAlpha_v003_smoke): cells.csv, errors.csv, results.mat, log.txt.
% Stops rather than overwrites existing output. Toolbox, exact references and the band-effective target come from this clone only.
if ~exist('SMOKE', 'var'), SMOKE = false; end
here = fileparts(mfilename('fullpath'));
outDir = fullfile(here, "out_highAlpha_v003" + string(repmat('_smoke', 1, SMOKE)));
if isfolder(outDir), error('runCheck:exists', 'Output folder exists; use a new version: %s', outDir); end
mkdir(outDir);
diary(fullfile(outDir, 'log.txt'));  cleanupDiary = onCleanup(@() diary('off'));
fprintf('runCheckShapeNoiseEndMatched_v003 %s  MATLAB %s  SMOKE=%d\n', string(datetime('now')), version, SMOKE);
paths = {fullfile(here, '..', 'toolbox'), here, fullfile(here, 'exactReference'), fullfile(here, 'stage1')};
addpath(paths{:});
fprintf('shapeNoise: %s\nestimateWhittleAlpha: %s\nexactRefSeries_v001: %s\n', which('shapeNoise'), which('estimateWhittleAlpha'), which('exactRefSeries_v001'));
pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', 6); end
wait(parfevalOnAll(pool, @addpath, 0, paths{:}));
fprintf('workers: %d\n', pool.NumWorkers);
tStart = tic;
if SMOKE
    [T, E] = checkShapeNoiseEndMatched_v003(Alphas=[3.25 6.75], Config=[120 1000; 300 3000], NRep=4);
else
    [T, E] = checkShapeNoiseEndMatched_v003();
end
elapsedMin = toc(tStart) / 60;
writetable(T, fullfile(outDir, 'cells.csv'));
writetable(E, fullfile(outDir, 'errors.csv'));
save(fullfile(outDir, 'results.mat'), 'T', 'E', 'elapsedMin', '-v7.3');
fprintf('cells: %d; series: %d; error rows: %d; elapsed %.1f min\n', height(T), sum(T.nRep), height(E), elapsedMin);
fprintf('DONE\n');
