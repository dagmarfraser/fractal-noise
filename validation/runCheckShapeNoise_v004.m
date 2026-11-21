% runCheckShapeNoise_v004.m  Runner for checkShapeNoisePhiMechanism_v004 and checkShapeNoiseLayers_v004.
% Run from a shell:  matlab -batch "SMOKE=false; run('<this file>')"     SMOKE=true: a small grid and a separate output folder.
% Output folder beside this file: out_v004 (or out_v004_smoke): phi_cells.csv, phi_errors.csv, layers_cells.csv, layers_errors.csv, results.mat, log.txt.
% Stops rather than overwrites existing output. Toolbox, exact references and the band-effective target come from this clone only.
if ~exist('SMOKE', 'var'), SMOKE = false; end
here = fileparts(mfilename('fullpath'));
outDir = fullfile(here, "out_v004" + string(repmat('_smoke', 1, SMOKE)));
if isfolder(outDir), error('runCheck:exists', 'Output folder exists; use a new version: %s', outDir); end
mkdir(outDir);
diary(fullfile(outDir, 'log.txt'));  cleanupDiary = onCleanup(@() diary('off'));
fprintf('runCheckShapeNoise_v004 %s  MATLAB %s  SMOKE=%d\n', string(datetime('now')), version, SMOKE);
paths = {fullfile(here, '..', 'toolbox'), here, fullfile(here, 'exactReference'), fullfile(here, 'stage1')};
addpath(paths{:});
fprintf('shapeNoise: %s\nestimateWhittleAlpha: %s\nnoiseXu: %s\n', which('shapeNoise'), which('estimateWhittleAlpha'), which('noiseXu'));

% The extended reference must reproduce exactRefSeries_v001 wherever both exist
for a = [-0.5 0.5 2.5 4.5 6.5]
    rng(1, 'twister');  x1 = exactRefSeries_v001(500, a, 1);
    rng(1, 'twister');  x2 = exactRefExtended_v004(500, a);
    if ~isequal(x1, x2), error('runCheck:referenceMismatch', 'exactRefExtended_v004 differs from exactRefSeries_v001 at alpha %g.', a); end
end
fprintf('exactRefExtended_v004 reproduces exactRefSeries_v001 at alpha -0.5, 0.5, 2.5, 4.5, 6.5.\n');

pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', 6); end
wait(parfevalOnAll(pool, @addpath, 0, paths{:}));
fprintf('workers: %d\n', pool.NumWorkers);
tStart = tic;
if SMOKE
    [P, PE] = checkShapeNoisePhiMechanism_v004(Alphas=[2.95 3], Config=[120 1000], NRep=4);
else
    [P, PE] = checkShapeNoisePhiMechanism_v004();
end
writetable(P, fullfile(outDir, 'phi_cells.csv'));  writetable(PE, fullfile(outDir, 'phi_errors.csv'));
phiMin = toc(tStart) / 60;
fprintf('Phi mechanism: %d rows, %.1f min\n', height(P), phiMin);
tStart = tic;
if SMOKE
    [T, E] = checkShapeNoiseLayers_v004(Alphas=[-4 -3 -2.5 2.95 3 8 10], Config=[120 1000; 300 3000], NRep=4);
else
    [T, E] = checkShapeNoiseLayers_v004();
end
writetable(T, fullfile(outDir, 'layers_cells.csv'));  writetable(E, fullfile(outDir, 'layers_errors.csv'));
layersMin = toc(tStart) / 60;
save(fullfile(outDir, 'results.mat'), 'P', 'PE', 'T', 'E', 'phiMin', 'layersMin', '-v7.3');
fprintf('Layers: %d cells, %d series, %.1f min\n', height(T), sum(T.nRep), layersMin);
fprintf('DONE\n');
