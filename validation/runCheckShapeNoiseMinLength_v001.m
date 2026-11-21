% runCheckShapeNoiseMinLength_v001.m  Rerun the v002 (low alpha) and v003 (high alpha) shapeNoise grids with MinLength = 120 (default 32).
% Run from a shell:  matlab -batch "SMOKE=false; run('<this file>')"     SMOKE=true: the runners' smoke grids, separate output folders.
% Uses checkShapeNoiseEndMatched_v005 and _v006, which are v002 and v003 with one added option (ShapeArgs); v002 and v003 are not edited.
% Before anything else it checks, on the smoke grids, that (1) v005/v006 with ShapeArgs = {} give exactly what v002/v003 give, and
% (2) both match the saved 18 Nov smoke outputs (out_lowAlpha_v002_smoke, out_highAlpha_v003_smoke) to 1e-9. Either failing stops the run.
% Output folders beside this file: out_lowAlpha_minL120_v001 and out_highAlpha_minL120_v001 (or _smoke): cells.csv, errors.csv, results.mat, log.txt.
% Stops rather than overwrites existing output.
if ~exist('SMOKE', 'var'), SMOKE = false; end
here = fileparts(mfilename('fullpath'));
suffix = string(repmat('_smoke', 1, SMOKE));
outLow = fullfile(here, "out_lowAlpha_minL120_v001" + suffix);
outHigh = fullfile(here, "out_highAlpha_minL120_v001" + suffix);
if isfolder(outLow) || isfolder(outHigh), error('runCheck:exists', 'Output folder exists; use a new version: %s', outLow); end
mkdir(outLow);  mkdir(outHigh);
diary(fullfile(outLow, 'log.txt'));  cleanupDiary = onCleanup(@() diary('off'));
fprintf('runCheckShapeNoiseMinLength_v001 %s  MATLAB %s  SMOKE=%d\n', string(datetime('now')), version, SMOKE);
paths = {fullfile(here, '..', 'toolbox'), here, fullfile(here, 'exactReference'), fullfile(here, 'stage1')};
addpath(paths{:});
fprintf('shapeNoise: %s\n', which('shapeNoise'));
pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', 6); end
wait(parfevalOnAll(pool, @addpath, 0, paths{:}));
fprintf('workers: %d\n', pool.NumWorkers);

% --- equivalence and reproduction checks on the smoke grids (always run) ---
smokeLow = {'Alphas', [-0.5 1.05], 'Config', [120 1000; 300 3000], 'NRep', 4};
smokeHigh = {'Alphas', [3.25 6.75], 'Config', [120 1000; 300 3000], 'NRep', 4};
[T2, E2] = checkShapeNoiseEndMatched_v002(smokeLow{:});
[T5, E5] = checkShapeNoiseEndMatched_v005(smokeLow{:});
if ~isequaln(T2, T5) || ~isequaln(E2, E5), error('runCheck:notEquivalent', 'v005 with ShapeArgs = {} differs from v002.'); end
[T3, E3] = checkShapeNoiseEndMatched_v003(smokeHigh{:});
[T6, E6] = checkShapeNoiseEndMatched_v006(smokeHigh{:});
if ~isequaln(T3, T6) || ~isequaln(E3, E6), error('runCheck:notEquivalent', 'v006 with ShapeArgs = {} differs from v003.'); end
fprintf('EQUIVALENT: v005 = v002 and v006 = v003 on the smoke grids (ShapeArgs = {}).\n');
reproduces(T5, fullfile(here, 'out_lowAlpha_v002_smoke', 'cells.csv'), 'low');
reproduces(T6, fullfile(here, 'out_highAlpha_v003_smoke', 'cells.csv'), 'high');

% --- the reruns with MinLength = 120 ---
args = {'ShapeArgs', {"MinLength", 120}};
tStart = tic;
if SMOKE
    [TL, EL] = checkShapeNoiseEndMatched_v005(smokeLow{:}, args{:});
else
    [TL, EL] = checkShapeNoiseEndMatched_v005(args{:});
end
minLow = toc(tStart) / 60;
writetable(TL, fullfile(outLow, 'cells.csv'));  writetable(EL, fullfile(outLow, 'errors.csv'));
save(fullfile(outLow, 'results.mat'), 'TL', 'EL', 'minLow', '-v7.3');
fprintf('low alpha, MinLength 120: cells %d; series %d; error rows %d; elapsed %.1f min\n', height(TL), sum(TL.nRep), height(EL), minLow);
tStart = tic;
if SMOKE
    [TH, EH] = checkShapeNoiseEndMatched_v006(smokeHigh{:}, args{:});
else
    [TH, EH] = checkShapeNoiseEndMatched_v006(args{:});
end
minHigh = toc(tStart) / 60;
writetable(TH, fullfile(outHigh, 'cells.csv'));  writetable(EH, fullfile(outHigh, 'errors.csv'));
save(fullfile(outHigh, 'results.mat'), 'TH', 'EH', 'minHigh', '-v7.3');
fprintf('high alpha, MinLength 120: cells %d; series %d; error rows %d; elapsed %.1f min\n', height(TH), sum(TH.nRep), height(EH), minHigh);
fprintf('DONE\n');

function reproduces(T, savedCsv, label)
% The rerun must reproduce the saved 18 Nov smoke output (numeric columns to 1e-9, text columns exactly).
S = readtable(savedCsv, 'TextType', 'string');
if height(S) ~= height(T) || ~isequal(sort(S.Properties.VariableNames), sort(T.Properties.VariableNames))
    error('runCheck:notReproduced', '%s: saved smoke output has a different shape.', label);
end
for v = string(T.Properties.VariableNames)
    a = T.(v);  b = S.(v);
    if isnumeric(a) || islogical(a)
        a = double(a);  b = double(b);
        bad = ~((isnan(a) & isnan(b)) | abs(a - b) <= 1e-9 * max(1, abs(b)));
    else
        bad = string(a) ~= string(b);
    end
    if any(bad), error('runCheck:notReproduced', '%s: column %s differs from %s.', label, v, savedCsv); end
end
fprintf('REPRODUCED: %s smoke grid matches %s.\n', label, savedCsv);
end
