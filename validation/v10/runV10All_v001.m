function runV10All_v001()
%RUNV10ALL_V001 Detached launcher for V10 (lean): exact-reference null, then the three real datasets.
%   matlab -batch "cd('<fractalnoise>/validation/v10'); runV10All_v001"
%   Null -> validation/v10/runs/null_v001/; residuals -> PowerLawSimulationPreReg/results/residualMFDFA/.
%   Hickman PLAC and HALO are run separately (never pooled).
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..', 'toolbox'));
pre = fullfile(here, '..', '..', '..', '..', 'PowerLawSimulationPreReg', 'src');
if ~isfile(fullfile(pre, 'residualMFDFA_v001.m')), error('v10:path', 'residualMFDFA_v001 not found in %s.', pre); end
addpath(pre);
fprintf('[%s] null start\n', char(datetime('now')));
runV10Null_v001(fullfile(here, 'runs', 'null_v001'), 200);
for ds = ["FraserN48" "Hickman PLAC" "Hickman HALO"]
    fprintf('[%s] %s start\n', char(datetime('now')), ds);
    residualMFDFA_v001(ds);
end
fprintf('[%s] ALL DONE\n', char(datetime('now')));
end
