function runV10Diag_v001()
%RUNV10DIAG_V001 Detached launcher: V10 null v002 (adds q>0 width), then residualMFDFA_v002 on the three datasets.
%   matlab -batch "cd('<fractalnoise>/validation/v10'); runV10Diag_v001"
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', '..', 'toolbox'));
pre = fullfile(here, '..', '..', '..', '..', 'PowerLawSimulationPreReg', 'src');
if ~isfile(fullfile(pre, 'residualMFDFA_v002.m')), error('v10:path', 'residualMFDFA_v002 not found in %s.', pre); end
addpath(pre);
fprintf('[%s] null v002 start\n', char(datetime('now')));
runV10Null_v002(fullfile(here, 'runs', 'null_v002'), 200);
for ds = ["FraserN48" "Hickman PLAC" "Hickman HALO"]
    fprintf('[%s] %s start\n', char(datetime('now')), ds);
    residualMFDFA_v002(ds);
end
fprintf('[%s] ALL DONE\n', char(datetime('now')));
end
