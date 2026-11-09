function runExtraSeeds_v001(seeds)
%RUNEXTRASEEDS_V001 Detached launcher: full stage-1 grid for extra master seeds, one after another.
%   Meant for a separate `matlab -batch` process so an MCP/desktop MATLAB restart cannot kill it
%   (the first attempt at seeds 20261109/20261110, run inside the MCP session, died after the manifest step).
%   Output: runs/stage1_full_seed<seed>_run2/ per seed (progress.log, estimates_*.csv, summary.csv).
%   Usage (Terminal or system()): matlab -batch "cd('<stage1 dir>'); runExtraSeeds_v001([20261109 20261110])"
arguments
    seeds (1,:) double = [20261109 20261110]
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'toolbox'));
for s = seeds
    out = fullfile(here, 'runs', sprintf('stage1_full_seed%d_run2', s));
    fprintf('[%s] seed %d -> %s\n', char(datetime('now')), s, out);
    runStage1_v002(MasterSeed=s, OutDir=string(out));
    fprintf('[%s] seed %d done\n', char(datetime('now')), s);
end
end
