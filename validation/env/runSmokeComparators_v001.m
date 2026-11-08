function T = runSmokeComparators_v001(opts)
%RUNSMOKECOMPARATORS_V001 One real call per Python/R comparator on shared exact fGn.
%   Input: exact fGn, H = 0.75 (alpha = 0.5), N = 4096 (exactRefSeries_v001).
%   Estimators should land near 0.5; the generator rows (colorednoise /
%   neurodsp / fbm) are scored on their own target alpha, given in the name.
%   Pass = ran and returned a finite value within +/-0.3 of its target
%   (one realisation, so this is a smoke test, not calibration). Errors if any fail.
arguments
    opts.PythonExe (1,1) string = "~/venvs/fnval/bin/python"
    opts.RscriptExe (1,1) string = "Rscript"
end
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "exactReference"));
py = regexprep(opts.PythonExe, "^~", getenv("HOME"));
x  = exactRefSeries_v001(4096, 0.5, 1, Stream=RandStream('threefry', Seed=11));
f  = fullfile(tempdir, "fn_smoke_fgn.csv");  writematrix(x, f);  c = onCleanup(@() delete(f));
[~, oP] = system(sprintf('"%s" "%s" "%s" 100 2>&1', py, fullfile(here, "smoke_comparators_v001.py"), f));
[~, oR] = system(sprintf('"%s" "%s" "%s" 2>&1', opts.RscriptExe, fullfile(here, "smoke_comparators_v001.R"), f));
L = splitlines(string(oP) + newline + string(oR));  L = L(startsWith(L, "RES|"));
p = split(extractAfter(L, "RES|"), "|", 2);
T = table(p(:,1), str2double(p(:,2)), p(:,3), VariableNames=["tool","value","status"]);
target = 0.5 * ones(height(T), 1);
target(contains(T.tool, "_a2")) = 2;
T.target = target;
T.pass = T.status == "ok" & isfinite(T.value) & abs(T.value - T.target) <= 0.3;
disp(T)
if height(T) ~= 15; error('smoke:count', 'Expected 15 RES lines, got %d. Raw output:\n%s\n%s', height(T), oP, oR); end
if ~all(T.pass); error('smoke:fail', '%d comparator(s) FAILED: %s', sum(~T.pass), strjoin(T.tool(~T.pass), ", ")); end
fprintf('All 15 comparators ran and landed near target.\n');
end
