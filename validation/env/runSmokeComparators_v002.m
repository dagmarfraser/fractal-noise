function S = runSmokeComparators_v002(opts)
%RUNSMOKECOMPARATORS_V002 Multi-seed smoke test of all 15 Python/R comparators.
%   K exact-fGn series (alpha = 0.5, N = 4096, exactRefSeries_v001), each with
%   its own seed from a threefry master stream; Python generators and R
%   surrogates are seeded per series from the same seed list. A tool passes
%   if it never errors, its median is within MedTol of target, and at least
%   FracMin of its K values are within Tol. Still a smoke test, not
%   calibration: it checks the tools behave sanely across seeds.
%   Supersedes v001 (single realisation, unseeded Python generators).
arguments
    opts.K (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.MasterSeed (1,1) double = 20261107
    opts.Tol (1,1) double = 0.3
    opts.MedTol (1,1) double = 0.15
    opts.FracMin (1,1) double = 0.9
    opts.PythonExe (1,1) string = "~/venvs/fnval/bin/python"
    opts.RscriptExe (1,1) string = "Rscript"
end
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "exactReference"));
py = regexprep(opts.PythonExe, "^~", getenv("HOME"));
ms = RandStream('threefry', Seed=opts.MasterSeed);
seeds = randi(ms, 2^31 - 1, opts.K, 1);
X = zeros(4096, opts.K);
for k = 1:opts.K; X(:, k) = exactRefSeries_v001(4096, 0.5, 1, Stream=RandStream('threefry', Seed=seeds(k))); end
fX = fullfile(tempdir, "fn_smoke_X.csv");  fS = fullfile(tempdir, "fn_smoke_seeds.csv");
writematrix(X, fX);  writematrix(seeds, fS);  c = onCleanup(@() delete(fX, fS));
[~, oP] = system(sprintf('"%s" "%s" "%s" 100 "%s" 2>&1', py, fullfile(here, "smoke_comparators_v002.py"), fX, fS));
[~, oR] = system(sprintf('"%s" "%s" "%s" "%s" 2>&1', opts.RscriptExe, fullfile(here, "smoke_comparators_v002.R"), fX, fS));
L = splitlines(string(oP) + newline + string(oR));  L = L(startsWith(L, "RES|"));
p = split(extractAfter(L, "RES|"), "|", 2);
R = table(p(:,1), str2double(p(:,2)), str2double(p(:,3)), p(:,4), VariableNames=["tool","k","value","status"]);
if height(R) ~= 15 * opts.K
    error('smoke:count', 'Expected %d RES lines, got %d. Raw output:\n%s\n%s', 15*opts.K, height(R), oP, oR);
end
tools = unique(R.tool, "stable");  n = numel(tools);
S = table(tools, zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1), false(n,1), ...
    VariableNames=["tool","target","median","p05","p95","fracWithinTol","nFail","pass"]);
for i = 1:n
    r = R(R.tool == tools(i), :);  t = 0.5 + 1.5 * contains(tools(i), "_a2");
    v = r.value(r.status == "ok");
    S.target(i) = t;  S.median(i) = median(v);  q = prctile(v, [5 95]);  S.p05(i) = q(1);  S.p95(i) = q(2);
    S.fracWithinTol(i) = mean(abs(r.value - t) <= opts.Tol);   % NaN (failed) counts as outside
    S.nFail(i) = sum(r.status ~= "ok");
    S.pass(i) = S.nFail(i) == 0 && abs(S.median(i) - t) <= opts.MedTol && S.fracWithinTol(i) >= opts.FracMin;
end
fprintf('K = %d seeds (master %d)\n', opts.K, opts.MasterSeed);  disp(S)
if ~all(S.pass)
    bad = R(R.status ~= "ok", :);  if ~isempty(bad); disp(bad(1:min(5,height(bad)), :)); end
    error('smoke:fail', '%d comparator(s) FAILED across seeds: %s', sum(~S.pass), strjoin(S.tool(~S.pass), ", "));
end
fprintf('All 15 comparators pass across %d seeds.\n', opts.K);
end
