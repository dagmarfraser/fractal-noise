function T = runMutationChecks(options)
%RUNMUTATIONCHECKS Do the tests catch deliberate bugs? Apply one small mutation at a time and check that the tests fail.
%
%   T = runMutationChecks()
%   T = runMutationChecks(Only="harmonicTemplate")
%
%   A passing suite shows the code agrees with the tests; it does not show the tests can fail. This script copies toolbox/ and
%   tests/ to a temporary folder, makes ONE textual change to one function (a mutant: drop the Hann weights, ignore an option,
%   loosen a bound, ...), runs the tests that cover that function and records whether any of them failed (the mutant is KILLED).
%   A mutant that survives means the tests cannot see that bug. The repository is never modified. First a baseline run with no
%   mutation must pass, so a kill cannot come from the temporary copy itself.
%
%   Each mutation's target text must occur exactly once in its file; otherwise the script errors rather than skipping it.
%
% OUTPUT:
%   T - table: id, file, change, killed (logical), nFailed, nTests.
%
% See also: runtests

arguments
    options.Only (1, 1) string = ""
end
here = fileparts(mfilename("fullpath"));  root = fileparts(here);
M = mutants();
if options.Only ~= "", M = M(contains([M.file], options.Only)); end
tests = ["testHarmonicTemplate" "testShapeNoiseGenericity" "testWrappers" "testHarmonicContent" "testShapeNoise"];

oldPath = path;  cleanupPath = onCleanup(@() path(oldPath));  %#ok<NASGU>
oldWarn = warning("on", "all");  cleanupWarn = onCleanup(@() warning(oldWarn));  %#ok<NASGU>   % warnings must be ON: several tests verify them

[ok0, n0] = runIn(root, "", "", "", tests);   % baseline: no mutation
if ~ok0, error("runMutationChecks:baseline", "The unmutated copy does not pass its own tests (%d tests); fix that first.", n0); end

id = strings(0, 1);  file = id;  change = id;  killed = false(0, 1);  nFailed = zeros(0, 1);  nTests = zeros(0, 1);
for k = 1:numel(M)
    m = M(k);
    [~, nf, nt] = runIn(root, m.file, m.from, m.to, m.tests);
    id(end + 1, 1) = m.id;  file(end + 1, 1) = m.file;  change(end + 1, 1) = m.change;  killed(end + 1, 1) = nf > 0;  nFailed(end + 1, 1) = nf;  nTests(end + 1, 1) = nt; %#ok<AGROW>
    fprintf("%-4s %-24s %-8s %s\n", m.id, m.file, ternaryStr(nf > 0, "killed", "SURVIVED"), m.change);
end
T = table(id, file, change, killed, nFailed, nTests);
fprintf("\n%d of %d mutants killed.\n", nnz(killed), numel(killed));
end

function s = ternaryStr(c, a, b)
if c, s = a; else, s = b; end
end

function [allPass, nFailed, nTests] = runIn(root, file, from, to, testNames)
tmp = tempname;  mkdir(tmp);  cleanup = onCleanup(@() rmdir(tmp, "s"));  %#ok<NASGU>
copyfile(fullfile(root, "toolbox"), fullfile(tmp, "toolbox"));
copyfile(fullfile(root, "tests"), fullfile(tmp, "tests"));
if file ~= ""
    f = fullfile(tmp, "toolbox", file + ".m");  txt = string(fileread(f));
    nOcc = count(txt, from);
    if nOcc ~= 1, error("runMutationChecks:target", "Mutation target occurs %d times (needs exactly 1) in %s: %s", nOcc, file, from); end
    txt = replace(txt, from, to);
    fid = fopen(f, "w");  fwrite(fid, char(txt));  fclose(fid);
end
p0 = path;  addpath(fullfile(tmp, "toolbox"));  restore = onCleanup(@() path(p0));  %#ok<NASGU>
rehash;  clear functions;
nFailed = 0;  nTests = 0;
for tn = testNames
    f = fullfile(tmp, "tests", tn + ".m");
    if ~isfile(f), continue, end
    res = [];
    evalc("res = runtests(f, ""Verbosity"", 0);");
    nTests = nTests + numel(res);  nFailed = nFailed + nnz([res.Failed]);
end
allPass = (nFailed == 0);
end

function M = mutants()
mk = @(id, file, tests, change, from, to) struct("id", id, "file", file, "tests", tests, "change", change, "from", from, "to", to);
HT = ["testHarmonicTemplate" "testWrappers"];
SN = ["testShapeNoise" "testShapeNoiseGenericity" "testWrappers"];
HC = "testHarmonicContent";
M = [
 mk("T1", "harmonicTemplate", HT, "drop the Hann weights", "han = 0.5 * (1 - cos(2 * pi * (0:win - 1)' / (win - 1)));", "han = ones(win, 1);")
 mk("T2", "harmonicTemplate", HT, "remove one harmonic from the design", "for h = 1:K, D = [D,", "for h = 1:max(K - 1, 0), D = [D,")
 mk("T3", "harmonicTemplate", HT, "no end-anchored final window", "if options.EndAnchored && starts(end) + win - 1 < N", "if false && starts(end) + win - 1 < N")
 mk("T4", "harmonicTemplate", HT, "crop an eighth of a window, not a quarter", "cl = max(round(win / 4), 1);", "cl = max(round(win / 8), 1);")
 mk("T5", "harmonicTemplate", HT, "no tooShort error", "if win > N || win < 2 * (2 * K + 2)", "if false")
 mk("T6", "harmonicTemplate", HT, "no overlap-add normalisation", "res(covered, :) = res(covered, :) ./ ww(covered);", "res(covered, :) = res(covered, :);")
 mk("T7", "harmonicTemplate", HT, "ignore the Trend option", "D = ones(win, 1);  if options.Trend, D = [D, tw]; end", "D = ones(win, 1);")
 mk("T8", "harmonicTemplate", HT, "ignore the Nyquist check", "if K * f0 >= fs / 2, error(", "if false, error(")
 mk("S1", "shapeNoise", SN, "ignore CheckUncovered (always check)", "if options.CheckUncovered, checkZeroTail(resid); end", "checkZeroTail(resid);")
 mk("S2", "shapeNoise", SN, "default band no longer warns", "warning(""shapeNoise:defaultAlphaBand"",", "warning(""shapeNoise:mutated"",")
 mk("S3", "shapeNoise", SN, "allow fLow = 0 in AlphaBand", "ab(1) > 0 && ab(2) > ab(1)", "ab(2) > ab(1)")
 mk("S4", "shapeNoise", SN, "allow fHigh = fs/2", "if band(2) >= fs/2 || band(2) <= band(1)", "if band(2) > fs/2 || band(2) <= band(1)")
 mk("S5", "shapeNoise", SN, "zero-tail threshold 3 to 4", "if lead >= 3 || trail >= 3", "if lead >= 4 || trail >= 4")
 mk("S6", "shapeNoise", SN, "kinematic band upper edge 20 to 25", "kin = @() [1, min(20, fs/2 - 1)];", "kin = @() [1, min(25, fs/2 - 1)];")
 mk("C1", "harmonicContent", HC, "count stops one harmonic late", "K = max(1, below - 1);", "K = max(1, below);")
 mk("C2", "harmonicContent", HC, "no between-harmonic noise reference", "ampNull(k) = sqrt(Am(2 * k - 1) * Am(2 * k + 1));", "ampNull(k) = 1;")
 mk("C4", "harmonicContent", HC, "noHarmonic flag always false", """noHarmonic"", ~(ratio(1) >= options.Cut)", """noHarmonic"", false")
 mk("C5", "harmonicContent", HC, "channels not averaged (first only)", "Am = mean(A, 1);", "Am = A(1, :);")
 mk("W1", "cetNoiseWrapper", "testWrappers", "no LastYear guard", "if lastYear > lastComplete, error(", "if false, error(")
 mk("W2", "cetNoiseWrapper", "testWrappers", "missing months not detected", "if any(R(:, 2:13) <= -99, ""all""), error(", "if false, error(")
 mk("W3", "shapeNoiseWrapper", "testWrappers", "template step never runs", "useTemplate = ~isempty(options.F0) || all(~isnan(options.F0Band));", "useTemplate = false;")
 ];
end
