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
tests = ["testHarmonicTemplate" "testShapeNoiseGenericity" "testWrappers" "testHarmonicContent" "testShapeNoise" "testSummariseShapeNoiseBatch" ...
         "testNoiseXu" "testEstimateWhittleAlpha" "testEstimateIRASA" "testCheckSpectralKnee" "testConvertScalingExponent"];

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
if isfolder(fullfile(root, "validation", "exactReference"))   % the exact-reference tests need it
    copyfile(fullfile(root, "validation", "exactReference"), fullfile(tmp, "validation", "exactReference"));
end
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
SB = "testSummariseShapeNoiseBatch";
CS = "testConvertScalingExponent";  NX = "testNoiseXu";  WA = "testEstimateWhittleAlpha";  IR = "testEstimateIRASA";  CK = "testCheckSpectralKnee";
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
 mk("S8", "shapeNoise", SN, "default MinLength back to 32", "mustBeGreaterThanOrEqual(options.MinLength, 5)} = 120", "mustBeGreaterThanOrEqual(options.MinLength, 5)} = 32")
 mk("C1", "harmonicContent", HC, "count stops one harmonic late", "K = max(1, below - 1);", "K = max(1, below);")
 mk("C2", "harmonicContent", HC, "no between-harmonic noise reference", "ampNull(k) = sqrt(Am(2 * k - 1) * Am(2 * k + 1));", "ampNull(k) = 1;")
 mk("C4", "harmonicContent", HC, "noHarmonic flag always false", """noHarmonic"", ~(ratio(1) >= options.Cut)", """noHarmonic"", false")
 mk("C5", "harmonicContent", HC, "channels not averaged (first only)", "Am = mean(A, 1);", "Am = A(1, :);")
 mk("W1", "cetNoiseWrapper", "testWrappers", "no LastYear guard", "if lastYear > lastComplete, error(", "if false, error(")
 mk("W2", "cetNoiseWrapper", "testWrappers", "missing months not detected", "if any(R(:, 2:13) <= -99, ""all""), error(", "if false, error(")
 mk("W3", "shapeNoiseWrapper", "testWrappers", "template step never runs", "useTemplate = ~isempty(options.F0) || all(~isnan(options.F0Band));", "useTemplate = false;")
 mk("W4", "shapeNoiseWrapper", "testWrappers", "default K restored (1)", "options.K double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []", "options.K double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = 1")
 mk("W5", "cetNoiseWrapper", "testWrappers", "default K restored (2)", "options.K double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []", "options.K double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = 2")
 mk("S7", "shapeNoise", SB, "no-draw warning no longer names MaxDraw", "of the target in MaxDraw = %d tries", "of the target in %d tries")
 mk("B1", "summariseShapeNoiseBatch", SB, "flipped inequality in the test", "isLow = pValue < options.Level;", "isLow = pValue > options.Level;")
 mk("B2", "summariseShapeNoiseBatch", SB, "minimum number of channels exclusive", "tested = nChannels >= options.MinChannels;", "tested = nChannels > options.MinChannels;")
 mk("B3", "summariseShapeNoiseBatch", SB, "wrong default minimum (60)", "mustBePositive} = 90", "mustBePositive} = 60")
 mk("B4", "summariseShapeNoiseBatch", SB, "wrong default level (0.05)", "mustBeLessThan(options.Level, 1)} = 0.01", "mustBeLessThan(options.Level, 1)} = 0.05")
 mk("B5", "summariseShapeNoiseBatch", SB, "wrong default margin (0.99)", "mustBeLessThan(options.Margin, 1)} = 0.95", "mustBeLessThan(options.Margin, 1)} = 0.99")
 mk("B16", "summariseShapeNoiseBatch", SB, "wrong default standard (0.95)", "mustBeLessThan(options.Standard, 1)} = 0.99", "mustBeLessThan(options.Standard, 1)} = 0.95")
 mk("B17", "summariseShapeNoiseBatch", SB, "warning tests the standard, not the margin", "pValue = binomialLowerTail(nAccepted, nChannels, options.Margin);", "pValue = binomialLowerTail(nAccepted, nChannels, options.Standard);")
 mk("B18", "summariseShapeNoiseBatch", SB, "reported standard result uses the margin", "pValueStandard = binomialLowerTail(nAccepted, nChannels, options.Standard);", "pValueStandard = binomialLowerTail(nAccepted, nChannels, options.Margin);")
 mk("B19", "summariseShapeNoiseBatch", SB, "belowStandard flag inverted", """belowStandard"", pValueStandard < options.Level", """belowStandard"", pValueStandard >= options.Level")
 mk("B6", "summariseShapeNoiseBatch", SB, "lower tail stops one count short", "i = (0:x)';", "i = (0:max(x - 1, 0))';")
 mk("B7", "summariseShapeNoiseBatch", SB, "tries percentile is the median", "sorted(ceil(0.9 * nChannels))", "sorted(ceil(0.5 * nChannels))")
 mk("B8", "summariseShapeNoiseBatch", SB, "only the first channel of a call counts", "accepted = [accepted; a(:)];", "accepted = [accepted; a(1)];")
 mk("B9", "summariseShapeNoiseBatch", SB, "wrong warning identifier", "warning(""shapeNoise:lowAcceptance"",", "warning(""shapeNoise:mutated"",")
 mk("B10", "summariseShapeNoiseBatch", SB, "fractional tries accepted", "any(d ~= round(d))", "false")
 mk("B11", "summariseShapeNoiseBatch", SB, "channel-count mismatch not detected", "if numel(a) ~= numel(d), error(", "if false, error(")
 mk("B12", "summariseShapeNoiseBatch", SB, "empty batch not detected", "if nInfos == 0, error(", "if false, error(")
 mk("B13", "summariseShapeNoiseBatch", SB, "untested batch reported as ok", "decision = ""tooFew"";", "decision = ""ok"";")
 mk("B14", "summariseShapeNoiseBatch", SB, "wrong binomial term", "(n - i) * log1p(-q)", "(n - i) * log(q)")
 mk("B15", "summariseShapeNoiseBatch", SB, "per-draw acceptance ignores the draws made", "nAccepted / sum(tries)", "nAccepted / nChannels")
 mk("X1", "noiseXu", NX, "exponent halved (d = alpha/4)", "d = alpha / 2;", "d = alpha / 4;")
 mk("X2", "noiseXu", NX, "Phi ignored in the recursion", "h = phi .* (k + d - 1) ./ k;", "h = (k + d - 1) ./ k;")
 mk("X3", "noiseXu", NX, "output not rescaled to the amplitude", "noise = noise * (amplitude / currentStd);", "noise = noise * amplitude;")
 mk("X4", "noiseXu", NX, "default Phi 1 instead of 0.99", "mustBeLessThanOrEqual(options.Phi, 1)} = 0.99", "mustBeLessThanOrEqual(options.Phi, 1)} = 1")
 mk("A1", "estimateWhittleAlpha", WA, "differencing order one too low", "kUsed = min(max(floor((info.alphaDFA + 1) / 2), 0), options.MaxK);", "kUsed = min(max(floor((info.alphaDFA - 1) / 2), 0), options.MaxK);")
 mk("A2", "estimateWhittleAlpha", WA, "Hann taper dropped", "w = 0.5 * (1 - cos(2*pi*(0:N-1)' / N));", "w = ones(N, 1);")
 mk("A3", "estimateWhittleAlpha", WA, "periodogram not recoloured", "Ib = I(m) ./ (2 - 2*cos(lam)).^kUsed;", "Ib = I(m);")
 mk("A4", "estimateWhittleAlpha", WA, "estimate on a bound not an error", "if flag ~= 1 || alphaEst - lo < 1e-3 || hi - alphaEst < 1e-3", "if flag ~= 1")
 mk("A5", "estimateWhittleAlpha", WA, "too few bins allowed (3, not 10)", "if info.nBins < 10", "if info.nBins < 3")
 mk("A6", "estimateWhittleAlpha", WA, "DFA exponent not converted (2H, not 2H - 1)", "a = 2 * p(1) - 1;", "a = 2 * p(1);")
 mk("P1", "estimateWhittleAlpha", WA, "peak check never warns", "if options.PeakCheck && info.peakRatio > info.peakThreshold", "if false && info.peakRatio > info.peakThreshold")
 mk("P2", "estimateWhittleAlpha", WA, "Hann correlation factor dropped (n not x 18/35)", "nEff = nUsed * 18 / 35;", "nEff = nUsed;")
 mk("P3", "estimateWhittleAlpha", WA, "the bin counts among its own neighbours (gap 0)", "W = 12;  gap = 2;", "W = 12;  gap = 0;")
 mk("P4", "estimateWhittleAlpha", WA, "threshold ignores the false-alarm rate", "pBin = falseAlarm / nnz(tested);", "pBin = 1e-2;")
 mk("P5", "shapeNoise", SN, "peak check repeated in every draw", "PeakCheck=false) - aRef;", "PeakCheck=true) - aRef;")
 mk("P6", "shapeNoise", SN, "peak hit not recorded", "flag = info.peakRatio > info.peakThreshold;  hz = NaN;", "flag = false;  hz = NaN;")
 mk("I1", "estimateIRASA", IR, "IRASA reverts to the pre-2.0 resampling bug", "isIrasa = resampling == ""irasa"";", "isIrasa = false;")
 mk("I2", "estimateIRASA", IR, "mean instead of median over Hset", "pFractal = median(geoMeans, 2, 'omitnan');", "pFractal = mean(geoMeans, 2, 'omitnan');")
 mk("I3", "estimateIRASA", IR, "band check fHigh*max(Hset) < fs/2 removed", "if options.Resampling == ""irasa"" && fHigh * max(options.Hset) >= fs / 2", "if false")
 mk("I4", "estimateIRASA", IR, "prewhitened copies not recoloured", "geoMeans(:, hi) = geoMeans(:, hi) ./ gainH;", "geoMeans(:, hi) = geoMeans(:, hi);")
 mk("I5", "estimateIRASA", IR, "noise variance halved", "noiseVar = trapz(fVec(validMask), pFractal(validMask)) * 2;", "noiseVar = trapz(fVec(validMask), pFractal(validMask));")
 mk("K1", "checkSpectralKnee", CK, "strong threshold -2.5 instead of -10", "if deltaAIC < -10", "if deltaAIC < -2.5")
 mk("V1", "convertScalingExponent", CS, "fGn H as alpha = 2H (not 2H - 1)", "a = 2 * x - 1 + 2 * fBm;", "a = 2 * x + 2 * fBm;")
 mk("V2", "convertScalingExponent", CS, "H ignores the class on output", "case ""H"",     y = (a + 1 - 2 * fBm) / 2;", "case ""H"",     y = (a + 1) / 2;")
 mk("V3", "convertScalingExponent", CS, "D ignores the class on input", "a = 3 - 2 * x + 2 * fBm;", "a = 3 - 2 * x;")
 mk("V4", "convertScalingExponent", CS, "DFA not converted (alpha = 2a)", "case ""DFA"",   a = 2 * x - 1;", "case ""DFA"",   a = 2 * x;")
 mk("V5", "convertScalingExponent", CS, "no class required for H or D", "if needsClass && options.Class == """"", "if false")
 mk("V6", "convertScalingExponent", CS, "class range not checked", "if any(a <= lo | a >= hi, ""all"")", "if false")
 mk("V7", "convertScalingExponent", CS, "D ignores the class on output", "case ""D"",     y = (3 - a + 2 * fBm) / 2;", "case ""D"",     y = (3 - a) / 2;")
 mk("V8", "convertScalingExponent", CS, "d taken as alpha (not alpha/2)", "case ""d"",     a = 2 * x;", "case ""d"",     a = x;")
 mk("K2", "checkSpectralKnee", CK, "knee model penalised as 2 parameters", "aicKnee = nWide * log(bestRss / nWide) + 2 * 3;", "aicKnee = nWide * log(bestRss / nWide) + 2 * 2;")
 mk("K3", "checkSpectralKnee", CK, "kneeInBand ignores the lower band edge", "result.kneeInBand = kneeHz >= fBandLow && kneeHz <= fHigh;", "result.kneeInBand = kneeHz <= fHigh;")
 mk("K4", "checkSpectralKnee", CK, "too-short residual allowed (8, not 80)", "if nSamples < 80", "if nSamples < 8")
 mk("K5", "checkSpectralKnee", CK, "Lorentzian exponent doubled", "log10(10 ^ (2 * th(2)) + f .^ 2);", "log10(10 ^ (2 * th(2)) + f .^ 4);")
 ];
end
