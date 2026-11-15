classdef testShapeNoise < matlab.unittest.TestCase
%TESTSHAPENOISE Unit tests for shapeNoise.m.
%
% Toolbox: fractalnoise
% Covers, for the default end-matched surrogate (v2.0.0): source-equivalence
% with endToEndSegment_v006 + endMatchedNoise_v001, derivative-noise fidelity
% and alpha selection on steep non-periodic noise, segment bookkeeping, two
% channels, reproducibility, and error paths (zero tail, PhaseSource,
% constant, alpha size, no segment). For Method="legacy" (the 1.x
% construction): source-equivalence with shapeXu_local, all PhaseSource
% carriers, edge cases, error paths, and its warning.

    methods (TestClassSetup)
        function addToolboxAndLegacyPaths(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
            addLegacySource(testCase);   % FRACTALNOISE_LEGACY_SRC
        end
    end

    methods (Test)
        function testMatchesSourceShapeXuLocalForSameSeed(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2, ...
                "generateCustomNoise_v004 not on path (set FRACTALNOISE_LEGACY_SRC); skipping source-equivalence test.");
            rng(11);
            resid = cumsum(randn(500, 1));
            resid = resid - mean(resid);
            fs = 100;
            alpha = 3.0;

            % Reconstruct source shapeXu_local inline (private local
            % function, not independently callable).
            rng(99);
            sigOld = std(resid, 0, 1);
            eEmpOld = abs(fft(resid));
            baseOld = generateCustomNoise_v004(numel(resid), alpha, sigOld, fs);
            fbOld = fft(baseOld);
            foutOld = eEmpOld .* exp(1i * angle(fbOld));
            nOld = real(ifft(foutOld));
            sNOld = std(nOld, 0, 1);
            if sNOld > 0
                nOld = nOld * (sigOld / sNOld);
            end

            rng(99);
            nNew = shapeNoise(resid, fs, alpha, Method="legacy");

            testCase.verifyEqual(nOld, nNew);
        end

        function testWhiteCarrierPreservesStd(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            rng(21);
            resid = cumsum(randn(300, 1));
            n = shapeNoise(resid, 100, 2.0, PhaseSource="white", Method="legacy");
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testCustomFunctionHandleCarrier(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            rng(31);
            resid = cumsum(randn(300, 1));
            customCarrier = @(M, fs, a) randn(M, 1) * 2;
            n = shapeNoise(resid, 100, 2.0, PhaseSource=customCarrier, Method="legacy");
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testZeroVarianceResidualReturnsZeros(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            n = shapeNoise(zeros(50, 1), 100, 2.0, Method="legacy");
            testCase.verifyTrue(all(n == 0));
        end

        function testUnknownPhaseSourceErrors(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            resid = randn(100, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource="bogus", Method="legacy"), ...
                "shapeNoise:UnknownPhaseSource");
        end

        function testInvalidCarrierLengthErrors(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            resid = randn(100, 1);
            badCarrier = @(M, fs, a) randn(M + 1, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource=badCarrier, Method="legacy"), ...
                "shapeNoise:InvalidCarrier");
        end

        function testPinknoiseWarnsOnAlphaMismatchWhenAvailable(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            testCase.assumeTrue(exist("pinknoise", "file") == 2 || exist("pinknoise", "builtin") ~= 0, ...
                "pinknoise not available on this MATLAB installation; skipping.");
            resid = randn(200, 1);
            testCase.verifyWarning(@() shapeNoise(resid, 100, 3.0, PhaseSource="pinknoise", Method="legacy"), ...
                "shapeNoise:PinknoiseAlphaMismatch");
        end

        function testDspInRangeAlphaRunsCleanly(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            testCase.assumeTrue(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise (DSP System Toolbox) not installed; skipping.");
            rng(41);
            resid = cumsum(randn(400, 1));
            resid = resid - mean(resid);
            n = testCase.verifyWarningFree(@() shapeNoise(resid, 100, 1.5, PhaseSource="dsp", Method="legacy"));
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testDspOutOfRangeAlphaClampsWithWarning(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            testCase.assumeTrue(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise (DSP System Toolbox) not installed; skipping.");
            rng(42);
            resid = cumsum(randn(400, 1));
            resid = resid - mean(resid);
            testCase.verifyWarning(@() shapeNoise(resid, 100, 4.0, PhaseSource="dsp", Method="legacy"), ...
                "shapeNoise:DspAlphaClamped");
            n = shapeNoise(resid, 100, 4.0, PhaseSource="dsp", Method="legacy");
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testDspUnavailableErrorsCleanly(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            testCase.assumeFalse(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise is installed on this machine; unavailable-path test not applicable.");
            resid = randn(200, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource="dsp", Method="legacy"), ...
                "shapeNoise:DspUnavailable");
        end

        function testMatchesEndMatchedSourceForSameSeed(testCase)
            % Default ("detect"): same segment as endToEndSegment_v006 and same draws as endMatchedNoise_v001.
            addPreRegSrc(testCase);
            rng(5);  R = [noiseXu(1500, 4.0, 1, 133), noiseXu(1500, 3.5, 1, 133)];
            [a0, a1] = endToEndSegment_v006(R, Rule="detect");
            sMaj = R(a0(1):a1(1), 1);  sMin = R(a0(2):a1(2), 2);
            aSeg = [alphaPWLW_v001(sMaj, 133, [1 20]), alphaPWLW_v001(sMin, 133, [1 20])];
            rng(77);  [oMaj, oMin] = endMatchedNoise_v001(sMaj, sMin, 133, aSeg(1), aSeg(2));
            rng(77);  [n, info] = shapeNoise(R, 133, aSeg, MinLength=5);
            testCase.verifyEqual([info.i0; info.i1], [a0; a1]);
            testCase.verifyEqual(n, [oMaj, oMin], "AbsTol", 1e-12);
        end

        function testSegmentRulesMatchSource(testCase)
            % All three rules give endToEndSegment_v006's segment; "detect" levels equal chi2inv's.
            addPreRegSrc(testCase);
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                ["shapeNoise:noAcceptedDraw", "shapeNoise:estimateFailed"]));   % only the segment is under test here
            rng(6);  R = [noiseXu(1200, 4.5, 1, 240), noiseXu(1200, 3.0, 1, 240)];
            args = {{"detect", {}, {}}, {"budget", {"SegmentBudget", 0.01}, {"Budget", 0.01}}, {"tol", {}, {"Tol", 1}}};
            for k = 1:numel(args)
                [rule, mine, theirs] = args{k}{:};
                [a0, a1, src] = endToEndSegment_v006(R, "Rule", rule, theirs{:});
                [~, info] = shapeNoise(R, 240, 3, "SegmentRule", rule, mine{:}, MinLength=5, MaxDraw=1);
                testCase.verifyEqual([info.i0; info.i1], [a0; a1], "rule " + rule);
                testCase.verifyEqual([info.cJump, info.cSlip, info.cSlipMax], [src.cJump, src.cSlip, src.cSlipMax], ...
                    "RelTol", 1e-12, "rule " + rule);
                if rule == "detect"
                    testCase.verifyEqual([info.qJ, info.qS], [src.qJ, src.qS], "RelTol", 1e-10);
                end
            end
        end

        function testBudgetRuleNeedsBudget(testCase)
            rng(16);  x = noiseXu(800, 3.0, 1, 100);
            testCase.verifyError(@() shapeNoise(x, 100, [], SegmentRule="budget"), "shapeNoise:segmentBudget");
        end

        function testDerivativeNoiseAndAlphaHeldOnSteepNoise(testCase)
            % The v2.0.0 acceptance criterion: derivative noise of the segment AND alpha, on steep non-periodic input.
            fs = 133;  nRep = 20;  r1 = nan(nRep, 1);  r2 = nan(nRep, 1);  dA = nan(nRep, 1);
            for r = 1:nRep
                rng(900 + r);  x = noiseXu(3000, 4.0, 1, fs);
                [n, info] = shapeNoise(x, fs, []);
                seg = x(info.i0:info.i1);
                r1(r) = std(diff(n)) / std(diff(seg));  r2(r) = std(diff(n, 2)) / std(diff(seg, 2));
                dA(r) = estimateWhittleAlpha(n, fs, 1, 20) - estimateWhittleAlpha(seg, fs, 1, 20);
            end
            testCase.verifyLessThan(abs(median(r1) - 1), 0.05, sprintf("velocity ratio median %.3f", median(r1)));
            testCase.verifyLessThan(abs(median(r2) - 1), 0.05, sprintf("acceleration ratio median %.3f", median(r2)));
            testCase.verifyLessThanOrEqual(max(abs(dA)), 0.3 + 1e-9, "every accepted draw is within AlphaTol");
        end

        function testSegmentBookkeeping(testCase)
            rng(8);  x = noiseXu(2000, 3.0, 1, 100);
            [n, info] = shapeNoise(x, 100, []);
            testCase.verifySize(n, [info.L, 1]);
            testCase.verifyEqual(info.i1 - info.i0 + 1, info.L);
            testCase.verifyEqual(info.Lfrac, info.L / 2000, "AbsTol", 1e-15);
            testCase.verifyEqual(info.SegmentRule, "detect");
            testCase.verifyLessThanOrEqual(info.cJump, info.qJ * (1 + 1e-12));
            testCase.verifyLessThanOrEqual(info.cSlipMax, info.qS * (1 + 1e-12));
            testCase.verifyEqual(std(n), std(x(info.i0:info.i1)), "RelTol", 1e-9);
            testCase.verifyTrue(info.accepted);
        end

        function testTwoChannelsShareLength(testCase)
            rng(9);  X = [noiseXu(1800, 4.5, 1, 133), noiseXu(1800, 3.0, 2, 133)];
            [n, info] = shapeNoise(X, 133, []);
            testCase.verifySize(n, [info.L, 2]);
            testCase.verifyEqual(info.i1 - info.i0, [1 1] * (info.L - 1));
        end

        function testReproducibleWithSeed(testCase)
            rng(10);  x = noiseXu(1200, 3.5, 1, 133);
            rng(3);  n1 = shapeNoise(x, 133, []);
            rng(3);  n2 = shapeNoise(x, 133, []);
            testCase.verifyEqual(n1, n2);
        end

        function testUncoveredZeroTailErrors(testCase)
            rng(11);  x = noiseXu(1000, 3.0, 1, 100);  x(end-19:end) = 0;
            testCase.verifyError(@() shapeNoise(x, 100, []), "shapeNoise:zeroTail");
        end

        function testPhaseSourceIsLegacyOnly(testCase)
            rng(12);  x = noiseXu(500, 2.0, 1, 100);
            testCase.verifyError(@() shapeNoise(x, 100, [], PhaseSource="white"), "shapeNoise:phaseSourceLegacyOnly");
        end

        function testConstantChannelErrors(testCase)
            testCase.verifyError(@() shapeNoise(ones(200, 1), 100, []), "shapeNoise:constant");
        end

        function testNoEndMatchedSegmentErrors(testCase)
            % A ramp: every segment's wrap step is (L-1) mean steps, far above any rule's level.
            rng(14);  x = (1:300)' + 0.01 * randn(300, 1);
            testCase.verifyError(@() shapeNoise(x, 100, []), "shapeNoise:noSegment");
        end

        function testSegmentTooShortForAlphaErrors(testCase)
            % 80 white samples: a segment is found, but too short for estimateWhittleAlpha.
            rng(15);  x = randn(80, 1);
            testCase.verifyError(@() shapeNoise(x, 100, []), "shapeNoise:segmentAlpha");
        end

        function testAlphaSizeErrors(testCase)
            rng(13);  X = noiseXu(600, 2.0, 1, 100);  X = [X, X(end:-1:1)];
            testCase.verifyError(@() shapeNoise(X, 100, [2 2 2]), "shapeNoise:alphaSize");
        end

        function testNoAcceptedDrawWarnsAndReturnsNaN(testCase)
            rng(14);  x = noiseXu(1500, 3.0, 1, 133);
            [n, info] = testCase.verifyWarning(@() shapeNoise(x, 133, 6.5, MaxDraw=3), "shapeNoise:noAcceptedDraw");
            testCase.verifyTrue(all(isnan(n)));
            testCase.verifyFalse(info.accepted);
            testCase.verifyEqual(info.nDraw, 3);
        end

        function testLegacyWarns(testCase)
            rng(15);  x = cumsum(randn(300, 1));
            testCase.verifyWarning(@() shapeNoise(x, 100, 2.0, Method="legacy"), "shapeNoise:legacyConstruction");
        end
    end
end

function addPreRegSrc(testCase)
% FRACTALNOISE_LEGACY_SRC holds the pre-reg src/functions; endToEndSegment_v006 lives one level up in src/.
d = string(getenv("FRACTALNOISE_LEGACY_SRC"));
if d ~= "" && isfile(fullfile(d, "..", "endToEndSegment_v006.m"))
    up = char(java.io.File(fullfile(d, "..")).getCanonicalPath());
    addpath(up);  testCase.addTeardown(@() rmpath(up));
end
testCase.assumeTrue(exist("endToEndSegment_v006", "file") == 2 && exist("endMatchedNoise_v001", "file") == 2 ...
    && exist("alphaPWLW_v001", "file") == 2 && exist("chi2inv", "file") == 2, ...
    "endToEndSegment_v006 / endMatchedNoise_v001 / alphaPWLW_v001 / chi2inv not on path (set FRACTALNOISE_LEGACY_SRC); skipping.");
end
