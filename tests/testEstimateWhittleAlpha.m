classdef testEstimateWhittleAlpha < matlab.unittest.TestCase
%TESTESTIMATEWHITTLEALPHA Unit tests for estimateWhittleAlpha.m.
%
% Toolbox: fractalnoise
% Covers: source-equivalence (alphaPWLW_v001), known-alpha recovery on noiseXu output and (if the
% validation tree is present) on exact integrated-fGn references, oracle K, fail-loud cases, and a
% toolbox-only-path check for this function (kept here so testToolboxSelfContained is not edited).

    properties (Constant)
        ToolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox")
        ExactRefDir = fullfile(fileparts(mfilename("fullpath")), "..", "validation", "exactReference")
    end

    methods (TestClassSetup)
        function addPaths(testCase)
            addpath(testCase.ToolboxDir);  testCase.addTeardown(@() rmpath(testCase.ToolboxDir));
            addLegacySource(testCase);   % FRACTALNOISE_LEGACY_SRC
        end
    end

    methods (Test)
        function testMatchesAlphaPWLWv001(testCase)
            testCase.assumeTrue(exist("alphaPWLW_v001", "file") == 2, "alphaPWLW_v001 not on path (set FRACTALNOISE_LEGACY_SRC); skipping.");
            rng(11);  x = noiseXu(2400, 4.5, 1, 240);
            [a1, k1] = estimateWhittleAlpha(x, 240, 1, 20);
            [a2, k2] = alphaPWLW_v001(x, 240, [1 20]);
            testCase.verifyEqual(a1, a2, "AbsTol", 1e-10);
            testCase.verifyEqual(k1, k2);
        end

        function testRecoversAlphaOnNoiseXu(testCase)
            % noiseXu (Phi 0.99) band-effective alpha over 1-20 Hz at fs 240 is within ~0.05 of nominal for these alphas.
            rng(3);
            for a = [2.5 4 5.5]
                est = arrayfun(@(~) estimateWhittleAlpha(noiseXu(3000, a, 1, 240), 240, 1, 20), 1:20);
                testCase.verifyLessThan(abs(median(est) - a), 0.15, sprintf("alpha %.1f: median %.3f", a, median(est)));
            end
        end

        function testRecoversAlphaOnExactReference(testCase)
            testCase.assumeTrue(isfolder(testCase.ExactRefDir), "validation/exactReference not present; skipping.");
            addpath(testCase.ExactRefDir);  testCase.addTeardown(@() rmpath(testCase.ExactRefDir));
            rng(5);
            for a = [4 6]
                X = exactRefSeries_v001(3000, a, 20);
                est = arrayfun(@(j) estimateWhittleAlpha(X(:, j), 240, 1, 20), 1:20);
                testCase.verifyLessThan(abs(median(est) - a), 0.15, sprintf("alpha %.1f: median %.3f", a, median(est)));
            end
        end

        function testOracleKIsUsed(testCase)
            rng(9);  x = noiseXu(2000, 3.5, 1, 100);
            [~, k] = estimateWhittleAlpha(x, 100, 1, 20, K=2);
            testCase.verifyEqual(k, 2);
        end

        function testErrorsOnConstantInput(testCase)
            testCase.verifyError(@() estimateWhittleAlpha(zeros(1000, 1), 240, 1, 20), "estimateWhittleAlpha:constant");
        end

        function testErrorsOnNaN(testCase)
            x = randn(1000, 1);  x(10) = NaN;
            testCase.verifyError(@() estimateWhittleAlpha(x, 240, 1, 20), "estimateWhittleAlpha:nonFinite");
        end

        function testErrorsOnBadBand(testCase)
            testCase.verifyError(@() estimateWhittleAlpha(randn(1000, 1), 40, 1, 20), "estimateWhittleAlpha:band");
        end

        function testErrorsOnShortSeries(testCase)
            testCase.verifyError(@() estimateWhittleAlpha(randn(80, 1), 240, 1, 100), "estimateWhittleAlpha:short");
        end

        function testWorksWithOnlyToolboxOnPath(testCase)
            originalPath = path();  testCase.addTeardown(@() path(originalPath));
            restoredefaultpath();  rehash();  addpath(testCase.ToolboxDir);
            rng(1);  a = estimateWhittleAlpha(noiseXu(1500, 3.0, 1, 100), 100, 1, 20);
            testCase.verifyTrue(isfinite(a));
        end
    end
end
