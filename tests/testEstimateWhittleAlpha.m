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

        function testDfaPreEstimateKnownAnswers(testCase)
            % The DFA-3 pre-estimate is returned as a spectral exponent (2H - 1): about 0 for white noise, about 2 for alpha-2 noise.
            rng(2);  [~, ~, iw] = estimateWhittleAlpha(randn(3000, 1), 240, 1, 20);
            rng(3);  [~, ~, ib] = estimateWhittleAlpha(noiseXu(3000, 2, 1, 240, Phi=1), 240, 1, 20);
            testCase.verifyEqual(iw.alphaDFA, 0, "AbsTol", 0.25);
            testCase.verifyEqual(ib.alphaDFA, 2, "AbsTol", 0.35);
        end

        function testDifferencingOrderAndAccuracyOnExactReferences(testCase)
            % Known truth: k = floor((alpha + 1) / 2) for every series, and the median estimate within 0.1 of alpha.
            testCase.assumeTrue(isfolder(testCase.ExactRefDir), "validation/exactReference not present; skipping.");
            addpath(testCase.ExactRefDir);  testCase.addTeardown(@() rmpath(testCase.ExactRefDir));
            for a = [0.5 2.5 4.5 6.5]
                X = exactRefSeries_v001(3000, a, 20);  est = zeros(20, 1);  k = est;
                for r = 1:20, [est(r), k(r)] = estimateWhittleAlpha(X(:, r), 240, 1, 20); end
                testCase.verifyEqual(k, repmat(floor((a + 1) / 2), 20, 1), sprintf("Differencing order at alpha %.1f.", a));
                testCase.verifyEqual(median(est), a, "AbsTol", 0.1, sprintf("Median estimate at alpha %.1f.", a));
            end
        end

        function testEstimateOnSearchBoundIsAnError(testCase)
            rng(5);  x = noiseXu(3000, 4, 1, 240);
            testCase.verifyError(@() estimateWhittleAlpha(x, 240, 1, 20, AlphaBounds=[-2 2]), "estimateWhittleAlpha:bound");
        end

        function testFewerThanTenBinsIsAnError(testCase)
            rng(4);  x = noiseXu(600, 2, 1, 240);   % 0.4 Hz resolution: 8 bins in [1 4] Hz, 12 in [1 6] Hz
            testCase.verifyError(@() estimateWhittleAlpha(x, 240, 1, 4), "estimateWhittleAlpha:band");
            [~, ~, info] = estimateWhittleAlpha(x, 240, 1, 6);
            testCase.verifyEqual(info.nBins, 12);
        end

        function testSpreadMatchesHannVarianceFactor(testCase)
            % Known answer: with the Hann taper the SD over exact references is sqrt(35/18) = 1.39 times the untapered
            % local Whittle SD for the band (Velasco 1999, p. 101); measured 1.40 here, 0.99 without the taper.
            testCase.assumeTrue(isfolder(testCase.ExactRefDir), "validation/exactReference not present; skipping.");
            addpath(testCase.ExactRefDir);  testCase.addTeardown(@() rmpath(testCase.ExactRefDir));
            X = exactRefSeries_v001(3000, 2.5, 100);  e = zeros(100, 1);
            for r = 1:100, e(r) = estimateWhittleAlpha(X(:, r), 240, 1, 20); end
            M = 3000 - 1;  f = (0:M - 1)' * 240 / M;  f = f(f >= 1 & f <= 20);  nu = log(f) - mean(log(f));
            ratio = std(e) * sqrt(sum(nu .^ 2));
            testCase.verifyGreaterThan(ratio, 1.2);
            testCase.verifyLessThan(ratio, 1.65);
        end

        function testHannTaperProtectsAgainstUnderDifferencing(testCase)
            % Known answer (validation/compareHannTaper_v001.m): with k forced one too low the Hann-tapered estimate stays within 0.03
            % of the truth (median of 50); without the taper it is 0.43-0.48 too low. Here 20 series, tolerance 0.1.
            testCase.assumeTrue(isfolder(testCase.ExactRefDir), "validation/exactReference not present; skipping.");
            addpath(testCase.ExactRefDir);  testCase.addTeardown(@() rmpath(testCase.ExactRefDir));
            for a = [2.5 4.5]
                X = exactRefSeries_v001(3000, a, 20);  e = zeros(20, 1);
                for r = 1:20, e(r) = estimateWhittleAlpha(X(:, r), 240, 1, 20, K=floor((a + 1) / 2) - 1); end
                testCase.verifyEqual(median(e), a, "AbsTol", 0.1, sprintf("k one too low at alpha %.1f.", a));
            end
        end
    end
end
