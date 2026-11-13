classdef testShapeNoise < matlab.unittest.TestCase
%TESTSHAPENOISE Unit tests for shapeNoise.m.
%
% Toolbox: fractalnoise
% Covers: source-equivalence (shapeXu_local, reconstructed inline since
% it is a private local function not directly callable), all PhaseSource
% carriers, edge cases, and error paths.

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
            nNew = shapeNoise(resid, fs, alpha, Prewhiten=0);   % the source is the Prewhiten=0 construction

            testCase.verifyEqual(nOld, nNew);
        end

        function testPrewhitenZeroIsLegacyConstruction(testCase)
            rng(12);  resid = cumsum(randn(500, 1));  resid = resid - mean(resid);
            rng(98);  sig = std(resid, 0, 1);
            nOld = real(ifft(abs(fft(resid)) .* exp(1i * angle(fft(noiseXu(500, 3.0, sig, 100))))));
            nOld = nOld * (sig / std(nOld, 0, 1));
            rng(98);  nNew = shapeNoise(resid, 100, 3.0, Prewhiten=0);
            testCase.verifyEqual(nNew, nOld);
        end

        function testPrewhitenAutoRemovesShallowTail(testCase)
            % Numbers behind the Prewhiten help text (validation/checkShapeNoiseMechanism_v001.m): at alpha 4 about
            % one surrogate in five is more than 0.5 shallower than its reference with Prewhiten=0, none with "auto".
            nRep = 40;  d0 = nan(nRep, 1);  dA = nan(nRep, 1);
            for r = 1:nRep
                rng(600 + r);  ref = noiseXu(3000, 4.0, 1, 133);  aRef = estimateWhittleAlpha(ref, 133, 1, 20);
                rng(700 + r);  d0(r) = estimateWhittleAlpha(shapeNoise(ref, 133, 4.0, Prewhiten=0), 133, 1, 20) - aRef;
                rng(700 + r);  dA(r) = estimateWhittleAlpha(shapeNoise(ref, 133, 4.0), 133, 1, 20) - aRef;
            end
            testCase.verifyLessThanOrEqual(mean(dA < -0.5), 0.05, "The default left shallow surrogates at alpha 4.");
            testCase.verifyLessThan(abs(median(dA)), 0.15);
            testCase.verifyGreaterThan(mean(d0 < -0.5), 0.08, "Prewhiten=0 no longer shows the shallow tail; update the help text.");
        end

        function testPrewhitenPreservesStdAndMean(testCase)
            rng(13);  resid = 3 + cumsum(cumsum(randn(800, 1)));
            for k = 1:3
                n = shapeNoise(resid, 100, 4.0, Prewhiten=k);
                testCase.verifyEqual(numel(n), numel(resid));
                testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
                testCase.verifyEqual(mean(n), mean(resid), "AbsTol", 1e-9 * std(resid));
            end
        end

        function testPrewhitenAutoBelowAlphaOneIsZero(testCase)
            rng(14);  resid = randn(400, 1);
            rng(97);  nAuto = shapeNoise(resid, 100, 0.5);
            rng(97);  nZero = shapeNoise(resid, 100, 0.5, Prewhiten=0);
            testCase.verifyEqual(nAuto, nZero);
        end

        function testInvalidPrewhitenErrors(testCase)
            resid = randn(200, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, Prewhiten=4), "shapeNoise:InvalidPrewhiten");
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, Prewhiten="bogus"), "shapeNoise:InvalidPrewhiten");
            testCase.verifyError(@() shapeNoise(resid, 100, NaN, PhaseSource="white"), "shapeNoise:PrewhitenAutoNeedsAlpha");
            testCase.verifyError(@() shapeNoise(randn(10, 1), 100, 2.0, Prewhiten=3), "shapeNoise:TooShortForPrewhiten");
        end

        function testWhiteCarrierPreservesStd(testCase)
            rng(21);
            resid = cumsum(randn(300, 1));
            n = shapeNoise(resid, 100, 2.0, PhaseSource="white");
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testCustomFunctionHandleCarrier(testCase)
            rng(31);
            resid = cumsum(randn(300, 1));
            customCarrier = @(M, fs, a) randn(M, 1) * 2;
            n = shapeNoise(resid, 100, 2.0, PhaseSource=customCarrier);
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testZeroVarianceResidualReturnsZeros(testCase)
            n = shapeNoise(zeros(50, 1), 100, 2.0);
            testCase.verifyTrue(all(n == 0));
        end

        function testUnknownPhaseSourceErrors(testCase)
            resid = randn(100, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource="bogus"), ...
                "shapeNoise:UnknownPhaseSource");
        end

        function testInvalidCarrierLengthErrors(testCase)
            resid = randn(100, 1);
            badCarrier = @(M, fs, a) randn(M + 1, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource=badCarrier), ...
                "shapeNoise:InvalidCarrier");
        end

        function testPinknoiseWarnsOnAlphaMismatchWhenAvailable(testCase)
            testCase.assumeTrue(exist("pinknoise", "file") == 2 || exist("pinknoise", "builtin") ~= 0, ...
                "pinknoise not available on this MATLAB installation; skipping.");
            resid = randn(200, 1);
            testCase.verifyWarning(@() shapeNoise(resid, 100, 3.0, PhaseSource="pinknoise"), ...
                "shapeNoise:PinknoiseAlphaMismatch");
        end

        function testDspInRangeAlphaRunsCleanly(testCase)
            testCase.assumeTrue(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise (DSP System Toolbox) not installed; skipping.");
            rng(41);
            resid = cumsum(randn(400, 1));
            resid = resid - mean(resid);
            n = testCase.verifyWarningFree(@() shapeNoise(resid, 100, 1.5, PhaseSource="dsp"));
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testDspOutOfRangeAlphaClampsWithWarning(testCase)
            testCase.assumeTrue(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise (DSP System Toolbox) not installed; skipping.");
            rng(42);
            resid = cumsum(randn(400, 1));
            resid = resid - mean(resid);
            testCase.verifyWarning(@() shapeNoise(resid, 100, 4.0, PhaseSource="dsp"), ...
                "shapeNoise:DspAlphaClamped");
            n = shapeNoise(resid, 100, 4.0, PhaseSource="dsp");
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testDspUnavailableErrorsCleanly(testCase)
            testCase.assumeFalse(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise is installed on this machine; unavailable-path test not applicable.");
            resid = randn(200, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource="dsp"), ...
                "shapeNoise:DspUnavailable");
        end
    end
end
