classdef testCheckSpectralKnee < matlab.unittest.TestCase
%TESTCHECKSPECTRALKNEE Unit tests for checkSpectralKnee.m.
%
% Toolbox: fractalnoise
% Covers: source-equivalence (checkLorentzianKnee_v001's per-trial core,
% reconstructed inline for comparison since the source is a script, not
% a callable function), default-starts behaviour, and error paths.

    methods (TestClassSetup)
        function addToolboxAndLegacyPaths(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
            addLegacySource(testCase);   % FRACTALNOISE_LEGACY_SRC
        end
    end

    methods (Test)
        function testMatchesSourcePerTrialCoreGivenSameStarts(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2 ...
                && exist("iraAlphaSigma_v003", "file") == 2, ...
                "Legacy source functions not on path (set FRACTALNOISE_LEGACY_SRC); skipping source-equivalence test.");

            fs = 133;
            fitBandLow = 1.0;
            wideFitLow = 0.15;
            fHigh = fs / 2.2;
            hset = 1.1:0.05:1.9;
            phiXu = 0.99;
            fKneeXu = (1 - phiXu) * fs / (2 * pi);

            rng(5);
            x = generateCustomNoise_v004(900, 2.5, 3.0, fs);
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("estimateIRASA:legacyResampling"));   % the source script is the legacy form
            sig = x - mean(x);

            lorentzLog10 = @(th, f) th(1) - (th(3) / 2) .* log10(10 ^ (2 * th(2)) + f .^ 2);
            fminOpts = optimset("Display", "off", "TolFun", 1e-10, "TolX", 1e-10, "MaxFunEvals", 10000);

            [~, ~, pFrac, fVec] = iraAlphaSigma_v003(sig, fs, wideFitLow, fHigh, hset);
            fVec = fVec(:); pFrac = pFrac(:);

            maskS = fVec >= fitBandLow & fVec <= fHigh & ~isnan(pFrac) & pFrac > 0;
            pCS = polyfit(log10(fVec(maskS)), log10(pFrac(maskS)), 1);
            alphaStdOld = -pCS(1);

            maskW = fVec >= wideFitLow & fVec <= fHigh & ~isnan(pFrac) & pFrac > 0;
            fW = fVec(maskW); logPW = log10(pFrac(maskW)); nW = numel(fW);
            pCW = polyfit(log10(fW), logPW, 1);
            bW = pCW(2); alphaW = -pCW(1);
            rssF = sum((logPW - (bW - alphaW * log10(fW))) .^ 2);
            aicF = nW * log(rssF / nW) + 2 * 2;

            sseFn = @(th) sum((logPW - lorentzLog10(th, fW)) .^ 2);
            starts = [bW, log10(fKneeXu), alphaW; bW, log10(wideFitLow), alphaW; bW, log10(fitBandLow), alphaW];
            bestRss = Inf; thBest = starts(1, :);
            for si = 1:size(starts, 1)
                [thTry, rssTry] = fminsearch(sseFn, starts(si, :), fminOpts);
                if rssTry < bestRss
                    bestRss = rssTry; thBest = thTry;
                end
            end
            kneeHzOld = 10 ^ thBest(2);
            daicOld = nW * log(bestRss / nW) + 2 * 3 - aicF;

            result = checkSpectralKnee(x, fs, FitBandLow=fitBandLow, WideBandLow=wideFitLow, ...
                Hset=hset, ExtraKneeStarts=fKneeXu, Resampling="legacy");

            testCase.verifyEqual(result.alphaStd, alphaStdOld, "AbsTol", 1e-9);
            testCase.verifyEqual(result.alphaWide, alphaW, "AbsTol", 1e-9);
            testCase.verifyEqual(result.kneeHz, kneeHzOld, "AbsTol", 1e-6);
            testCase.verifyEqual(result.deltaAIC, daicOld, "AbsTol", 1e-6);
        end

        function testDefaultArgsCallRuns(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2, ...
                "generateCustomNoise_v004 not on path (set FRACTALNOISE_LEGACY_SRC); skipping.");
            rng(5);
            x = generateCustomNoise_v004(900, 2.5, 3.0, 133);
            result = checkSpectralKnee(x, 133);
            testCase.verifyEqual(result.nSamples, 900);
            testCase.verifyTrue(ismember(result.verdict, ["strong", "marginal", "none"]));
            testCase.verifyTrue(isfinite(result.alphaStd));
            testCase.verifyTrue(isfinite(result.kneeHz));
        end

        function testDefaultRunsOnToolboxNoiseWithIrasaBand(testCase)
            % Needs no legacy code: default call on noiseXu output, band capped for corrected IRASA.
            rng(21);  x = noiseXu(1500, 2.0, 1.0, 100);
            result = testCase.verifyWarningFree(@() checkSpectralKnee(x, 100));
            testCase.verifyEqual(result.fHigh, 0.95 * 100 / (2 * 1.9), "RelTol", 1e-12);
            testCase.verifyTrue(ismember(result.verdict, ["strong", "marginal", "none"]));
            testCase.verifyTrue(isfinite(result.alphaStd) && isfinite(result.kneeHz));
        end

        function testBandAboveIrasaLimitErrors(testCase)
            rng(22);  x = noiseXu(1500, 2.0, 1.0, 100);
            testCase.verifyError(@() checkSpectralKnee(x, 100, FitBandHigh=45), "estimateIRASA:irasaBand");
        end

        function testLegacyResamplingKeepsOldBandAndWarns(testCase)
            rng(23);  x = noiseXu(1500, 2.0, 1.0, 100);
            result = testCase.verifyWarning(@() checkSpectralKnee(x, 100, Resampling="legacy"), "estimateIRASA:legacyResampling");
            testCase.verifyEqual(result.fHigh, 100 / 2.2, "RelTol", 1e-12);
        end

        function testTooShortResidualErrors(testCase)
            testCase.verifyError(@() checkSpectralKnee(randn(50, 1), 133), ...
                "checkSpectralKnee:TooShort");
        end

        function testKneeInBandLogicConsistentWithBounds(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2, ...
                "generateCustomNoise_v004 not on path (set FRACTALNOISE_LEGACY_SRC); skipping.");
            rng(8);
            x = generateCustomNoise_v004(1200, 2.0, 2.0, 100);
            result = checkSpectralKnee(x, 100);
            expected = result.kneeHz >= result.fBandLow && result.kneeHz <= result.fHigh;
            testCase.verifyEqual(result.kneeInBand, expected);
        end

        function testBreakFreeNoiseIsPurePenalty(testCase)
            % No knee: the knee model fits no better, so deltaAIC is exactly the extra parameter's penalty, +2.
            rng(8);  r = checkSpectralKnee(noiseXu(4000, 3.8, 3, 133, Phi=1), 133);
            testCase.verifyEqual(string(r.verdict), "none");
            testCase.verifyEqual(r.deltaAIC, 2, "AbsTol", 1e-9);
        end

        function testKneeLocatedBelowAndInsideBand(testCase)
            % Exact Lorentzian spectra (alpha 4.5): knee 0.28 Hz is below the 1 Hz band edge, knee 3 Hz is inside it.
            ff = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox", "external", "fftnoise");  addpath(ff);  testCase.addTeardown(@() rmpath(ff));
            rng(1);  below = checkSpectralKnee(lorentzianSeries(4000, 133, 4.5, 0.28), 133, WideBandLow=0.05);
            rng(2);  inside = checkSpectralKnee(lorentzianSeries(4000, 133, 4.5, 3.0), 133);
            testCase.verifyFalse(below.kneeInBand);
            testCase.verifyEqual(below.kneeHz, 0.28, "AbsTol", 0.1);
            testCase.verifyTrue(inside.kneeInBand);
            testCase.verifyEqual(inside.kneeHz, 3.0, "AbsTol", 0.5);
            testCase.verifyEqual(string(inside.verdict), "strong");
        end

        function testMarginalVerdictBand(testCase)
            % A weak sub-floor knee in a short record: deltaAIC between -10 and -2 is "marginal" (here about -5.5).
            ff = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox", "external", "fftnoise");  addpath(ff);  testCase.addTeardown(@() rmpath(ff));
            fs = 133;  N = 500;  f = (0:N - 1)' * fs / N;  f(f > fs / 2) = f(f > fs / 2) - fs;
            rng(300);  lz = fftnoise((0.5 ^ 2 + f .^ 2) .^ (-3 / 4));  lz = lz / std(lz);
            rng(600);  bf = noiseXu(N, 1.5, 1, fs, Phi=1);
            r = checkSpectralKnee(0.6 * lz + bf, fs);
            testCase.verifyGreaterThan(r.deltaAIC, -10);
            testCase.verifyLessThan(r.deltaAIC, -2);
            testCase.verifyEqual(string(r.verdict), "marginal");
        end
    end
end

function x = lorentzianSeries(N, fs, alpha, kneeHz)
% Series with an exact Lorentzian spectrum: magnitude (knee^2 + f^2)^(-alpha/4) at every FFT bin, random phases (fftnoise).
f = (0:N - 1)' * fs / N;
f(f > fs / 2) = f(f > fs / 2) - fs;
x = fftnoise((kneeHz ^ 2 + f .^ 2) .^ (-alpha / 4));
x = (x - mean(x)) / std(x);
end
