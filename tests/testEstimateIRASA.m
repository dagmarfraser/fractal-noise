classdef testEstimateIRASA < matlab.unittest.TestCase
%TESTESTIMATEIRASA Unit tests for estimateIRASA.m.
%
% Toolbox: fractalnoise
% Covers: output arities, Hset override, the PSD and Prewhiten options (v2), the CompareFieldTrip
% degraded-mode path, and the ComparePMTM resampling-free comparison
% (including pRaw plumbing source-equivalence against a direct pmtm
% call). All of the above use the toolbox's own noiseXu as the test
% signal generator, so they run standalone on any machine with only the
% toolbox directory on path -- no dependency on PowerLawSimulationPreReg
% legacy code or its machine-specific paths.
%
% testMatchesSourceForSameSeed is the one deliberate exception: its job
% is specifically to confirm extraction fidelity against
% iraAlphaSigma_v003/generateCustomNoise_v004, so it stays coupled to
% the legacy source and self-filters via assumeTrue when that source
% isn't on this machine's path.

    properties (Constant)
        LegacyFunctionsDir = "/Users/d.s.fraser/Library/CloudStorage/Dropbox/Brain2Bee/PowerLawSimulationPreReg/src/functions"
        FieldTripDir = "/Users/d.s.fraser/Library/CloudStorage/Dropbox/MATLAB/fieldtrip-20251201"
    end

    methods (TestClassSetup)
        function addToolboxAndLegacyPaths(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
            if isfolder(testCase.LegacyFunctionsDir)
                addpath(testCase.LegacyFunctionsDir);
                testCase.addTeardown(@() rmpath(testCase.LegacyFunctionsDir));
            end
            if isempty(which("ft_defaults")) && isfolder(testCase.FieldTripDir)
                addpath(testCase.FieldTripDir);
                evalc("ft_defaults");
            end
        end
    end

    methods (Test)
        function testMatchesSourceForSameSeed(testCase)
            % Genuine source-equivalence test: this one must stay coupled
            % to the legacy functions, since its job is confirming
            % estimateIRASA matches iraAlphaSigma_v003 exactly. Every
            % other test below only needs *a* noise signal with known
            % alpha, which the toolbox's own noiseXu already provides --
            % see the vexed-about-hardcoded-paths discussion this
            % session, 2026-07-15.
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2 ...
                && exist("iraAlphaSigma_v003", "file") == 2, ...
                "Legacy source functions not on path; skipping source-equivalence test.");
            rng(23);
            x = generateCustomNoise_v004(3000, 3.5, 5.0, 120);
            fLow = 1.0; fHigh = 30.0; fs = 120;

            [a1, s1, p1, f1] = iraAlphaSigma_v003(x, fs, fLow, fHigh);
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("estimateIRASA:legacyResampling"));
            [a2, s2, p2, f2] = estimateIRASA(x, fs, fLow, fHigh, Resampling="legacy");   % the source is the legacy form

            testCase.verifyEqual(a1, a2);
            testCase.verifyEqual(s1, s2);
            testCase.verifyTrue(isequaln(p1, p2), ...
                "pFractal must match source including NaN positions (isequaln, not isequal).");
            testCase.verifyEqual(f1, f2);
        end

        function testTwoOutputArity(testCase)
            rng(1);
            x = noiseXu(2000, 2.0, 3.0, 100);
            [a, s] = estimateIRASA(x, 100, 1.0, 25.0);
            testCase.verifyTrue(isscalar(a));
            testCase.verifyTrue(isscalar(s));
            testCase.verifyTrue(isfinite(a));
            testCase.verifyTrue(isfinite(s));
        end

        function testCustomHsetRuns(testCase)
            rng(2);
            x = noiseXu(2000, 2.0, 3.0, 100);
            [a, s] = estimateIRASA(x, 100, 1.0, 25.0, Hset=1.2:0.1:1.8);
            testCase.verifyTrue(isfinite(a));
            testCase.verifyTrue(isfinite(s));
        end

        function testDefaultFtInfoWhenNotRequested(testCase)
            rng(4);
            x = noiseXu(1500, 2.0, 3.0, 100);
            [~, ~, ~, ~, ftInfo] = estimateIRASA(x, 100, 1.0, 25.0);
            testCase.verifyFalse(ftInfo.available);
            testCase.verifyEqual(ftInfo.message, "CompareFieldTrip not requested");
        end

        function testCompareFieldTripDegradesGracefullyWhenUnavailable(testCase)
            testCase.assumeTrue(isempty(which("ft_defaults")), ...
                "FieldTrip is on the path on this machine; unavailable-path test not applicable.");
            rng(6);
            x = noiseXu(1500, 2.0, 3.0, 100);
            testCase.verifyWarning(@() estimateIRASA(x, 100, 1.0, 25.0, CompareFieldTrip=true), ...
                "estimateIRASA:NoFieldTrip");
            [~, ~, ~, ~, ftInfo] = estimateIRASA(x, 100, 1.0, 25.0, CompareFieldTrip=true); %#ok<ASGLU>
        end

        function testCompareFieldTripSuccessPathWhenAvailable(testCase)
            testCase.assumeTrue(~isempty(which("ft_defaults")), ...
                "FieldTrip not on path on this machine; success-path test not applicable.");
            rng(9);
            x = noiseXu(3000, 2.5, 4.0, 120);
            [alphaHB, ~, ~, ~, ftInfo] = testCase.verifyWarningFree(...
                @() estimateIRASA(x, 120, 1.0, 30.0, CompareFieldTrip=true));
            testCase.verifyTrue(ftInfo.available);
            testCase.verifyTrue(isfinite(ftInfo.alphaFT));
            testCase.verifyTrue(isfinite(ftInfo.gapHBFT));
            testCase.verifyEqual(ftInfo.gapHBFT, abs(alphaHB - ftInfo.alphaFT), "AbsTol", 1e-9);
            testCase.verifyEqual(ftInfo.message, "");
            % Sanity, not precision: two independent-ish IRASA implementations
            % on the same signal should not be wildly discrepant.
            testCase.verifyLessThan(ftInfo.gapHBFT, 1.0);
        end

        function testComparePMTMSuccessPath(testCase)
            % No external dependency (unlike CompareFieldTrip): should
            % always succeed when requested, regardless of machine.
            rng(9);
            x = noiseXu(3000, 2.5, 4.0, 120);
            [alphaHB, ~, ~, ~, ~, pmInfo] = testCase.verifyWarningFree(...
                @() estimateIRASA(x, 120, 1.0, 30.0, ComparePMTM=true));
            testCase.verifyTrue(pmInfo.available);
            testCase.verifyTrue(isfinite(pmInfo.alphaPM));
            testCase.verifyTrue(isfinite(pmInfo.gapHBPM));
            testCase.verifyEqual(pmInfo.gapHBPM, abs(alphaHB - pmInfo.alphaPM), "AbsTol", 1e-9);
            testCase.verifyEqual(pmInfo.message, "");
            % Sanity, not precision: HB and PM share the same PSD estimator
            % and differ only in the resampling-median step, so at this
            % alpha (comfortably inside the [-2,6] verified range) they
            % should not be wildly discrepant -- see
            % README_IRASAThreeWayAgreement_v001.md Section 3, alpha=2-3 rows.
            testCase.verifyLessThan(pmInfo.gapHBPM, 1.0);
        end

        function testDefaultPmInfoWhenNotRequested(testCase)
            rng(4);
            x = noiseXu(1500, 2.0, 3.0, 100);
            [~, ~, ~, ~, ~, pmInfo] = estimateIRASA(x, 100, 1.0, 25.0);
            testCase.verifyFalse(pmInfo.available);
            testCase.verifyEqual(pmInfo.message, "ComparePMTM not requested");
        end

        function testPRawMatchesDirectPmtmCall(testCase)
            % Source-equivalence for the pRaw plumbing added alongside
            % ComparePMTM: the h=1 PSD returned internally must be bit-
            % identical to calling pmtm(x,4,[],fs) directly -- this is
            % what ComparePMTM's alphaPM is fit on, so if this drifts,
            % ComparePMTM silently stops being the resampling-free
            % estimator it claims to be.
            rng(11);
            x = randn(2048, 1);
            fs = 100;
            [pDirect, fDirect] = pmtm(x, 4, [], fs);

            [~, ~, ~, fVec, ~, pmInfo] = estimateIRASA(x, fs, 1.0, 25.0, ComparePMTM=true); %#ok<ASGLU>
            testCase.verifyEqual(fVec, fDirect);

            fitMask = fDirect >= 1.0 & fDirect <= 25.0 & pDirect > 0;
            pFit = polyfit(log10(fDirect(fitMask)), log10(pDirect(fitMask)), 1);
            alphaDirect = -pFit(1);
            testCase.verifyEqual(pmInfo.alphaPM, alphaDirect, "AbsTol", 1e-12);
        end

        function testNewOptionDefaultsAreUnchanged(testCase)
            rng(12);  x = noiseXu(2000, 3.0, 1.0, 240);
            [a1, s1, p1, f1] = estimateIRASA(x, 240, 1, 20);
            [a2, s2, p2, f2] = estimateIRASA(x, 240, 1, 20, PSD="pmtm", Prewhiten=0);
            testCase.verifyEqual([a1 s1], [a2 s2]);
            testCase.verifyTrue(isequaln(p1, p2));
            testCase.verifyEqual(f1, f2);
        end

        function testPrewhitenFixesBreakFreeOvershoot(testCase)
            % checkIRASAOptions_v001: break-free alpha 4.5 at fs 240, N 3000 reads
            % about +2.1 high by default and within 0.1 with Prewhiten=2.
            % Measured with the pre-2.0 resampling, so it is pinned to Resampling="legacy".
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("estimateIRASA:legacyResampling"));
            rng(13);
            X = cell2mat(arrayfun(@(~) noiseXu(3000, 4.5, 1, 240, Phi=1), 1:10, 'UniformOutput', false));
            aDef = arrayfun(@(j) estimateIRASA(X(:, j), 240, 1, 20, Resampling="legacy"), 1:10);
            aPw  = arrayfun(@(j) estimateIRASA(X(:, j), 240, 1, 20, Prewhiten=2, Resampling="legacy"), 1:10);
            testCase.verifyGreaterThan(median(aDef), 5.5, "Default overshoot on break-free noise has changed; update the RANGE help text.");
            testCase.verifyLessThan(abs(median(aPw) - 4.5), 0.25);
        end

        function testWelchRunsAndReturnsWelchGrid(testCase)
            rng(14);  x = noiseXu(3000, 2.0, 1.0, 240);
            [a, s, p, f] = estimateIRASA(x, 240, 1, 20, PSD="welch");
            testCase.verifyTrue(isfinite(a) && isfinite(s));
            testCase.verifyEqual(numel(p), numel(f));
            testCase.verifyEqual(numel(f), 1024 / 2 + 1);   % nfft = 1024 for a 1024-sample segment
        end

        function testPrewhitenWithFieldTripErrors(testCase)
            testCase.verifyError(@() estimateIRASA(randn(500, 1), 100, 1, 20, Prewhiten=1, CompareFieldTrip=true), ...
                "estimateIRASA:prewhitenFieldTrip");
        end

        function testResamplingDefaultIsIrasa(testCase)
            rng(15);  x = noiseXu(2000, 3.0, 1.0, 240);
            [a1, s1, p1, f1] = testCase.verifyWarningFree(@() estimateIRASA(x, 240, 1, 20));
            [a2, s2, p2, f2] = estimateIRASA(x, 240, 1, 20, Resampling="irasa");
            testCase.verifyEqual([a1 s1], [a2 s2]);
            testCase.verifyTrue(isequaln(p1, p2));
            testCase.verifyEqual(f1, f2);
        end

        function testLegacyResamplingWarns(testCase)
            rng(17);  x = noiseXu(2000, 2.0, 1.0, 240);
            testCase.verifyWarning(@() estimateIRASA(x, 240, 1, 20, Resampling="legacy"), "estimateIRASA:legacyResampling");
        end

        function testIrasaResamplingRemovesPeak(testCase)
            % A 10 Hz sine on alpha-1.5 noise: the default removes the peak, "legacy" keeps it (documented).
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("estimateIRASA:legacyResampling"));
            fs = 240;  N = 3000;  t = (0:N-1)' / fs;
            for psd = ["pmtm", "welch"]
                dI = zeros(1, 10);  dT = zeros(1, 10);
                for j = 1:10
                    rng(100 + j);  x = noiseXu(N, 1.5, 1, fs) + 0.5 * sin(2*pi*10*t + 2*pi*rand);
                    [~, ~, pI, f] = estimateIRASA(x, fs, 1, 20, PSD=psd);
                    [~, ~, pT] = estimateIRASA(x, fs, 1, 20, PSD=psd, Resampling="legacy");
                    dI(j) = testEstimateIRASA.peakDb(f, pI);  dT(j) = testEstimateIRASA.peakDb(f, pT);
                end
                testCase.verifyLessThan(median(dI), 4, "The default left the 10 Hz peak in pFractal (" + psd + ").");
                testCase.verifyGreaterThan(median(dT), 10, "Resampling=""legacy"" now removes peaks (" + psd + "); update the help text.");
            end
        end

        function testIrasaResamplingAgreesOnPeakFreeNoise(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("estimateIRASA:legacyResampling"));
            fs = 240;  aI = zeros(1, 10);  aT = zeros(1, 10);
            for j = 1:10
                rng(200 + j);  x = noiseXu(3000, 2.0, 1, fs);
                aI(j) = estimateIRASA(x, fs, 1, 20);
                aT(j) = estimateIRASA(x, fs, 1, 20, Resampling="legacy");
            end
            testCase.verifyLessThan(abs(median(aI) - 2.0), 0.2);
            testCase.verifyLessThan(abs(median(aI) - median(aT)), 0.15);
        end

        function testIrasaResamplingWithPrewhitenOnSteepBreakFree(testCase)
            % Per-h recolouring: break-free alpha 4.5 must come back with Prewhiten=2 under "irasa", both PSDs.
            for psd = ["pmtm", "welch"]
                a = zeros(1, 10);
                for j = 1:10
                    rng(300 + j);  x = noiseXu(3000, 4.5, 1, 240, Phi=1);
                    a(j) = estimateIRASA(x, 240, 1, 20, PSD=psd, Prewhiten=2, Resampling="irasa");
                end
                testCase.verifyLessThan(abs(median(a) - 4.5), 0.3, "PSD = " + psd);
            end
        end

        function testIrasaResamplingMasksAboveValidRange(testCase)
            rng(16);  x = noiseXu(2000, 2.0, 1, 240);
            [~, s, p, f] = estimateIRASA(x, 240, 1, 20, Resampling="irasa");
            testCase.verifyTrue(all(isnan(p(f * 1.9 >= 120))));
            testCase.verifyTrue(all(isfinite(p(f >= 1 & f <= 20))));
            testCase.verifyTrue(isfinite(s));
        end

        function testIrasaResamplingBandTooHighErrors(testCase)
            testCase.verifyError(@() estimateIRASA(randn(2000, 1), 100, 1, 30, Resampling="irasa"), "estimateIRASA:irasaBand");
        end
    end

    methods (Static)
        function d = peakDb(f, p)
            % Peak in 9.5-10.5 Hz over the median of 7-8.5 and 12-14 Hz, in dB.
            nb = (f >= 7 & f <= 8.5) | (f >= 12 & f <= 14);
            d = 10 * log10(max(p(f >= 9.5 & f <= 10.5)) / median(p(nb), 'omitnan'));
        end
    end
end
