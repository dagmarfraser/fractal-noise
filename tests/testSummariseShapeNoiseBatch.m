classdef testSummariseShapeNoiseBatch < matlab.unittest.TestCase
%TESTSUMMARISESHAPENOISEBATCH Unit tests for summariseShapeNoiseBatch.m and the shapeNoise no-accepted-draw message.
%
% Toolbox: fractalnoise
% The binomial answers below were computed independently of the function, in exact rational arithmetic (Python fractions).
% Random walks are used as the integration reference because a cumulative sum of white noise has a known exponent (alpha 2).

    methods (TestClassSetup)
        function addToolboxPath(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
        end
    end

    methods (Static)
        function infos = batch(n, nAccepted)
            % n single-channel infos, the first nAccepted of them accepted
            infos = cell(n, 1);
            for k = 1:n
                infos{k} = struct("accepted", k <= nAccepted, "nDraw", 1 + 19 * (k > nAccepted));
            end
        end
    end

    methods (Test)
        function testExactLowerTailProbabilities(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:lowAcceptance"));
            known = [42 39 0.00858385379973348; 42 40 0.066182322212583; 50 46 0.00159617309073287; 50 47 0.0138172708306004; ...
                100 95 0.00343232158775452; 100 96 0.0183740364446497; 200 192 0.00101255715052078];
            for r = 1:size(known, 1)
                S = summariseShapeNoiseBatch(testCase.batch(known(r, 1), known(r, 2)));
                testCase.verifyEqual(S.pValue, known(r, 3), "RelTol", 1e-10);
            end
        end

        function testCriticalCountsAtDefaults(testCase)
            % largest accepted count that warns, from the exact binomial, Standard 0.99, Level 0.01
            critical = [42 39; 50 46; 100 95; 200 193];
            for r = 1:size(critical, 1)
                n = critical(r, 1);  c = critical(r, 2);
                S = testCase.verifyWarning(@() summariseShapeNoiseBatch(testCase.batch(n, c)), "shapeNoise:lowAcceptance");
                testCase.verifyEqual(S.decision, "low");
                S = testCase.verifyWarningFree(@() summariseShapeNoiseBatch(testCase.batch(n, c + 1)));
                testCase.verifyEqual(S.decision, "ok");
            end
        end

        function testDefaultsAreTheDesignStandard(testCase)
            S = summariseShapeNoiseBatch(testCase.batch(50, 50));
            testCase.verifyEqual([S.Standard, S.Level, S.MinChannels], [0.99, 0.01, 42]);
        end

        function testPValueRisesWithAcceptance(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:lowAcceptance"));
            p = arrayfun(@(x) summariseShapeNoiseBatch(testCase.batch(100, x)).pValue, 80:100);
            testCase.verifyGreaterThan(diff(p), 0);
            testCase.verifyEqual(p(end), 1, "AbsTol", 1e-12);
        end

        function testSilentAtOrAboveStandard(testCase)
            for x = [99 100]
                S = testCase.verifyWarningFree(@() summariseShapeNoiseBatch(testCase.batch(100, x)));
                testCase.verifyEqual(S.decision, "ok");
                testCase.verifyEqual(S.acceptance, x / 100);
            end
        end

        function testBelowMinimumReportsWithoutTesting(testCase)
            S = testCase.verifyWarningFree(@() summariseShapeNoiseBatch(testCase.batch(41, 0)));
            testCase.verifyEqual(S.decision, "tooFew");
            testCase.verifyTrue(isnan(S.pValue));
            testCase.verifyEqual(S.acceptance, 0);
            testCase.verifyEqual(S.nChannels, 41);
        end

        function testMinimumIsInclusive(testCase)
            S = testCase.verifyWarning(@() summariseShapeNoiseBatch(testCase.batch(42, 0)), "shapeNoise:lowAcceptance");
            testCase.verifyEqual(S.decision, "low");
        end

        function testMinChannelsOption(testCase)
            testCase.verifyWarning(@() summariseShapeNoiseBatch(testCase.batch(10, 0), MinChannels=10), "shapeNoise:lowAcceptance");
            S = testCase.verifyWarningFree(@() summariseShapeNoiseBatch(testCase.batch(10, 0), MinChannels=11));
            testCase.verifyEqual(S.decision, "tooFew");
        end

        function testLevelOption(testCase)
            testCase.verifyWarningFree(@() summariseShapeNoiseBatch(testCase.batch(100, 96)));                       % p = 0.0184
            testCase.verifyWarning(@() summariseShapeNoiseBatch(testCase.batch(100, 96), Level=0.05), "shapeNoise:lowAcceptance");
        end

        function testStandardOption(testCase)
            % Standard 0.95, n 100: the largest accepted count that warns is 88
            testCase.verifyWarning(@() summariseShapeNoiseBatch(testCase.batch(100, 88), Standard=0.95), "shapeNoise:lowAcceptance");
            testCase.verifyWarningFree(@() summariseShapeNoiseBatch(testCase.batch(100, 89), Standard=0.95));
        end

        function testTriesKnownAnswers(testCase)
            infos = arrayfun(@(d) struct("accepted", true, "nDraw", d), (1:10)', "UniformOutput", false);
            S = summariseShapeNoiseBatch(infos);
            testCase.verifyEqual(S.triesMedian, 5.5);
            testCase.verifyEqual(S.triesP90, 9);
            infos = arrayfun(@(d) struct("accepted", true, "nDraw", d), (100:-1:1)', "UniformOutput", false);
            S = summariseShapeNoiseBatch(infos);
            testCase.verifyEqual(S.triesP90, 90);
        end

        function testPerDrawAcceptanceKnownAnswer(testCase)
            % accepted after 1, 2 and 3 draws, and one channel that used all 20: 3 acceptances in 26 draws
            infos = {struct("accepted", true, "nDraw", 1), struct("accepted", true, "nDraw", 2), struct("accepted", true, "nDraw", 3), struct("accepted", false, "nDraw", 20)};
            S = summariseShapeNoiseBatch(infos);
            testCase.verifyEqual(S.perDrawAcceptance, 3 / 26, "AbsTol", 1e-12);
        end

        function testCountsChannelsNotCalls(testCase)
            infos = repmat({struct("accepted", [true false], "nDraw", [1 20])}, 25, 1);
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:lowAcceptance"));
            S = summariseShapeNoiseBatch(infos);
            testCase.verifyEqual([S.nInfos, S.nChannels, S.nAccepted], [25, 50, 25]);
            testCase.verifyEqual(S.acceptance, 0.5);
        end

        function testStructArrayEqualsCellArray(testCase)
            infos = testCase.batch(60, 59);
            testCase.verifyEqual(summariseShapeNoiseBatch([infos{:}]), summariseShapeNoiseBatch(infos));
        end

        function testBadInputsError(testCase)
            good = struct("accepted", true, "nDraw", 1);
            testCase.verifyError(@() summariseShapeNoiseBatch({}), "summariseShapeNoiseBatch:empty");
            testCase.verifyError(@() summariseShapeNoiseBatch(5), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({good, 5}), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({struct("accepted", true)}), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({struct("nDraw", 1)}), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({struct("accepted", 1, "nDraw", 1)}), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({struct("accepted", [true true], "nDraw", 1)}), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({struct("accepted", true, "nDraw", 0)}), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({struct("accepted", true, "nDraw", 1.5)}), "summariseShapeNoiseBatch:badInfo");
            testCase.verifyError(@() summariseShapeNoiseBatch({struct("accepted", true, "nDraw", NaN)}), "summariseShapeNoiseBatch:badInfo");
        end

        function testBadOptionsError(testCase)
            infos = testCase.batch(50, 50);
            testCase.verifyError(@() summariseShapeNoiseBatch(infos, Standard=1), ?MException);
            testCase.verifyError(@() summariseShapeNoiseBatch(infos, Standard=0), ?MException);
            testCase.verifyError(@() summariseShapeNoiseBatch(infos, Level=0), ?MException);
            testCase.verifyError(@() summariseShapeNoiseBatch(infos, Level=1), ?MException);
            testCase.verifyError(@() summariseShapeNoiseBatch(infos, MinChannels=0), ?MException);
            testCase.verifyError(@() summariseShapeNoiseBatch(infos, MinChannels=1.5), ?MException);
        end

        function testLegacyInfoIsRefused(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:legacyConstruction"));
            rng(3);  [~, info] = shapeNoise(cumsum(randn(400, 1)), 133, 2, Method="legacy");
            testCase.verifyError(@() summariseShapeNoiseBatch({info}), "summariseShapeNoiseBatch:badInfo");
        end

        function testWarningSaysWhatItMeansAndWhatToDo(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:lowAcceptance"));
            lastwarn("");  summariseShapeNoiseBatch(testCase.batch(50, 40));
            [msg, id] = lastwarn;
            testCase.verifyEqual(string(id), "shapeNoise:lowAcceptance");
            for phrase = ["40 of 50", "not a clean power law", "AlphaTol", "MaxDraw", "Inspect the residual"]
                testCase.verifySubstring(string(msg), phrase);
            end
        end

        function testVeryLowAcceptanceWarnsAndNeverErrors(testCase)
            S = testCase.verifyWarning(@() summariseShapeNoiseBatch(testCase.batch(500, 0)), "shapeNoise:lowAcceptance");
            testCase.verifyEqual(S.acceptance, 0);
            testCase.verifyEqual(S.decision, "low");
        end

        function testSilentOnRandomWalks(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture(["shapeNoise:noAcceptedDraw", "shapeNoise:estimateFailed"]));
            rng(20261118, "twister");  infos = cell(60, 1);
            for k = 1:60
                [~, infos{k}] = shapeNoise(cumsum(randn(1000, 1)), 133, [], AlphaBand=[1 20]);
            end
            S = testCase.verifyWarningFree(@() summariseShapeNoiseBatch(infos));
            testCase.verifyEqual(S.decision, "ok");
            testCase.verifyGreaterThanOrEqual(S.acceptance, 0.99);
        end

        function testWarnsOnBatchWithUnreachableTargets(testCase)
            % A fifth of the channels are asked for an exponent their amplitude spectrum cannot give
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture(["shapeNoise:noAcceptedDraw", "shapeNoise:estimateFailed"]));
            rng(20261119, "twister");  infos = cell(60, 1);
            for k = 1:60
                x = cumsum(randn(1000, 1));
                if mod(k, 5) == 0, [~, infos{k}] = shapeNoise(x, 133, 6.5, AlphaBand=[1 20], MaxDraw=5);
                else, [~, infos{k}] = shapeNoise(x, 133, [], AlphaBand=[1 20]);
                end
            end
            S = testCase.verifyWarning(@() summariseShapeNoiseBatch(infos), "shapeNoise:lowAcceptance");
            testCase.verifyEqual(S.decision, "low");
            testCase.verifyLessThanOrEqual(S.nAccepted, 48);   % the 12 unreachable channels can never be accepted; an ordinary channel may also miss
        end

        function testNoAcceptedDrawMessageNamesTheLimits(testCase)
            rng(7, "twister");  x = cumsum(randn(1000, 1));
            lastwarn("");  warnState = warning("off", "shapeNoise:noAcceptedDraw");  restore = onCleanup(@() warning(warnState));
            [~, info] = shapeNoise(x, 133, 6.5, AlphaBand=[1 20], MaxDraw=3);
            [msg, id] = lastwarn;
            testCase.verifyFalse(info.accepted);
            testCase.verifyEqual(string(id), "shapeNoise:noAcceptedDraw");
            testCase.verifySubstring(string(msg), "AlphaTol = 0.30");
            testCase.verifySubstring(string(msg), "MaxDraw = 3");
        end
    end
end
