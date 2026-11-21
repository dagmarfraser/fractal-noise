classdef testHarmonicContent < matlab.unittest.TestCase
    % Known-answer, property and fail-loud tests for harmonicContent.
    methods (TestClassSetup)
        function addPath(testCase)
            here = fileparts(mfilename("fullpath"));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, "..", "toolbox")));
        end
    end
    methods (Static)
        function x = drawing(nh, seed)
            rng(seed);  fs = 100;  t = (0:1/fs:30)';  amps = [2 0.8 0.4 0.2 0.1];  x = zeros(size(t));
            for h = 1:nh, x = x + amps(h) * sin(2 * pi * h * 1.2 * t + 2 * pi * rand); end
            x = x + cumsum(randn(size(t))) * 0.002;
        end
    end
    methods (Test)
        function testCountsTheTrueHarmonics(testCase)
            for nh = [1 2 3 4]
                K = harmonicContent(testHarmonicContent.drawing(nh, 10 + nh), 100, 1.2);
                testCase.verifyEqual(K, nh, "true harmonics: " + nh);
            end
        end
        function testRatioIsLargeForHarmonicsAndNearOneBeyond(testCase)
            [~, info] = harmonicContent(testHarmonicContent.drawing(3, 21), 100, 1.2);
            testCase.verifyGreaterThan(min(info.ratio(1:3)), 10);
            testCase.verifyLessThan(median(info.ratio(5:end)), 3);
            testCase.verifyEqual(info.firstBelow, 4);
        end
        function testNoHarmonicIsFlaggedAndFloorsAtOne(testCase)
            rng(5);  x = cumsum(randn(3001, 1)) * 0.01;
            [K, info] = harmonicContent(x, 100, 1.2);
            testCase.verifyEqual(K, 1);
            testCase.verifyTrue(info.noHarmonic);
        end
        function testCutChangesTheCount(testCase)
            x = testHarmonicContent.drawing(4, 31);
            [~, info] = harmonicContent(x, 100, 1.2);
            hi = 0.5 * (info.ratio(4) + info.ratio(5));
            testCase.verifyEqual(harmonicContent(x, 100, 1.2, Cut=hi), 4);
            testCase.verifyGreaterThan(harmonicContent(x, 100, 1.2, Cut=1.0), 4);
        end
        function testIsScaleInvariant(testCase)
            x = testHarmonicContent.drawing(3, 41);
            [~, a] = harmonicContent(x, 100, 1.2);  [~, b] = harmonicContent(7 * x, 100, 1.2);
            testCase.verifyEqual(a.ratio, b.ratio, RelTol=1e-8);
        end
        function testTwoChannelsAverageAmplitudes(testCase)
            x1 = testHarmonicContent.drawing(2, 51);  x2 = testHarmonicContent.drawing(2, 52);
            [K, info] = harmonicContent([x1 x2], 100, 1.2);
            testCase.verifyEqual(K, 2);
            testCase.verifyEqual(numel(info.ratio), info.Kmax);
        end
        function testF0IsEstimatedWithABand(testCase)
            x = testHarmonicContent.drawing(3, 61);
            [K, info] = harmonicContent(x, 100, [], F0Band=[0.8 1.6]);
            testCase.verifyEqual(info.f0, 1.2, AbsTol=0.03);
            testCase.verifyEqual(K, 3);
        end
        function testMonthlyUnits(testCase)
            rng(3);  t = (0:1/12:60)';  x = 5 * sin(2 * pi * t) + 1.5 * sin(4 * pi * t + 1) + randn(size(t)) * 0.3;
            testCase.verifyEqual(harmonicContent(x, 12, 1), 2);
        end
        function testKmaxIsLimitedByNyquist(testCase)
            [~, info] = harmonicContent(testHarmonicContent.drawing(2, 71), 100, 10);
            testCase.verifyEqual(info.Kmax, 4);   % (k + 0.5) * 10 < 50
        end
        function testErrorPaths(testCase)
            x = testHarmonicContent.drawing(2, 81);
            testCase.verifyError(@() harmonicContent(x, 100, []), "harmonicContent:f0Band");
            testCase.verifyError(@() harmonicContent([x; NaN], 100, 1.2), "harmonicContent:nonFinite");
            testCase.verifyError(@() harmonicContent(x(1:150), 100, 1.2), "harmonicContent:tooShort");
            testCase.verifyError(@() harmonicContent(x, 100, 49), "harmonicContent:nyquist");
        end
        function testChannelsAreAveragedNotFirstOnly(testCase)
            x1 = cumsum(randn(3001, 1)) * 0.01;                    % channel 1: no harmonic at all
            rng(7);  x2 = testHarmonicContent.drawing(3, 91);      % channel 2: three harmonics
            testCase.verifyEqual(harmonicContent(x1, 100, 1.2), 1);
            testCase.verifyEqual(harmonicContent([x1 x2], 100, 1.2), 3);
        end
    end
end
