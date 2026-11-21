classdef testWrappers < matlab.unittest.TestCase
    % Tests for shapeNoiseWrapper and cetNoiseWrapper.
    methods (TestClassSetup)
        function addPath(testCase)
            here = fileparts(mfilename("fullpath"));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(here, "..", "toolbox")));
        end
    end
    methods (Test)
        function testNoTemplateEqualsShapeNoise(testCase)
            rng(3);  r = cumsum(cumsum(randn(1500, 1))) * 0.01;  r = r - mean(r);
            rng(11);  a = shapeNoise(r, 100, [], AlphaBand=[1 20]);
            rng(11);  b = shapeNoiseWrapper(r, 100, AlphaBand=[1 20]);
            testCase.verifyEqual(b, a);
        end
        function testTemplateStepRunsAndRemovesPeriodicPart(testCase)
            rng(4);  t = (0:1/100:20)';  e = cumsum(cumsum(randn(numel(t), 1))) * 0.002;
            x = 2 * sin(2 * pi * 0.8 * t) + e;
            [n, info] = shapeNoiseWrapper(x, 100, F0=0.8, K=1, AlphaBand=[1 20]);
            testCase.verifyNotEmpty(info.template);
            testCase.verifyEqual(info.template.f0, 0.8);
            testCase.verifyLessThan(std(info.residual), 0.5 * std(x));
            testCase.verifyEqual(size(n, 1), info.shape.L);
            testCase.verifyTrue(all(isfinite(n), "all"));
        end
        function testF0BandIsEstimated(testCase)
            rng(5);  t = (0:1/100:20)';  x = 2 * sin(2 * pi * 0.8 * t) + cumsum(cumsum(randn(numel(t), 1))) * 0.0003;
            [~, info] = shapeNoiseWrapper(x, 100, F0Band=[0.5 1.5], K=1, AlphaBand=[1 20]);   % a narrow band: on steep noise the highest peak can sit at a wide band's low edge
            testCase.verifyEqual(info.template.f0, 0.8, AbsTol=0.05);
        end
        function testMonthlySeriesLoads(testCase)
            c = cetNoiseWrapper(LoadOnly=true);
            testCase.verifyEqual(c.fs, 12);
            testCase.verifyEqual(c.year(1), 1659);
            testCase.verifyEqual(numel(c.value), numel(c.year));
            testCase.verifyEqual(mod(numel(c.value), 12), 0);
            testCase.verifyTrue(all(isfinite(c.value)) && all(c.value > -20 & c.value < 30));
            testCase.verifyEqual(c.month(1:12)', 1:12);
        end
        function testLastYearIsTheLastCompleteYear(testCase)
            c = cetNoiseWrapper(LoadOnly=true);
            testCase.verifyEqual(c.year(end), 2025);          % 2026 has no November or December in the bundled file
            testCase.verifyError(@() cetNoiseWrapper(LoadOnly=true, LastYear=2026), "cetNoiseWrapper:lastYear");
        end
        function testYearRangeRestriction(testCase)
            c = cetNoiseWrapper(LoadOnly=true, FirstYear=1900, LastYear=1999);
            testCase.verifyEqual(numel(c.value), 1200);
            testCase.verifyError(@() cetNoiseWrapper(LoadOnly=true, FirstYear=2200, LastYear=2300), "cetNoiseWrapper:lastYear");
        end
        function testKnownMonthValue(testCase)
            c = cetNoiseWrapper(LoadOnly=true, FirstYear=2025, LastYear=2025);   % Jan 2025 3.4, Jul 18.4 (file as downloaded 2026-11-17)
            testCase.verifyEqual(c.value([1 7]), [3.4; 18.4]);
        end
        function testCetSurrogateRunsInMonthlyUnits(testCase)
            rng(7);  [n, info] = cetNoiseWrapper(FirstYear=1900, LastYear=2025, K=2);
            testCase.verifyEqual(info.cet.fs, 12);
            testCase.verifyEqual(info.template.f0, 1);
            testCase.verifyTrue(all(isfinite(n), "all"));
            testCase.verifyEqual(size(n, 1), info.shape.L);
        end
        function testBadBandIsAnError(testCase)
            testCase.verifyError(@() cetNoiseWrapper(FirstYear=1900, LastYear=2000, K=2, AlphaBand=[0.5 7]), "shapeNoise:alphaBand");
        end
        function testMissingMonthInsideTheChosenYearsIsAnError(testCase)
            f = string(tempname) + ".txt";  cleanup = onCleanup(@() delete(f));  %#ok<NASGU>
            fid = fopen(f, "w");
            fprintf(fid, "header\n  Year    Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep    Oct    Nov    Dec  Annual\n");
            fprintf(fid, "  1999    1.0    2.0    3.0    4.0    5.0    6.0    7.0    8.0    9.0   10.0   11.0   12.0   6.50\n");
            fprintf(fid, "  2000    1.0    2.0    3.0  -99.9    5.0    6.0    7.0    8.0    9.0   10.0   11.0   12.0  -99.90\n");
            fprintf(fid, "  2001    1.0    2.0    3.0    4.0    5.0    6.0    7.0    8.0    9.0   10.0   11.0   12.0   6.50\n");
            fclose(fid);
            testCase.verifyError(@() cetNoiseWrapper(File=f, FirstYear=1999, LastYear=2001, LoadOnly=true), "cetNoiseWrapper:missing");
            c = cetNoiseWrapper(File=f, FirstYear=2001, LastYear=2001, LoadOnly=true);
            testCase.verifyEqual(c.value', 1:12);
        end
        function testKRequiredWhenATemplateIsStripped(testCase)
            % K has no default: a template step without K is an error, not a silent K = 1 and not a warning.
            rng(4);  t = (0:1/100:20)';  x = 2 * sin(2 * pi * 0.8 * t) + cumsum(cumsum(randn(numel(t), 1))) * 0.002;
            testCase.verifyError(@() shapeNoiseWrapper(x, 100, F0=0.8, AlphaBand=[1 20]), "shapeNoiseWrapper:KRequired");
            testCase.verifyError(@() shapeNoiseWrapper(x, 100, F0Band=[0.5 1.5], AlphaBand=[1 20]), "shapeNoiseWrapper:KRequired");
            testCase.verifyError(@() shapeNoiseWrapper(x, 100, F0=0.8, K=[], AlphaBand=[1 20]), "shapeNoiseWrapper:KRequired");
        end
        function testKRequiredMessagePointsToHarmonicContent(testCase)
            x = sin(2 * pi * (0:1/100:20)');
            try
                shapeNoiseWrapper(x, 100, F0=1, AlphaBand=[1 20]);
                testCase.verifyFail("No error without K.");
            catch err
                testCase.verifyEqual(string(err.identifier), "shapeNoiseWrapper:KRequired");
                testCase.verifySubstring(err.message, "harmonicContent");
            end
            try
                cetNoiseWrapper(FirstYear=1900, LastYear=1950);
                testCase.verifyFail("No error without K.");
            catch err
                testCase.verifyEqual(string(err.identifier), "cetNoiseWrapper:KRequired");
                testCase.verifySubstring(err.message, "harmonicContent");
            end
        end
        function testExplicitKIsExactlyTheTwoStepPath(testCase)
            % Known answer: the wrapper with K given is harmonicTemplate with that K, then shapeNoise, bit for bit.
            rng(4);  t = (0:1/100:20)';  x = 2 * sin(2 * pi * 0.8 * t) + 0.6 * sin(2 * pi * 1.6 * t) + cumsum(cumsum(randn(numel(t), 1))) * 0.002;
            res = harmonicTemplate(x, 100, 0.8, 2, Window=4);
            rng(21);  a = shapeNoise(res, 100, [], AlphaBand=[1 20]);
            rng(21);  [b, info] = shapeNoiseWrapper(x, 100, F0=0.8, K=2, AlphaBand=[1 20]);
            testCase.verifyEqual(info.template.K, 2);
            testCase.verifyEqual(info.residual, res);
            testCase.verifyEqual(b, a);
        end
        function testCetKRequired(testCase)
            testCase.verifyError(@() cetNoiseWrapper(), "cetNoiseWrapper:KRequired");
            testCase.verifyError(@() cetNoiseWrapper(FirstYear=1900, LastYear=2025), "cetNoiseWrapper:KRequired");
            c = cetNoiseWrapper(LoadOnly=true, FirstYear=1900, LastYear=1901);    % loading alone needs no K
            testCase.verifyEqual(numel(c.value), 24);
        end
        function testCetKChosenWithHarmonicContent(testCase)
            % Known answer for the README workflow (bundled file, 1900-2025): ratios 73.83, 6.00, 0.56 for harmonics 1-3, so K = 2.
            c = cetNoiseWrapper(LoadOnly=true, FirstYear=1900, LastYear=2025);
            [K, hc] = harmonicContent(c.value, 12, 1);
            testCase.verifyEqual(K, 2);
            testCase.verifyEqual(hc.firstBelow, 3);
            testCase.verifyFalse(hc.noHarmonic);
            testCase.verifyEqual(hc.ratio(1:3), [73.83 6.00 0.56], AbsTol=0.01);
            rng(7);  [a, info] = cetNoiseWrapper(FirstYear=1900, LastYear=2025, K=K);
            rng(7);  b = shapeNoiseWrapper(c.value, 12, F0=1, K=2, Window=4, AlphaBand=[0.5 5]);
            testCase.verifyEqual(info.template.K, 2);
            testCase.verifyEqual(a, b);
        end
    end
end
