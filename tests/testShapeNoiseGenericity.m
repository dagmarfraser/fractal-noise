classdef testShapeNoiseGenericity < matlab.unittest.TestCase
%TESTSHAPENOISEGENERICITY Tests for the v2.0.0 generality options of shapeNoise: AlphaBand handling and CheckUncovered.
%
% Toolbox: fractalnoise
    methods (TestClassSetup)
        function addToolbox(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
        end
    end

    methods (Test)
        function testDefaultBandWarnsAndKinematicIsSilent(testCase)
            rng(1);  r = noiseXu(2000, 3, 1, 100);  r = r - mean(r);
            testCase.verifyWarning(@() shapeNoise(r, 100, []), "shapeNoise:defaultAlphaBand");
            testCase.verifyWarningFree(@() shapeNoise(r, 100, [], AlphaBand="kinematic"));
        end

        function testKinematicKeywordEqualsTheDefaultBand(testCase)
            rng(2);  r = noiseXu(2000, 3, 1, 100);  r = r - mean(r);
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:defaultAlphaBand"));
            rng(10);  [a, ia] = shapeNoise(r, 100, []);
            rng(10);  [b, ib] = shapeNoise(r, 100, [], AlphaBand="kinematic");
            rng(10);  [c, ~] = shapeNoise(r, 100, [], AlphaBand=[1 20]);
            testCase.verifyEqual(a, b);  testCase.verifyEqual(b, c);  testCase.verifyEqual(ia.alphaTarget, ib.alphaTarget);
        end

        function testExplicitBandInAnotherUnit(testCase)
            rng(3);  fs = 12;  r = noiseXu(1200, 1.5, 1, fs);  r = r - mean(r);   % monthly data, fs = 12 per year
            [n, info] = shapeNoise(r, fs, [], AlphaBand=[0.2 5]);
            testCase.verifyEqual(size(n, 1), info.L);
            testCase.verifyTrue(all(isfinite(n)));
        end

        function testInvalidBandsAreErrors(testCase)
            rng(4);  r = noiseXu(1500, 3, 1, 100);  r = r - mean(r);
            testCase.verifyError(@() shapeNoise(r, 100, [], AlphaBand=[5 1]), "shapeNoise:alphaBand");
            testCase.verifyError(@() shapeNoise(r, 100, [], AlphaBand=[1 60]), "shapeNoise:alphaBand");   % above Nyquist (50)
            testCase.verifyError(@() shapeNoise(r, 100, [], AlphaBand="fast"), "shapeNoise:alphaBand");
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:defaultAlphaBand"));
            testCase.verifyError(@() shapeNoise(r, 2, []), "shapeNoise:alphaBand");   % the kinematic default is empty at fs = 2
        end

        function testCheckUncoveredOptionControlsTheZeroTailError(testCase)
            rng(5);  r = noiseXu(1500, 2, 1, 100);  r = r - mean(r);  r(1:10) = 0;   % leading zeros
            testCase.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture("shapeNoise:defaultAlphaBand"));
            testCase.verifyError(@() shapeNoise(r, 100, []), "shapeNoise:zeroTail");
            testCase.verifyError(@() shapeNoise(r, 100, [], CheckUncovered=true), "shapeNoise:zeroTail");
            n = shapeNoise(r, 100, [], CheckUncovered=false);
            testCase.verifyTrue(all(isfinite(n(:))) || any(isnan(n(:))));   % runs; a NaN channel is the documented no-draw outcome
        end
    end
end
