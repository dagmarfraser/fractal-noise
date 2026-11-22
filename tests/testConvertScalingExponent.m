classdef testConvertScalingExponent < matlab.unittest.TestCase
%TESTCONVERTSCALINGEXPONENT Known answers, round trips and fail-loud errors for convertScalingExponent.
%
% Toolbox: fractalnoise

    methods (TestClassSetup)
        function addToolbox(testCase)
            oldPath = path;  testCase.addTeardown(@() path(oldPath));
            addpath(fullfile(fileparts(mfilename("fullpath")), "..", "toolbox"));
        end
    end

    methods (Test)
        function testKnownAnswersFGn(testCase)
            c = @(x, f, t) convertScalingExponent(x, f, t, Class="fGn");
            testCase.verifyEqual(c(0.72, "H", "alpha"), 0.44, "AbsTol", 1e-12);
            testCase.verifyEqual(c(0.5, "H", "alpha"), 0, "AbsTol", 1e-12, "H 0.5 is white noise.");
            testCase.verifyEqual(c(0.39, "alpha", "H"), 0.695, "AbsTol", 1e-12);
            testCase.verifyEqual(c(1.3, "D", "H"), 0.7, "AbsTol", 1e-12);
            testCase.verifyEqual(c(1.3, "D", "alpha"), 0.4, "AbsTol", 1e-12);
            testCase.verifyEqual(c(0.4, "alpha", "D"), 1.3, "AbsTol", 1e-12);
            testCase.verifyEqual(c(0.65, "H", "D"), 1.35, "AbsTol", 1e-12);
        end

        function testKnownAnswersFBm(testCase)
            c = @(x, f, t) convertScalingExponent(x, f, t, Class="fBm");
            testCase.verifyEqual(c(0.5, "H", "alpha"), 2, "AbsTol", 1e-12, "Brownian motion: H 0.5, alpha 2.");
            testCase.verifyEqual(c(2, "alpha", "D"), 1.5, "AbsTol", 1e-12, "Brownian graph: D 1.5.");
            testCase.verifyEqual(c(1.2, "D", "alpha"), 2.6, "AbsTol", 1e-12);
            testCase.verifyEqual(c(1.2, "D", "H"), 0.8, "AbsTol", 1e-12);
        end

        function testClassFreeConversions(testCase)
            testCase.verifyEqual(convertScalingExponent(0.9, "DFA", "alpha"), 0.8, "AbsTol", 1e-12);
            testCase.verifyEqual(convertScalingExponent(4.5, "alpha", "DFA"), 2.75, "AbsTol", 1e-12, "Steep alpha has a DFA exponent.");
            testCase.verifyEqual(convertScalingExponent(0.4, "d", "alpha"), 0.8, "AbsTol", 1e-12);
            testCase.verifyEqual(convertScalingExponent(-0.5, "alpha", "d"), -0.25, "AbsTol", 1e-12);
            testCase.verifyEqual(convertScalingExponent(3.7, "alpha", "alpha"), 3.7);
        end

        function testRoundTripsAndShape(testCase)
            names = ["alpha", "d", "DFA", "H", "D"];
            for cls = ["fGn", "fBm"]
                a = (-0.95:0.05:0.95)' + 2 * (cls == "fBm");
                a = reshape(a, [], 3);  % keeps the array shape
                for f = names
                    v = convertScalingExponent(a, "alpha", f, Class=cls);
                    testCase.verifySize(v, size(a));
                    for t = names
                        back = convertScalingExponent(convertScalingExponent(v, f, t, Class=cls), t, "alpha", Class=cls);
                        testCase.verifyEqual(back, a, "AbsTol", 1e-12, cls + ": " + f + " -> " + t + " -> alpha");
                    end
                end
            end
        end

        function testHAndDNeedAClass(testCase)
            testCase.verifyError(@() convertScalingExponent(0.7, "H", "alpha"), "convertScalingExponent:classRequired");
            testCase.verifyError(@() convertScalingExponent(0.4, "alpha", "D"), "convertScalingExponent:classRequired");
        end

        function testOutOfRangeIsAnError(testCase)
            e = "convertScalingExponent:outOfRange";
            testCase.verifyError(@() convertScalingExponent(1.5, "alpha", "H", Class="fGn"), e, "alpha 1.5 is not fGn.");
            testCase.verifyError(@() convertScalingExponent(0.5, "alpha", "H", Class="fBm"), e, "alpha 0.5 is not fBm.");
            testCase.verifyError(@() convertScalingExponent(1, "alpha", "H", Class="fGn"), e, "alpha 1 is a boundary.");
            testCase.verifyError(@() convertScalingExponent(-1, "alpha", "D", Class="fGn"), e, "alpha -1 is a boundary.");
            testCase.verifyError(@() convertScalingExponent(3, "alpha", "H", Class="fBm"), e, "alpha 3 is a boundary.");
            testCase.verifyError(@() convertScalingExponent(1.2, "H", "alpha", Class="fGn"), e, "H above 1.");
            testCase.verifyError(@() convertScalingExponent(0, "H", "alpha", Class="fGn"), e, "H 0.");
            testCase.verifyError(@() convertScalingExponent(2.1, "D", "alpha", Class="fBm"), e, "D above 2.");
            testCase.verifyError(@() convertScalingExponent(1, "D", "alpha", Class="fBm"), e, "D 1.");
            testCase.verifyError(@() convertScalingExponent([0.3 1.4], "alpha", "H", Class="fGn"), e, "One bad element is enough.");
        end

        function testBadInputsRejected(testCase)
            testCase.verifyError(@() convertScalingExponent(NaN, "alpha", "d"), "MATLAB:validators:mustBeFinite");
            testCase.verifyError(@() convertScalingExponent(0.3, "beta", "d"), "MATLAB:validators:mustBeMember");
            testCase.verifyError(@() convertScalingExponent(0.3, "alpha", "H", Class="fgn"), "MATLAB:validators:mustBeMember");
        end

        function testAgreesWithTheEstimatorsOnNoiseXu(testCase)
            % noiseXu builds alpha from d = alpha/2; the Whittle alpha and the DFA pre-estimate of one long series convert to the
            % same H within sampling error.
            rng(7);
            x = noiseXu(2 ^ 15, 0.6, 1, 100, Phi=1);
            [a, ~, info] = estimateWhittleAlpha(x, 100, 0.05, 5);
            hW = convertScalingExponent(a, "alpha", "H", Class="fGn");
            hD = convertScalingExponent((info.alphaDFA + 1) / 2, "DFA", "H", Class="fGn");
            testCase.verifyEqual(hW, 0.8, "AbsTol", 0.05);
            testCase.verifyEqual(hD, 0.8, "AbsTol", 0.05);
        end
    end
end
