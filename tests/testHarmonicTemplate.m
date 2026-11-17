classdef testHarmonicTemplate < matlab.unittest.TestCase
%TESTHARMONICTEMPLATE Unit tests for harmonicTemplate.m (new in v2.0.0).
%
% Toolbox: fractalnoise
% Covers: recovery of known sinusoids, harmonics above K left in, the identity x = resid + fit, end-anchored coverage,
% parity with the PowerLawSimulationPreReg runner's templateSubtract_local (copied below), channel independence, a
% non-Hz unit (monthly data), f0 estimation, and the error paths.

    methods (TestClassSetup)
        function addToolbox(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
        end
    end

    methods (Test)
        function testRemovesKnownFundamental(testCase)
            rng(1);  fs = 100;  t = (0:1/fs:40)';  noise = 0.2 * randn(size(t));
            x = 3 * sin(2*pi*0.5*t + 0.4) + noise;
            [r, ~, info] = harmonicTemplate(x, fs, 0.5, 1);
            testCase.verifyLessThan(abs(std(r) / std(noise(info.index(1):info.index(2))) - 1), 0.10);
        end

        function testRemovesThreeHarmonicsWhenKIsThree(testCase)
            rng(2);  fs = 120;  t = (0:1/fs:60)';  noise = 0.3 * randn(size(t));
            x = 2*sin(2*pi*0.6*t) + 1.0*sin(2*pi*1.2*t + 1) + 0.7*cos(2*pi*1.8*t + 2) + noise;
            [r, ~, info] = harmonicTemplate(x, fs, 0.6, 3);
            testCase.verifyLessThan(abs(std(r) / std(noise(info.index(1):info.index(2))) - 1), 0.10);
        end

        function testHarmonicsAboveKStayInTheResidual(testCase)
            rng(3);  fs = 120;  t = (0:1/fs:60)';  noise = 0.1 * randn(size(t));
            x = 2*sin(2*pi*0.6*t) + 1.0*sin(2*pi*1.2*t + 1) + noise;   % second harmonic present, K = 1 leaves it
            [r, ~, info] = harmonicTemplate(x, fs, 0.6, 1);
            testCase.verifyGreaterThan(std(r), 5 * std(noise(info.index(1):info.index(2))));
        end

        function testResidualPlusFitIsTheRecord(testCase)
            rng(4);  x = cumsum(randn(2000, 2)) + 5 * sin(2*pi*0.5*(0:1999)'/100);
            [r, f, info] = harmonicTemplate(x, 100, 0.5, 2);
            testCase.verifyEqual(r + f, x(info.index(1):info.index(2), :), "AbsTol", 1e-10);
        end

        function testEndAnchoredCoversTheLastSample(testCase)
            rng(5);  x = randn(1350, 1);   % window 400 samples, hop 200: windows end at 1200 without the anchor, an uncovered tail of 150 > the 100-sample crop
            [~, ~, a] = harmonicTemplate(x, 100, 1, 1, EndAnchored=false);
            [~, ~, b] = harmonicTemplate(x, 100, 1, 1, EndAnchored=true);
            testCase.verifyFalse(a.Covered);
            testCase.verifyTrue(b.Covered);
            testCase.verifyEqual(b.nWindows, a.nWindows + 1);
        end

        function testMatchesProjectRunnerForTwoAxes(testCase)
            rng(6);  fs = 240;  N = 2400;  t = (0:N-1)' / fs;  f0 = 0.5;
            x = 40*cos(2*pi*f0*t) + cumsum(randn(N, 1)) * 0.05;  y = 20*sin(2*pi*f0*t + 0.3) + cumsum(randn(N, 1)) * 0.05;
            [r, ~, info] = harmonicTemplate([x y], fs, f0, 4);
            [rx, ry] = runnerTemplateSubtract(x, y, fs, f0, 4, 4, true);
            testCase.verifyEqual(r(:, 1), rx, "AbsTol", 1e-9);
            testCase.verifyEqual(r(:, 2), ry, "AbsTol", 1e-9);
            testCase.verifyEqual(size(r, 1), numel(rx));
            testCase.verifyEqual(info.index, [round(info.WindowSamples/4), round(info.WindowSamples/4) + numel(rx) - 1]);
        end

        function testChannelsAreIndependent(testCase)
            rng(7);  fs = 100;  t = (0:1/fs:40)';  x = sin(2*pi*0.5*t) + 0.2*randn(size(t));  y = 2*cos(2*pi*0.5*t) + 0.3*randn(size(t));
            both = harmonicTemplate([x y], fs, 0.5, 2);
            testCase.verifyEqual(both(:, 1), harmonicTemplate(x, fs, 0.5, 2), "AbsTol", 1e-10);
            testCase.verifyEqual(both(:, 2), harmonicTemplate(y, fs, 0.5, 2), "AbsTol", 1e-10);
        end

        function testWorksInAnotherUnitMonthlyData(testCase)
            rng(8);  fs = 12;  m = (0:12*100-1)';  noise = 0.4 * randn(size(m));   % 100 years, fs in cycles per year
            x = 8 * sin(2*pi*m/12) + 1.5 * sin(4*pi*m/12 + 1) + noise;
            [r, ~, info] = harmonicTemplate(x, fs, 1, 2, Window=10);
            testCase.verifyLessThan(abs(std(r) / std(noise(info.index(1):info.index(2))) - 1), 0.10);
        end

        function testEstimatesF0InsideTheGivenBand(testCase)
            rng(9);  fs = 100;  t = (0:1/fs:60)';  x = 3*sin(2*pi*0.7*t) + 0.5*randn(size(t));
            [~, ~, info] = harmonicTemplate(x, fs, [], 1, F0Band=[0.2 3]);
            testCase.verifyLessThan(abs(info.f0 - 0.7) / 0.7, 0.03);
        end

        function testErrorPaths(testCase)
            x = randn(500, 1);
            testCase.verifyError(@() harmonicTemplate(x, 100, [], 1), "harmonicTemplate:f0Band");
            testCase.verifyError(@() harmonicTemplate(x, 100, 0.01, 1), "harmonicTemplate:tooShort");   % 4 cycles of 0.01 Hz = 400 s
            testCase.verifyError(@() harmonicTemplate(x, 100, 30, 3), "harmonicTemplate:nyquist");
            testCase.verifyError(@() harmonicTemplate([x; NaN], 100, 1, 1), "harmonicTemplate:nonFinite");
        end
    end
end

function [resX, resY] = runnerTemplateSubtract(x, y, fs, f0, nH, nCW, endAnchored)
% Copied from PowerLawSimulationPreReg/src/runLoopClosureFftnoise_v028.m (templateSubtract_local), returning the cropped residual.
    N   = numel(x);
    win = round(nCW/f0*fs);
    hop = round(win/2);
    win = min(win, N);
    win = max(win, round(2/f0*fs));
    resX = zeros(N,1); resY = zeros(N,1); ww = zeros(N,1);
    starts = 1:hop:(N - win + 1);
    if win > N, starts = []; end
    if endAnchored && ~isempty(starts) && starts(end) + win - 1 < N, starts(end+1) = N - win + 1; end
    for s = starts
        e   = s + win - 1;
        tw  = (0:win-1)'/fs;
        han = 0.5*(1 - cos(2*pi*(0:win-1)'/(win-1)));
        D   = [ones(win,1), tw];
        for h = 1:nH
            D(:,end+1) = cos(2*pi*h*f0*tw);
            D(:,end+1) = sin(2*pi*h*f0*tw);
        end
        Dw        = D .* han;
        bX        = Dw \ (x(s:e) .* han);
        bY        = Dw \ (y(s:e) .* han);
        resX(s:e) = resX(s:e) + (x(s:e) - D*bX) .* han;
        resY(s:e) = resY(s:e) + (y(s:e) - D*bY) .* han;
        ww(s:e)   = ww(s:e) + han;
    end
    ok = ww > 0;
    resX(ok) = resX(ok) ./ ww(ok);
    resY(ok) = resY(ok) ./ ww(ok);
    cl   = max(round(win/4), 1);
    cs   = cl;
    ce   = min(max(N - cl, cs + 100), N);
    resX = resX(cs:ce);
    resY = resY(cs:ce);
end
