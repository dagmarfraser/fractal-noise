function [N, D] = calibratePeakCheck_v002()
%CALIBRATEPEAKCHECK_V002 False alarms and detection of the estimateWhittleAlpha peak check (bin over its neighbours, warning
%   "estimateWhittleAlpha:peakInBand", PeakFalseAlarm 1e-4 per call). The rule is set from theory, not from these runs.
%   N (null, peak-free): exact references (k-fold integrated fGn, no break) at alpha -0.75 to 6.5, noiseXu with its default
%   break (Phi 0.99) and without (Phi 1), and white noise; fs x N in {120, 240, 300 Hz} x {1000, 3000}; band 1-20 Hz; 200 series
%   per cell (100 for noiseXu). Counts the calls that would warn.
%   D (detection): exact references at alpha 0.5, 2.5, 4.5 and 6.5, fs 240 N 3000 and fs 133 N 1000, 50 series, plus a sine at
%   2, 10 or 18 Hz whose amplitude is a fraction of the series' SD. For each: share caught, median Whittle bias, and the largest
%   bias among the series the check misses (the bias that passes silently).
%   Master seed 1. Read-only on the toolbox. Prints both tables and returns them.
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "toolbox"), fullfile(here, "exactReference"));
cfg = [120 1000; 120 3000; 240 1000; 240 3000; 300 1000; 300 3000];
rows = {};
for c = 1:size(cfg, 1)
    fs = cfg(c, 1);  n = cfg(c, 2);
    gens = {"exact", [-0.75 -0.25 0.5 1.5 2.5 3.5 4.5 5.5 6.5], 200; "noiseXu Phi 0.99", [0.5 2.5 4.5 6.5], 100; ...
            "noiseXu Phi 1", [0.5 2.5 4.5 6.5], 100; "white", 0, 200};
    for g = 1:size(gens, 1)
        for a = gens{g, 2}
            rng(1e6 * c + 1e4 * g + round(100 * (a + 1)));
            nRep = gens{g, 3};
            switch gens{g, 1}
                case "exact", X = exactRefSeries_v001(n, a, nRep);
                case "noiseXu Phi 0.99", X = cell2mat(arrayfun(@(~) noiseXu(n, a, 1, fs), 1:nRep, "UniformOutput", false));
                case "noiseXu Phi 1", X = cell2mat(arrayfun(@(~) noiseXu(n, a, 1, fs, Phi=1), 1:nRep, "UniformOutput", false));
                case "white", X = randn(n, nRep);
            end
            hit = false(nRep, 1);  nErr = 0;  worst = 0;
            for r = 1:nRep
                try
                    [~, ~, info] = estimateWhittleAlpha(X(:, r), fs, 1, 20, PeakCheck=false);
                    hit(r) = info.peakRatio > info.peakThreshold;  worst = max(worst, info.peakRatio / info.peakThreshold);
                catch
                    nErr = nErr + 1;
                end
            end
            rows(end + 1, :) = {gens{g, 1}, fs, n, a, nRep, nErr, nnz(hit), worst}; %#ok<AGROW>
        end
    end
end
N = cell2table(rows, "VariableNames", ["generator", "fs", "N", "alpha", "series", "otherErrors", "falseAlarms", "worstRatioOverThreshold"]);
nCalls = sum(N.series - N.otherErrors);
fprintf("NULL: %d calls, %d false alarms (expected %.2f at 1e-4 per call); worst ratio/threshold %.2f\n", nCalls, sum(N.falseAlarms), 1e-4 * nCalls, max(N.worstRatioOverThreshold));
disp(groupsummary(N, "generator", "sum", ["series", "falseAlarms"]));

rows = {};
for cf = [240 3000; 133 1000]'
    fs = cf(1);  n = cf(2);  t = (0:n - 1)' / fs;
    for a = [0.5 2.5 4.5 6.5]
        rng(round(100 * a) + fs);  X = exactRefSeries_v001(n, a, 50);
        for f0 = [2 10 18]
            for amp = [0.0003 0.001 0.003 0.01 0.03 0.1 0.3]
                det = false(50, 1);  bias = NaN(50, 1);
                for r = 1:50
                    x = X(:, r) + amp * std(X(:, r)) * sin(2 * pi * f0 * t + 2 * pi * r / 50);
                    try
                        [e, ~, info] = estimateWhittleAlpha(x, fs, 1, 20, PeakCheck=false, AlphaBounds=[-20 30]);
                        det(r) = info.peakRatio > info.peakThreshold;  bias(r) = e - a;
                    catch
                    end
                end
                missed = ~det & ~isnan(bias);
                worstMissed = NaN;  if any(missed), worstMissed = max(abs(bias(missed))); end
                rows(end + 1, :) = {fs, n, a, f0, amp, mean(det), median(bias, "omitnan"), worstMissed}; %#ok<AGROW>
            end
        end
    end
end
D = cell2table(rows, "VariableNames", ["fs", "N", "alpha", "sineHz", "sineAmpOverSD", "caughtShare", "medianWhittleBias", "worstBiasMissed"]);
disp(D)
end
