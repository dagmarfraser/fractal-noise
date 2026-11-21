function T = checkWhittlePeak_v001()
%CHECKWHITTLEPEAK_V001 What a rhythmic peak inside the band does to estimateWhittleAlpha, and to estimateIRASA on the same series.
%   Exact references (alpha 2.5 and 4.5, fs 240, N 3000, 20 series each, master seed 1) plus a 10 Hz sine whose amplitude is a
%   fraction of the series' SD (0, 0.01, 0.03, 0.1, 0.3, 1). Band 1-20 Hz. estimateIRASA uses Prewhiten = k from the Whittle call
%   on the peak-free series. Errors are counted, never dropped. Returns the table and prints it.
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "toolbox"), fullfile(here, "exactReference"));
fs = 240;  t = (0:2999)' / fs;  rows = {};
for a = [2.5 4.5]
    rng(1);  X = exactRefSeries_v001(3000, a, 20);  k = floor((a + 1) / 2);
    for amp = [0 0.01 0.03 0.1 0.3 1]
        w = NaN(20, 1);  ir = w;  wErr = 0;  iErr = 0;
        for r = 1:20
            x = X(:, r) + amp * std(X(:, r)) * sin(2 * pi * 10 * t);
            try, w(r) = estimateWhittleAlpha(x, fs, 1, 20); catch, wErr = wErr + 1; end
            try, ir(r) = estimateIRASA(x, fs, 1, 20, Prewhiten=k); catch, iErr = iErr + 1; end
        end
        rows(end + 1, :) = {a, amp, median(w, "omitnan"), wErr, median(ir, "omitnan"), iErr}; %#ok<AGROW>
    end
end
T = cell2table(rows, "VariableNames", ["alpha", "sineAmpOverSD", "whittleMedian", "whittleErrors", "irasaMedian", "irasaErrors"]);
disp(T)
end
