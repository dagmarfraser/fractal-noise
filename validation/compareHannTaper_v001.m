function T = compareHannTaper_v001()
%COMPAREHANNTAPER_V001 estimateWhittleAlpha with its Hann taper against the same code without it (a temporary copy).
%   Cases: exact references (no break) at alpha 0.5-6.5 including the edges of each differencing order; the same with the
%   differencing order forced one too low (as when the DFA pre-estimate is wrong); noiseXu with its default break; a 10 Hz peak.
%   fs 240, N 3000, band 1-20 Hz, 50 series per cell. Returns a table and prints it. Read-only on the toolbox.
here = fileparts(mfilename("fullpath"));  tb = fullfile(here, "..", "toolbox");
addpath(tb, fullfile(here, "exactReference"));
d = tempname;  mkdir(d);  cleanup = onCleanup(@() rmdir(d, "s"));
s = fileread(fullfile(tb, "estimateWhittleAlpha.m"));
s2 = strrep(s, "w = 0.5 * (1 - cos(2*pi*(0:N-1)' / N));", "w = ones(N, 1);");
assert(~strcmp(s, s2), "taper line not found");
fid = fopen(fullfile(d, "estimateWhittleAlpha.m"), "w");  fwrite(fid, s2);  fclose(fid);
rows = {};  n = 50;  fs = 240;  t = (0:2999)' / fs;
cases = {"exact", [0.5 1.95 2.05 2.5 2.95 3.95 4.05 4.5 5.95 6.5]; "exact, k one too low", [2.5 4.5 6.5]; ...
         "noiseXu Phi 0.99 (break)", [2.5 4.5 6.5]; "exact + 10 Hz peak", [2.5 4.5]};
for c = 1:size(cases, 1)
    for a = cases{c, 2}
        rng(round(1000 * a) + c);
        if startsWith(cases{c, 1}, "noiseXu")
            X = cell2mat(arrayfun(@(~) noiseXu(3000, a, 1, fs), 1:n, "UniformOutput", false));  truth = NaN;
        else
            X = exactRefSeries_v001(3000, a, n);  truth = a;
            if contains(cases{c, 1}, "peak"), X = X + 2 * std(X(:)) * sin(2 * pi * 10 * t); end
        end
        kTrue = min(max(floor((a + 1) / 2), 0), 3);
        for variant = ["Hann", "none"]
            if variant == "none", addpath(d, "-begin"); else, rmpath(d); addpath(tb, "-begin"); end
            e = NaN(n, 1);  nErr = 0;
            for r = 1:n
                try
                    if contains(cases{c, 1}, "too low"), e(r) = estimateWhittleAlpha(X(:, r), fs, 1, 20, K=max(kTrue - 1, 0));
                    else, e(r) = estimateWhittleAlpha(X(:, r), fs, 1, 20); end
                catch
                    nErr = nErr + 1;
                end
            end
            rows(end + 1, :) = {string(cases{c, 1}), a, variant, median(e, "omitnan") - truth, std(e, "omitnan"), nErr, median(e, "omitnan")}; %#ok<AGROW>
        end
        rmpath(d);
    end
end
T = cell2table(rows, "VariableNames", ["case", "alpha", "taper", "medianBias", "sd", "errors", "medianEstimate"]);
disp(T)
end
