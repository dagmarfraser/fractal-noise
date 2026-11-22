% checkWhittleBandEdge_v001.m  Local Whittle bias against the upper band edge at the monthly-CET length (N 4404, fs 12, Nyquist 6).
% Exact references (fGn / fBm by exactRefSeries_v001), 50 series per alpha, master seed 1, peak check off. Bands in cycles per
% year: 0.02-5, 0.02-2, 0.02-1 and the cetNoiseWrapper default 0.5-5. Read-only on the toolbox.
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "toolbox"), fullfile(here, "exactReference"));
rng(1);  warning("off", "all");
for a = [-0.5 0.9 1.5 2.5]
    X = exactRefSeries_v001(4404, a, 50);
    for b = [0.02 5; 0.02 2; 0.02 1; 0.5 5]'
        e = arrayfun(@(r) estimateWhittleAlpha(X(:, r), 12, b(1), b(2), PeakCheck=false), 1:50);
        fprintf("alpha %5.2f band %.2f-%.1f: bias %+.3f sd %.3f\n", a, b(1), b(2), mean(e) - a, std(e));
    end
end
