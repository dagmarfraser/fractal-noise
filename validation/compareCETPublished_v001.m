% compareCETPublished_v001.m  Our CET estimates on the same footing as the two published analyses (read-only on the toolbox).
%   Rypdal, Ostvand & Rypdal 2013 (JGR Atmospheres 118, 7046-7062, doi 10.1002/jgrd.50399), section 4.1: monthly CET 1659-2011,
%     deseasonalised by subtracting the climatology (each calendar month's mean), then linear (P1) detrending; H 0.63 (WVA),
%     0.64 (DFA8), 0.66 (FA, P1); best estimate H = 0.65 +/- 0.07 (beta 0.30). Undetrended FA H 0.73; DFA0 H 0.79.
%   Gonzalez-Hidalgo, Pena-Angulo & Begueria 2020 (Cuadernos de Investigacion Geografica 46(2), 345-369, doi 10.18172/cig.4377),
%     p. 364: annual CET 1659-2017; H 0.72 (R/S), 0.78 (corrected R/S, Weron 2002), 0.80 (spectral regression, Robinson 1994).
% Ours: local Whittle over 0.02-1 cycles/yr (monthly) or 0.02-0.45 (annual), with the toolbox's DFA-3 pre-estimate, with and
% without P1 detrending, alpha converted with convertScalingExponent (fGn).
here = fileparts(mfilename("fullpath"));  addpath(fullfile(here, "..", "toolbox"));
warning("off", "estimateWhittleAlpha:peakInBand");
fprintf("%-44s %8s %6s %8s %6s\n", "series", "alpha", "H", "DFA3 a", "H");
for span = [1659 2011; 1659 2017; 1659 2025]'
    cet = cetNoiseWrapper(LoadOnly=true, FirstYear=span(1), LastYear=span(2));
    an = cet.value;  for m = 1:12, s = cet.month == m;  an(s) = an(s) - mean(an(s)); end
    t = (1:numel(an))';  anP1 = an - polyval(polyfit(t, an, 1), t);
    yrs = unique(cet.year);  ann = arrayfun(@(y) mean(cet.value(cet.year == y)), yrs);
    tY = (1:numel(ann))';  annP1 = ann - polyval(polyfit(tY, ann, 1), tY);
    rows = {"monthly anomalies", an, 12, [0.02 1]; "monthly anomalies, P1 detrended", anP1, 12, [0.02 1];
            "annual means", ann, 1, [0.02 0.45]; "annual means, P1 detrended", annP1, 1, [0.02 0.45]};
    for i = 1:size(rows, 1)
        [a, ~, info] = estimateWhittleAlpha(rows{i, 2}, rows{i, 3}, rows{i, 4}(1), rows{i, 4}(2));
        fprintf("%d-%d %-34s %8.3f %6.3f %8.3f %6.3f\n", span(1), span(2), rows{i, 1}, a, ...
            convertScalingExponent(a, "alpha", "H", Class="fGn"), info.alphaDFA, convertScalingExponent(info.alphaDFA, "alpha", "H", Class="fGn"));
    end
end
