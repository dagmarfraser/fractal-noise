% exploreHadCET_v001.m  What the toolbox gets from the bundled monthly Central England Temperature series (read-only).
% Annual cycle removed with harmonicTemplate (K from harmonicContent), then the spectral exponent of the residual over a
% short-term band (0.5-5 cycles/year, the cetNoiseWrapper default) and a long-term band (0.02-0.5 cycles/year: 2 to 50 years),
% by estimateWhittleAlpha (local Whittle; d = alpha/2), with its DFA-3 pre-estimate and the peak check, for 1659-2025 and 1900-2025.
here = fileparts(mfilename("fullpath"));  addpath(fullfile(here, "..", "toolbox"));
for span = [1659 2025; 1900 2025]'
    cet = cetNoiseWrapper(LoadOnly=true, FirstYear=span(1), LastYear=span(2));
    [K, hc] = harmonicContent(cet.value, 12, 1);
    res = harmonicTemplate(cet.value, 12, 1, K);
    fprintf("\n%d-%d: %d months, K = %d (ratios %s), residual SD %.2f degC\n", span(1), span(2), numel(cet.value), K, mat2str(round(hc.ratio(1:4), 1)), std(res));
    for band = [0.5 5; 0.02 0.5]'
        lastwarn("");
        [a, k, info] = estimateWhittleAlpha(res, 12, band(1), band(2));
        [~, wid] = lastwarn;
        fprintf("  band %.2f-%.2f cyc/yr: alpha %.3f (d %.3f), k %d, bins %d, DFA alpha %.3f, peak ratio %.1f / threshold %.1f at %.3f cyc/yr %s\n", ...
            band(1), band(2), a, a / 2, k, info.nBins, info.alphaDFA, info.peakRatio, info.peakThreshold, info.peakHz, wid);
    end
end

% For comparison with the climate literature: anomalies by removing each calendar month's mean over the span (no windowed template,
% so variability slower than the template's 4-year windows is kept).
for span = [1659 2025; 1900 2025]'
    cet = cetNoiseWrapper(LoadOnly=true, FirstYear=span(1), LastYear=span(2));
    an = cet.value;  for mo = 1:12, sel = cet.month == mo;  an(sel) = an(sel) - mean(an(sel)); end
    fprintf("\n%d-%d, monthly-climatology anomalies: SD %.2f degC\n", span(1), span(2), std(an));
    for band = [0.5 5; 0.02 0.5]'
        lastwarn("");  [a, k, info] = estimateWhittleAlpha(an, 12, band(1), band(2));  [~, wid] = lastwarn;
        fprintf("  band %.2f-%.2f cyc/yr: alpha %.3f (d %.3f, H %.3f), k %d, DFA alpha %.3f %s\n", band(1), band(2), a, a / 2, (a + 1) / 2, k, info.alphaDFA, wid);
    end
end
