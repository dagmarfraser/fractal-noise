% checkShapeNoiseEndMatched_v001.m  shapeNoise "legacy" (1.x, xu carrier) against "endMatched" (v2.0.0 default).
% Reference: noiseXu, default Phi, 100 seeds per cell, alpha 0.5:0.5:6, band 1-20 Hz (estimateWhittleAlpha).
% Configurations: fs 133 N 3000 (the example) and fs 240 N 1000 (a short trial).
% Each surrogate is scored against the series it was built from: legacy against the whole reference,
% endMatched against its fairShare segment (target alpha = the segment's own, alpha = []).
%   gap   = alpha(surrogate) - alpha(source): median / p10 / p90 / share below -0.5
%   d1,d2 = std(diff(surrogate, k)) / std(diff(source, k)), k = 1, 2: median (derivative noise)
%   endMatched also: median L/N and share with no accepted draw (NaN, counted, not scored).
% Whittle and shapeNoise failures are counted and printed, never hidden.
cfg = [133 3000; 240 1000];  alphas = 0.5:0.5:6;  nRep = 100;  fLow = 1;  fHigh = 20;
w = warning("off", "shapeNoise:legacyConstruction");  cw = onCleanup(@() warning(w));
warning("off", "shapeNoise:noAcceptedDraw");
q = @(v, p) v(min(numel(v), max(1, round(p * numel(v)))));
for ci = 1:size(cfg, 1)
    fs = cfg(ci, 1);  N = cfg(ci, 2);  nFail = 0;
    fprintf("\nfs %d, N %d\n%-5s | %-44s | %-44s | %-12s\n", fs, N, "alpha", ...
        "legacy: gap med/p10/p90/sh<-0.5, d1, d2", "endMatched: gap med/p10/p90/sh<-0.5, d1, d2", "L/N, noAcc");
    for ai = 1:numel(alphas)
        a = alphas(ai);  G = nan(nRep, 2);  D1 = nan(nRep, 2);  D2 = nan(nRep, 2);  Lf = nan(nRep, 1);  noAcc = 0;
        for r = 1:nRep
            rng(30000 * ci + 100 * ai + r);  ref = noiseXu(N, a, 1, fs);
            try
                s = shapeNoise(ref, fs, a, Method="legacy");
                G(r, 1) = estimateWhittleAlpha(s, fs, fLow, fHigh) - estimateWhittleAlpha(ref, fs, fLow, fHigh);
                D1(r, 1) = std(diff(s)) / std(diff(ref));  D2(r, 1) = std(diff(s, 2)) / std(diff(ref, 2));
                [s, info] = shapeNoise(ref, fs, []);  Lf(r) = info.Lfrac;
                if ~info.accepted, noAcc = noAcc + 1;  continue, end
                seg = ref(info.i0:info.i1);
                G(r, 2) = estimateWhittleAlpha(s, fs, fLow, fHigh) - info.alphaTarget;
                D1(r, 2) = std(diff(s)) / std(diff(seg));  D2(r, 2) = std(diff(s, 2)) / std(diff(seg, 2));
            catch err, nFail = nFail + 1;  fprintf("FAIL alpha %.1f rep %d: %s\n", a, r, err.message);
            end
        end
        out = strings(1, 2);
        for m = 1:2
            v = sort(G(~isnan(G(:, m)), m));
            out(m) = sprintf("%5.2f/%5.2f/%5.2f/%.2f, %4.2f, %4.2f", median(v), q(v, 0.1), q(v, 0.9), mean(v < -0.5), ...
                median(D1(:, m), "omitnan"), median(D2(:, m), "omitnan"));
        end
        fprintf("%-5.1f | %-44s | %-44s | %4.2f, %.2f\n", a, out(1), out(2), median(Lf, "omitnan"), noAcc / nRep);
    end
    fprintf("Failures (Whittle or shapeNoise): %d of %d seeds\n", nFail, numel(alphas) * nRep);
end
