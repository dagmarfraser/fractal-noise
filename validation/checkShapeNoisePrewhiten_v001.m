% checkShapeNoisePrewhiten_v001.m  shapeNoise Prewhiten=0 against Prewhiten="auto" across alpha and record length.
% Gap = estimateWhittleAlpha(surrogate) - estimateWhittleAlpha(reference), "xu" carrier, same carrier for both
% constructions (rng reset before each shapeNoise call), 100 seeds per cell, band 1-20 Hz.
% Configurations: fs 133 N 3000 (the example) and fs 240 N 1000 (a short trial). Reference: noiseXu, default Phi.
% Whittle failures are counted and printed, never hidden.
cfg = [133 3000; 240 1000];  alphas = 0.5:0.5:6;  nRep = 100;  fLow = 1;  fHigh = 20;
for ci = 1:size(cfg, 1)
    fs = cfg(ci, 1);  N = cfg(ci, 2);  nFail = 0;
    fprintf("\nfs %d, N %d: gap median / p10 / p90 / share<-0.5\n%-6s | %-2s | %-32s | %-32s\n", fs, N, "alpha", "k", "Prewhiten=0", "Prewhiten=""auto""");
    for ai = 1:numel(alphas)
        a = alphas(ai);  D = nan(nRep, 2);
        for r = 1:nRep
            rng(20000 * ci + 100 * ai + r);  ref = noiseXu(N, a, 1, fs);  s = rng;
            try
                aRef = estimateWhittleAlpha(ref, fs, fLow, fHigh);
                rng(s);  D(r, 1) = estimateWhittleAlpha(shapeNoise(ref, fs, a, Prewhiten=0), fs, fLow, fHigh) - aRef;
                rng(s);  D(r, 2) = estimateWhittleAlpha(shapeNoise(ref, fs, a), fs, fLow, fHigh) - aRef;
            catch err, nFail = nFail + 1; fprintf("FAIL alpha %.1f rep %d: %s\n", a, r, err.message);
            end
        end
        out = strings(1, 2);
        for m = 1:2
            v = sort(D(~isnan(D(:, m)), m));  q = @(p) v(min(numel(v), max(1, round(p * numel(v)))));
            out(m) = sprintf("%5.2f / %5.2f / %5.2f / %.2f", median(v), q(0.1), q(0.9), mean(v < -0.5));
        end
        fprintf("%-6.1f | %-2d | %-32s | %-32s\n", a, min(3, max(0, floor((a + 1) / 2))), out(1), out(2));
    end
    fprintf("Whittle failures: %d of %d seeds\n", nFail, numel(alphas) * nRep);
end
