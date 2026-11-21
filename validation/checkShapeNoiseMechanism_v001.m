% checkShapeNoiseMechanism_v001.m  Why are some shapeNoise surrogates shallow, and does a differenced-domain swap cure it?
% Hypothesis: |fft(ref)| of a steep, non-periodic series is dominated above the lowest bins by end-mismatch leakage
% (f^-2). The surrogate inherits that magnitude. It stays steep only when the carrier's phases in those bins are the
% coherent phases of the carrier's own end jump; a carrier whose ends happen to nearly match supplies random phases
% there, and the leakage magnitude becomes real f^-2 noise.
% Part A (alpha 4, 300 seeds): shortfall (Whittle surrogate - reference) by tercile of the carrier's end mismatch
%         |c(end) - c(1)| / std(c), and by tercile of the reference's.
% Part B (alpha 3:0.5:5, 100 seeds, same seeds and carriers as checkShapeNoiseAlpha4_v002 Part 2): the current
%         construction against a differenced-domain swap (difference both k times, swap, integrate k times,
%         k = min(3, floor((alpha + 1) / 2))).
% fs 133, N 3000, band 1-20 Hz. Whittle failures are counted and printed, never hidden.
fs = 133;  N = 3000;  fLow = 1;  fHigh = 20;
swap = @(r, c) real(ifft(abs(fft(r)) .* exp(1i * angle(fft(c)))));
jump = @(x) abs(x(end) - x(1)) / std(x);

nA = 300;  d = nan(nA, 1);  jc = nan(nA, 1);  jr = nan(nA, 1);  nFail = 0;
for r = 1:nA
    rng(9000 + r);
    ref = noiseXu(N, 4.0, 1.0, fs);  c = noiseXu(N, 4.0, 1.0, fs);  sur = swap(ref, c);
    jc(r) = jump(c);  jr(r) = jump(ref);
    try, d(r) = estimateWhittleAlpha(sur, fs, fLow, fHigh) - estimateWhittleAlpha(ref, fs, fLow, fHigh);
    catch err, nFail = nFail + 1; fprintf("FAIL A rep %d: %s\n", r, err.message); end
end
fprintf("PART A: alpha 4, %d seeds, Whittle failures %d\n", nA, nFail);
for v = ["carrier", "reference"]
    j = jc;  if v == "reference", j = jr; end
    [~, o] = sort(j);  e = round(linspace(0, nA, 4));
    for t = 1:3
        ix = o(e(t)+1:e(t+1));  dd = d(ix);  dd = dd(~isnan(dd));
        fprintf("%-9s end mismatch tercile %d (%.2f to %.2f SD): median gap %6.2f, share below -0.5: %.2f\n", ...
            v, t, min(j(ix)), max(j(ix)), median(dd), mean(dd < -0.5));
    end
end

alphas = 3:0.5:5;  nRep = 100;  nFailB = 0;
fprintf("\nPART B: surrogate minus reference (Whittle), %d seeds\n%-6s | %-2s | %-30s | %-30s\n", nRep, "target", "k", ...
    "current: med / p10 / share<-0.5", "differenced: med / p10 / share<-0.5");
for ai = 1:numel(alphas)
    a = alphas(ai);  k = min(3, floor((a + 1) / 2));  D = nan(nRep, 2);
    for r = 1:nRep
        rng(7000 + 1000 * ai + r);
        ref = noiseXu(N, a, 1.0, fs);  c = noiseXu(N, a, std(ref), fs);   % same draw order as shapeNoise(ref, fs, a)
        s1 = swap(ref, c);
        s2 = swap(diff(ref, k), diff(c, k));  for q = 1:k, s2 = cumsum([0; s2]); end
        try
            aR = estimateWhittleAlpha(ref, fs, fLow, fHigh);
            D(r, 1) = estimateWhittleAlpha(s1, fs, fLow, fHigh) - aR;
            D(r, 2) = estimateWhittleAlpha(s2, fs, fLow, fHigh) - aR;
        catch err, nFailB = nFailB + 1; fprintf("FAIL B alpha %.1f rep %d: %s\n", a, r, err.message);
        end
    end
    out = strings(1, 2);
    for m = 1:2
        s = sort(D(~isnan(D(:, m)), m));
        out(m) = sprintf("%6.2f / %6.2f / %.2f", median(s), s(max(1, round(0.1 * numel(s)))), mean(s < -0.5));
    end
    fprintf("%-6.1f | %-2d | %-30s | %-30s\n", a, k, out(1), out(2));
end
fprintf("Whittle failures in Part B: %d of %d seeds\n", nFailB, numel(alphas) * nRep);
