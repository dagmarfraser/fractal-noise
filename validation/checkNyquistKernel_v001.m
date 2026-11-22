% checkNyquistKernel_v001.m  Is the local Whittle bias near Nyquist (checkWhittleBandEdge_v001.m) a property of the estimator's
% power-law kernel lambda^(-alpha), or of the test series? Read-only on the toolbox. Master seed 1.
%
% The estimator models the recoloured periodogram as G*lambda^(-alpha) (lambda in rad/sample). That is exact only near zero
% frequency. Three generators with the same low-frequency alpha differ near Nyquist:
%   exact   - fGn / k-fold integrated fGn (exactRefSeries_v001): spectrum ~ lambda^(-alpha) at low frequency only;
%   arfima  - ARFIMA(0,d,0), k-fold integrated: spectrum exactly |2 sin(lambda/2)|^(-alpha);
%   xu      - noiseXu with Phi = 1 (no break): the toolbox's own power-law generator.
% Two kernels are compared in a local copy of the estimator (same DFA-3 pre-estimate, k differences, Hann taper, exact
% recolouring): "lambda" (the toolbox's) and "sin", |2 sin(lambda/2)|^(-alpha). The local "lambda" copy must reproduce
% estimateWhittleAlpha exactly (self-check, fails loud).
% Settings: monthly CET length (N 4404, fs 12, bands 0.5-5, 0.02-5, 0.02-1 cycles/yr) and a movement-like case (N 3000,
% fs 240, band 1-20 Hz, upper edge 1/6 of Nyquist). 50 series per cell.
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "toolbox"), fullfile(here, "exactReference"));
rng(1);  warning("off", "all");
cases = {4404, 12, [0.5 5]; 4404, 12, [0.02 5]; 4404, 12, [0.02 1]; 3000, 240, [1 20]};
alphas = [-0.5 0.3 0.9 1.5 2.5 4.5];  gens = ["exact", "arfima", "xu"];  nRep = 50;
rows = {};  nCheck = 0;
for c = 1:size(cases, 1)
    [n, fs, band] = cases{c, :};
    for g = gens
        for a = alphas
            X = makeSeries(g, n, a, fs, nRep);
            e = NaN(nRep, 2);
            for r = 1:nRep
                e(r, 1) = localWhittle(X(:, r), fs, band, "lambda");
                e(r, 2) = localWhittle(X(:, r), fs, band, "sin");
                if r <= 2
                    ref = estimateWhittleAlpha(X(:, r), fs, band(1), band(2), PeakCheck=false);
                    if abs(ref - e(r, 1)) > 1e-6, error("checkNyquistKernel:selfCheck", "Local copy differs from toolbox: %g vs %g.", e(r, 1), ref); end
                    nCheck = nCheck + 1;
                end
            end
            rows(end + 1, :) = {n, fs, sprintf("%g-%g", band), g, a, mean(e(:, 1)) - a, mean(e(:, 2)) - a, std(e(:, 1)), std(e(:, 2))}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, "VariableNames", ["N", "fs", "band", "generator", "alpha", "biasLambda", "biasSin", "sdLambda", "sdSin"]);
fprintf("Self-check: %d series, local lambda kernel equals estimateWhittleAlpha.\n", nCheck);
disp(T)
fprintf("\nLargest |bias| per setting and kernel:\n");
T.absLambda = abs(T.biasLambda);  T.absSin = abs(T.biasSin);
disp(groupsummary(T, ["fs", "band", "generator"], "max", ["absLambda", "absSin"]))
writetable(T, fullfile(here, "checkNyquistKernel_v001.csv"));

function X = makeSeries(g, n, a, fs, nRep)
switch g
    case "exact"
        if mod(a + 1, 2) == 0, a = a + 0.01; end
        X = exactRefSeries_v001(n, a, nRep);
    case "xu"
        X = cell2mat(arrayfun(@(~) noiseXu(n, a, 1, fs, Phi=1), 1:nRep, "UniformOutput", false));
    case "arfima"
        k = floor((a + 1) / 2);  d = a / 2 - k;       % stationary part d in [-0.5, 0.5)
        burn = 20000;  L = n + burn;  nfft = 2 ^ nextpow2(2 * L);
        psi = ones(L, 1);  for j = 2:L, psi(j) = psi(j - 1) * (j - 2 + d) / (j - 1); end
        Psi = fft(psi, nfft);  X = zeros(n, nRep);
        for r = 1:nRep
            y = real(ifft(Psi .* fft(randn(L, 1), nfft)));  y = y(burn + 1:L);
            for j = 1:k, y = cumsum(y); end
            X(:, r) = y;
        end
end
end

function aEst = localWhittle(x, fs, band, kernel)
% Copy of estimateWhittleAlpha's estimate (default options, no peak check) with a choice of spectral kernel.
x = double(x(:));  x = x - mean(x);
aDFA = dfaLocal(x, 3);  kUsed = min(max(floor((aDFA + 1) / 2), 0), 3);
xd = x;  for j = 1:kUsed, xd = diff(xd); end
N = numel(xd);  w = 0.5 * (1 - cos(2*pi*(0:N-1)' / N));
I = abs(fft((xd - mean(xd)) .* w)).^2 / sum(w.^2);
f = (0:N-1)' * fs / N;  m = f >= band(1) & f <= band(2) & f < fs/2;
lam = 2*pi*f(m)/fs;  Ib = I(m) ./ (2 - 2*cos(lam)).^kUsed;
if kernel == "lambda", base = lam; else, base = 2 * sin(lam / 2); end
obj = @(a) log(mean(base.^a .* Ib)) - a * mean(log(base));
aEst = fminbnd(obj, -2, 10, optimset('TolX', 1e-6));
end

function a = dfaLocal(x, order)
% Same as the toolbox's dfaAlphaLocal (DFA of the given order, 15 log-spaced scales from 10 to N/10).
y = cumsum(x - mean(x));  n = numel(y);
s = unique(round(exp(linspace(log(10), log(floor(n/10)), 15))));  F = NaN(size(s));
for i = 1:numel(s)
    nW = floor(n / s(i));  Y = reshape(y(1:nW * s(i)), s(i), nW);
    t = (1:s(i))' / s(i);  [Q, ~] = qr(t .^ (0:order), 0);
    R = Y - Q * (Q' * Y);  F(i) = sqrt(mean(R(:).^2));
end
p = polyfit(log(s), log(F), 1);  a = 2 * p(1) - 1;
end
