% compareScalingExponents_v001.m  Hurst exponent, DFA exponent and Higuchi (Hausdorff) fractal dimension as spectral alpha.
% Read-only on the toolbox. Master seed 1.
%
% Conversions (Gaussian, self-affine, single scaling regime; e.g. Mandelbrot & Van Ness 1968, Voss 1985, Eke et al. 2002):
%   fGn-like series (stationary, -1 < alpha < 1):  alpha = 2H - 1,  d = H - 1/2
%   fBm-like series (1 < alpha < 3):               alpha = 2H + 1  (H of the increments' sum, i.e. of the series itself)
%   DFA fluctuation exponent a (any order):        alpha = 2a - 1  (a = H for fGn, H + 1 for fBm)
%   Higuchi fractal dimension D of a trace:        D = 2 - H of that trace, so alpha = 5 - 2D for an fBm-like trace;
%                                                  for an fGn-like series D is taken on its profile (cumsum), giving alpha = 3 - 2D
%   Rescaled range (R/S) estimates H of an fGn-like series, alpha = 2H - 1 (known small-sample upward bias at low H).
% Part A checks each route on exact references (alpha known) at the length of monthly CET 1659-2025.
% Part B applies them to CET monthly-climatology anomalies and to annual-mean anomalies, and converts published CET / land
% temperature H values to alpha. Published values are from the abstracts named below; method details not verified here.
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "toolbox"), fullfile(here, "exactReference"));
rng(1);

% Part A: exact references
cet = cetNoiseWrapper(LoadOnly=true, FirstYear=1659, LastYear=2025);
n = numel(cet.value);  nRep = 50;
alphas = [-0.5 0 0.3 0.6 0.9 1.5 2 2.5];
rows = {};
for a = alphas
    X = exactRefSeries_v001(n, a, nRep);
    est = NaN(nRep, 4);
    for r = 1:nRep
        x = X(:, r);
        est(r, :) = allRoutes(x, a < 1);
    end
    rows(end + 1, :) = [{a}, num2cell(mean(est, 1, "omitnan") - a), num2cell(std(est, 0, 1, "omitnan"))]; %#ok<AGROW>
end
A = cell2table(rows, "VariableNames", ["alphaTrue", "biasWhittle", "biasDFA1", "biasRS", "biasHiguchi", ...
    "sdWhittle", "sdDFA1", "sdRS", "sdHiguchi"]);
fprintf("\nPart A: exact references, N = %d, %d series each. Bias and SD of alpha from each route (R/S only for alpha < 1).\n", n, nRep);
disp(A)

% Part B: CET
an = cet.value;  for mo = 1:12, sel = cet.month == mo;  an(sel) = an(sel) - mean(an(sel)); end
yrs = unique(cet.year);  annual = arrayfun(@(y) mean(cet.value(cet.year == y)), yrs);  annual = annual - mean(annual);
eM = allRoutes(an, true);  eA = allRoutes(annual, true);
fprintf("\nPart B: CET 1659-2025 as alpha (H = (alpha + 1)/2 for fGn-like series)\n");
fprintf("  %-34s %8s %8s %8s %8s\n", "series", "Whittle", "DFA-1", "R/S", "Higuchi");
fprintf("  %-34s %8.2f %8.2f %8.2f %8.2f\n", "monthly anomalies (N " + numel(an) + ")", eM);
fprintf("  %-34s %8.2f %8.2f %8.2f %8.2f\n", "annual-mean anomalies (N " + numel(annual) + ")", eA);
fprintf("  (Whittle band: monthly 0.02-5 cycles/yr; annual 0.02-0.45 cycles/yr)\n");

pub = {"Gonzalez Hidalgo et al. 2020, annual CET, lower", 0.72; "Gonzalez Hidalgo et al. 2020, annual CET, upper", 0.80;
       "Rypdal et al. 2013, local/regional land", 0.65; "Rypdal et al. 2013, global land", 0.75; "Rypdal et al. 2013, ocean-dominated", 1.0};
fprintf("\n  Published H converted (fGn-like assumption): alpha = 2H - 1, d = H - 1/2, profile D = 2 - H\n");
for i = 1:size(pub, 1)
    H = pub{i, 2};
    fprintf("  %-50s H %.2f -> alpha %.2f, d %.2f, D %.2f\n", pub{i, 1}, H, 2 * H - 1, H - 0.5, 2 - H);
end

function e = allRoutes(x, fgnLike)
% alpha from local Whittle (toolbox), DFA-1, R/S (fGn-like only) and Higuchi D.
x = x(:);  n = numel(x);
if n > 1000, fs = 12; band = [0.02 5]; else, fs = 1; band = [0.02 0.45]; end
w = warning("off", "all");  c = onCleanup(@() warning(w));
e = NaN(1, 4);
e(1) = estimateWhittleAlpha(x, fs, band(1), band(2), PeakCheck=false);
e(2) = 2 * dfa1(x) - 1;
if fgnLike, e(3) = 2 * rescaledRange(x) - 1; end
if fgnLike, e(4) = 3 - 2 * higuchiD(cumsum(x - mean(x))); else, e(4) = 5 - 2 * higuchiD(x); end
end

function a = dfa1(x)
% DFA of order 1: slope of log F(s) on log s, s from 10 to N/4, non-overlapping windows, linear detrend.
y = cumsum(x - mean(x));  n = numel(y);
s = unique(round(logspace(log10(10), log10(n / 4), 20)));  F = NaN(size(s));
for i = 1:numel(s)
    m = floor(n / s(i));  Y = reshape(y(1:m * s(i)), s(i), m);  t = (1:s(i))';  G = [ones(s(i), 1) t];
    R = Y - G * (G \ Y);  F(i) = sqrt(mean(R(:) .^ 2));
end
p = polyfit(log(s), log(F), 1);  a = p(1);
end

function H = rescaledRange(x)
% Classic R/S: mean over non-overlapping blocks of range of cumulative deviation / SD, slope over block sizes 10 to N/4.
n = numel(x);  s = unique(round(logspace(log10(10), log10(n / 4), 20)));  RS = NaN(size(s));
for i = 1:numel(s)
    m = floor(n / s(i));  Y = reshape(x(1:m * s(i)), s(i), m);
    Z = cumsum(Y - mean(Y, 1), 1);  RS(i) = mean((max(Z, [], 1) - min(Z, [], 1)) ./ std(Y, 1, 1));
end
p = polyfit(log(s), log(RS), 1);  H = p(1);
end

function D = higuchiD(x)
% Higuchi (1988) fractal dimension, k = 1..16.
n = numel(x);  K = 1:16;  L = NaN(size(K));
for j = 1:numel(K)
    k = K(j);  Lm = NaN(k, 1);
    for m = 1:k
        idx = m:k:n;  nI = numel(idx) - 1;
        Lm(m) = sum(abs(diff(x(idx)))) * (n - 1) / (nI * k) / k;
    end
    L(j) = mean(Lm);
end
p = polyfit(log(K), log(L), 1);  D = -p(1);
end
