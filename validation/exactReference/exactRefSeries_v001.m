function [X, info] = exactRefSeries_v001(N, alpha, nRep, opts)
%EXACTREFSERIES_V001 Exact reference series with a known spectral exponent.
%   [X, info] = exactRefSeries_v001(N, alpha, nRep) returns k-fold cumulative
%   sums of exact fGn (fgnDaviesHarte_v001), with alpha = (2H - 1) + 2k and
%   H in (0,1). k = 0 gives stationary fGn (alpha in (-1,1)); k = 1 is fBm
%   (alpha in (1,3)); k = 2, 3 reach alpha in (3,7).
%   Odd-integer alpha (1, 3, 5, 7) is a degenerate boundary (H -> 0 or 1)
%   with no exact stationary-increment representative: it errors loudly
%   rather than silently nudging alpha. Use alpha +/- a small offset and
%   report it.
%   info: struct with H, k, alpha and the analytic spectral density handle
%   info.S(lambda), lambda in (0, pi] rad/sample (pseudo-spectrum for k > 0).
arguments
    N (1,1) double {mustBeInteger, mustBeGreaterThan(N, 1)}
    alpha (1,1) double {mustBeGreaterThan(alpha, -1), mustBeLessThan(alpha, 7)}
    nRep (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.Stream = RandStream.getGlobalStream
end
k = floor((alpha + 1) / 2);
H = (alpha + 1 - 2 * k) / 2;
if H <= 0 || H >= 1 || abs(H) < 1e-12
    error('exactRefSeries:oddAlpha', ['alpha = %g is an odd-integer boundary (H = %g): ' ...
        'no exact fGn-based representative. Use e.g. alpha = %g or %g and report the offset.'], ...
        alpha, H, alpha - 0.01, alpha + 0.01);
end
X = fgnDaviesHarte_v001(N, H, nRep, Stream=opts.Stream);
for j = 1:k; X = cumsum(X, 1); end
info = struct('H', H, 'k', k, 'alpha', alpha, 'S', @(lam) fgnSpectrum(lam, H) ./ (2 - 2*cos(lam)).^k);
end

function S = fgnSpectrum(lam, H)
% Exact fGn spectral density (Sinai 1976; e.g. Beran 1994 eq. 2.17), unit variance,
% normalised so that integral over (-pi, pi] equals 1. Infinite sum truncated at |j| <= J
% with an integral tail correction.
J  = 200;  lam = lam(:);
cH = gamma(2*H + 1) * sin(pi*H) / (2*pi);
j  = -J:J;
s  = sum(abs(lam + 2*pi*j).^(-2*H - 1), 2);
tail = ((2*pi*(J + 0.5) + lam).^(-2*H) + (2*pi*(J + 0.5) - lam).^(-2*H)) / (2*pi*2*H);
S  = 2 * cH * (1 - cos(lam)) .* (s + tail);
end
