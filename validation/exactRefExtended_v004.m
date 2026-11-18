function [x, info] = exactRefExtended_v004(N, alpha)
%EXACTREFEXTENDED_V004 Exact fGn-based reference at any non-odd-integer alpha, outside the (-1, 7) range of exactRefSeries_v001.
%   [x, info] = exactRefExtended_v004(N, alpha) writes alpha = (2H - 1) + 2k with H in (0, 1) and k = floor((alpha + 1) / 2).
%   k > 0: k-fold cumulative sum of exact fGn (as exactRefSeries_v001; k = 4, 5 reach alpha in (7, 11)).
%   k = 0: exact fGn.
%   k < 0: |k|-fold first difference of exact fGn of length N + |k| (alpha in (-5, -1)).
%   The fGn comes from fgnDaviesHarte_v001 (circulant embedding, eigenvalues checked on every call). Summing or differencing is an exact
%   linear operation in exact arithmetic; in double precision a k-fold sum loses high-frequency content to round-off at large k, which
%   this function does not hide: the layer A spectrum check in checkShapeNoiseLayers_v004 measures it.
%   Odd-integer alpha (H -> 0 or 1) has no exact representative and errors ("exactRefExtended:oddAlpha").
%   info: H, k, alpha, and S(lambda), the analytic (pseudo-)spectrum, lambda in (0, pi] rad/sample.
arguments
    N (1, 1) double {mustBeInteger, mustBeGreaterThan(N, 1)}
    alpha (1, 1) double {mustBeGreaterThan(alpha, -5), mustBeLessThan(alpha, 11)}
end
k = floor((alpha + 1) / 2);
H = (alpha + 1 - 2 * k) / 2;
if H < 1e-9 || H > 1 - 1e-9
    error('exactRefExtended:oddAlpha', 'alpha = %g is an odd integer (H = %g): no exact fGn-based representative.', alpha, H);
end
x = fgnDaviesHarte_v001(N + max(-k, 0), H, 1);
if k > 0
    for j = 1:k, x = cumsum(x, 1); end
elseif k < 0
    x = diff(x, -k, 1);
end
info = struct('H', H, 'k', k, 'alpha', alpha, 'S', @(lam) fgnSpectrum(lam, H) ./ (2 - 2 * cos(lam(:))) .^ k);
end

function S = fgnSpectrum(lam, H)
% Exact fGn spectral density, unit variance; copied unchanged from exactRefSeries_v001 (a local function there).
J  = 200;  lam = lam(:);
cH = gamma(2*H + 1) * sin(pi*H) / (2*pi);
j  = -J:J;
s  = sum(abs(lam + 2*pi*j).^(-2*H - 1), 2);
tail = ((2*pi*(J + 0.5) + lam).^(-2*H) + (2*pi*(J + 0.5) - lam).^(-2*H)) / (2*pi*2*H);
S  = 2 * cH * (1 - cos(lam)) .* (s + tail);
end
