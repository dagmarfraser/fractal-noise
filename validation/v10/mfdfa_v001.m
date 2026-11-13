function [h, dh, info] = mfdfa_v001(x, options)
%MFDFA_V001 Multifractal DFA (Kantelhardt 2002; Ihlen 2012 layout) with optional prewhitening.
%   [h, dh, info] = mfdfa_v001(x) returns generalised Hurst exponents h(q) for q = -5:5 and the width
%   dh = h(qmin) - h(qmax). A monofractal (linear Gaussian) series has flat h(q), dh ~ 0 up to
%   finite-size spread; calibrate that spread on exact references (runV10Null_v001), never assume 0.
%   Steep spectra: DFA order m only resolves h < m + 1 (Kiyono 2015). Use K (difference K times
%   first; every h(q) shifts by -K, dh is unchanged) to bring the series into the fGn range.
%   Options: Q (default -5:5), Order m (default 2), K (default 0), Scales (default 15 log-spaced,
%   16 .. N/5 after differencing). Segments are taken from both ends (2*floor(N/s) per scale).
%   Errors on non-finite or constant input, or fewer than 4 usable scales.
arguments
    x (:,1) double
    options.Q (1,:) double = -5:5
    options.Order (1,1) double {mustBeInteger, mustBePositive} = 2
    options.K (1,1) double {mustBeInteger, mustBeNonnegative} = 0
    options.Scales (1,:) double = []
end
if ~all(isfinite(x)), error('mfdfa:nonFinite', 'x contains NaN or Inf.'); end
for j = 1:options.K, x = diff(x); end
if std(x) == 0, error('mfdfa:constant', 'x is constant after differencing.'); end
N = numel(x);  q = options.Q;  m = options.Order;
s = options.Scales;
if isempty(s), s = unique(round(logspace(log10(16), log10(floor(N / 5)), 15))); end
s = s(s >= m + 2 & s <= floor(N / 4));
if numel(s) < 4, error('mfdfa:scales', 'Only %d usable scales for N = %d.', numel(s), N); end
Y = cumsum(x - mean(x));
Fq = NaN(numel(q), numel(s));
for si = 1:numel(s)
    n = s(si);  ns = floor(N / n);  t = (1:n)';
    [Qb, ~] = qr(t .^ (0:m), 0);                         % orthonormal polynomial basis, order m
    segs = [reshape(Y(1:ns*n), n, ns), reshape(Y(N-ns*n+1:N), n, ns)];
    F2 = mean((segs - Qb * (Qb' * segs)) .^ 2, 1);      % detrended variance per segment
    F2 = F2(F2 > 0);
    for qi = 1:numel(q)
        if q(qi) == 0, Fq(qi, si) = exp(0.5 * mean(log(F2)));
        else,           Fq(qi, si) = mean(F2 .^ (q(qi) / 2)) ^ (1 / q(qi)); end
    end
end
h = NaN(1, numel(q));
for qi = 1:numel(q)
    p = polyfit(log(s), log(Fq(qi, :)), 1);  h(qi) = p(1);
end
dh = h(1) - h(end);
info = struct('q', q, 'scales', s, 'Fq', Fq, 'order', m, 'K', options.K, 'h2', h(q == 2));
end
