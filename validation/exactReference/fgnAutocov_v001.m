function g = fgnAutocov_v001(k, H)
%FGNAUTOCOV_V001 Exact autocovariance of unit-variance fGn at integer lags k.
%   g(k) = 0.5 * (|k+1|^2H - 2|k|^2H + |k-1|^2H)  (Mandelbrot & Van Ness 1968)
k = abs(k);
g = 0.5 * (abs(k + 1).^(2*H) - 2 * k.^(2*H) + abs(k - 1).^(2*H));
end
