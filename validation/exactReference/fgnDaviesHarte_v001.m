function X = fgnDaviesHarte_v001(N, H, nRep, opts)
%FGNDAVIESHARTE_V001 Exact fractional Gaussian noise by circulant embedding.
%   X = fgnDaviesHarte_v001(N, H, nRep) returns an N-by-nRep matrix of unit-
%   variance fGn with Hurst exponent H in (0,1), spectral exponent
%   alpha = 2H - 1 at low frequency. Exact (Davies & Harte 1987; Dietrich &
%   Newsam 1997): the covariance of every column equals fgnAutocov_v001
%   exactly, PROVIDED the circulant eigenvalues are non-negative. That is
%   checked on every call and errors loudly if violated (never clipped).
%   Columns come in independent pairs (real/imag parts of one FFT).
arguments
    N (1,1) double {mustBeInteger, mustBeGreaterThan(N, 1)}
    H (1,1) double {mustBeGreaterThan(H, 0), mustBeLessThan(H, 1)}
    nRep (1,1) double {mustBeInteger, mustBePositive} = 1
    opts.Stream = RandStream.getGlobalStream
    opts.EigTol (1,1) double = 1e-10   % relative tolerance for round-off negatives
end
g   = fgnAutocov_v001(0:N-1, H);
c   = [g, g(end-1:-1:2)];                 % first row of the 2(N-1) circulant
M   = numel(c);
lam = real(fft(c));
if min(lam) < -opts.EigTol * max(lam)
    error('fgnDaviesHarte:negEig', ['Circulant embedding not non-negative definite ' ...
        '(min eig %.3g, H = %.4f, N = %d): exact simulation impossible here.'], min(lam), H, N);
end
lam = max(lam, 0);                        % only round-off-sized negatives reach here
nPair = ceil(nRep / 2);
Z  = randn(opts.Stream, M, nPair) + 1i * randn(opts.Stream, M, nPair);
W  = fft(sqrt(lam(:) / M) .* Z);
W  = W(1:N, :);
X  = zeros(N, 2 * nPair);
X(:, 1:2:end) = real(W);  X(:, 2:2:end) = imag(W);   % independent pair per FFT
X  = X(:, 1:nRep);
end
