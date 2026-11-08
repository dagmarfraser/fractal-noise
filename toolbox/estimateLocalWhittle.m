function [alphaEst, kUsed, info] = estimateLocalWhittle(x, fs, fLow, fHigh, options)
%ESTIMATELOCALWHITTLE Prewhitened, band-limited local Whittle estimate of a spectral exponent.
%
% Toolbox wrap of PowerLawSimulationPreReg's alphaPWLW_v001.m and dfaAlpha_v001.m (2026-11-08); the numerics
% are unchanged from source. Built because whole-record multitaper spectra, and IRASA computed on them
% (estimateIRASA), are biased at steep alpha when the noise has no low-frequency break: on exact
% k-fold-integrated fGn they read -2 to +4 off at alpha 4-6.5, depending on fs and N. This estimator recovers
% the band-limited exponent within +/-0.11 on exact references and on noiseXu output (fs 133/240 Hz,
% N 1000/3000, alpha 2.5-6.5). See the fractal-noise validation plan v008, findings F2-F4.
%
% METHOD:
%   1. Pre-estimate alpha with DFA of order DFAOrder (default 3; valid for alpha < 2*DFAOrder + 1, Kiyono 2015).
%   2. Difference x k times, k = floor((alphaDFA + 1)/2) clamped to [0, MaxK] (or the K you pass). This removes
%      the low-frequency power that would otherwise leak into the fitting band ("prewhitening").
%   3. Hann-tapered periodogram of the differenced series; divide by |1 - exp(-i*lambda)|^(2k) to recolour.
%   4. Local Whittle objective (Robinson 1995) for S(lambda) ~ G*lambda^(-alpha) over [fLow, fHigh] Hz,
%      minimised with fminbnd over AlphaBounds. A minimum on a bound is an error, not an estimate.
%   The target is the in-band slope of the true spectrum (what IRASA over the same band also targets).
%
% SYNTAX:
%   alphaEst = estimateLocalWhittle(x, fs, fLow, fHigh)
%   [alphaEst, kUsed, info] = estimateLocalWhittle(x, fs, fLow, fHigh, K=2)
%
% INPUTS:
%   x      - Time series, real vector (finite, non-constant)
%   fs     - Sampling rate in Hz
%   fLow   - Lower edge of the fitting band, Hz
%   fHigh  - Upper edge of the fitting band, Hz (must be below fs/2)
%   K           - Name-value. Differencing order; [] (default) chooses it from the DFA pre-estimate.
%   DFAOrder    - Name-value. Order of the DFA pre-estimate, default 3.
%   MaxK        - Name-value. Largest automatic k, default 3.
%   AlphaBounds - Name-value. Search interval for alpha, default [-2 10].
%
% OUTPUTS:
%   alphaEst - Spectral exponent (S ~ f^-alpha) over [fLow, fHigh]
%   kUsed    - Differencing order used
%   info     - struct: alphaDFA (pre-estimate), nBins (periodogram bins in the band), objective (at the minimum)
%
% ERRORS (fail loud): non-finite or constant input, fewer than 10 bins in the band, series too short for DFA,
%   optimiser on a bound or not converged.
%
% Requires base MATLAB only. Reference: Robinson, P. M. (1995). Gaussian semiparametric estimation of long range
% dependence. Annals of Statistics, 23(5), 1630-1661. Kiyono, K. (2015). Establishing a direct connection between
% detrended fluctuation analysis and Fourier analysis. Physical Review E, 92, 042925.

arguments
    x {mustBeNumeric, mustBeReal, mustBeVector}
    fs (1,1) double {mustBePositive}
    fLow (1,1) double {mustBePositive}
    fHigh (1,1) double {mustBePositive}
    options.K double {mustBeScalarOrEmpty, mustBeNonnegative, mustBeInteger} = []
    options.DFAOrder (1,1) double {mustBeInteger, mustBePositive} = 3
    options.MaxK (1,1) double {mustBeInteger, mustBeNonnegative} = 3
    options.AlphaBounds (1,2) double = [-2 10]
end
x = double(x(:));
if ~all(isfinite(x)), error('estimateLocalWhittle:nonFinite', 'x contains NaN or Inf.'); end
if std(x) == 0, error('estimateLocalWhittle:constant', 'x is constant (e.g. an empty residual); no exponent exists.'); end
if fHigh <= fLow || fHigh >= fs/2, error('estimateLocalWhittle:band', 'Need fLow < fHigh < fs/2 (got [%g %g] Hz, fs %g).', fLow, fHigh, fs); end
x = x - mean(x);
info.alphaDFA = dfaAlphaLocal(x, options.DFAOrder);
if isempty(options.K)
    kUsed = min(max(floor((info.alphaDFA + 1) / 2), 0), options.MaxK);
else
    kUsed = options.K;
end
xd = x;  for j = 1:kUsed, xd = diff(xd); end
N = numel(xd);  w = 0.5 * (1 - cos(2*pi*(0:N-1)' / N));
I = abs(fft((xd - mean(xd)) .* w)).^2 / sum(w.^2);
f = (0:N-1)' * fs / N;  m = f >= fLow & f <= fHigh & f < fs/2;
info.nBins = nnz(m);
if info.nBins < 10
    error('estimateLocalWhittle:band', 'Only %d periodogram bins in [%g %g] Hz (N = %d, fs = %g).', info.nBins, fLow, fHigh, N, fs);
end
lam = 2*pi*f(m)/fs;  Ib = I(m) ./ (2 - 2*cos(lam)).^kUsed;
obj = @(a) log(mean(lam.^a .* Ib)) - a * mean(log(lam));
lo = options.AlphaBounds(1);  hi = options.AlphaBounds(2);
[alphaEst, info.objective, flag] = fminbnd(obj, lo, hi, optimset('TolX', 1e-6));
if flag ~= 1 || alphaEst - lo < 1e-3 || hi - alphaEst < 1e-3
    error('estimateLocalWhittle:bound', 'Optimiser on its bound or not converged (alpha = %.4f, bounds [%g %g], flag %d).', alphaEst, lo, hi, flag);
end
end

function a = dfaAlphaLocal(x, order)
% DFA of the given order as a spectral exponent alpha = 2*Hfluct - 1 (copy of dfaAlpha_v001; MATLAB twin of the
% fractal-noise stage-1 R dfaN): 15 log-spaced box sizes 10..N/10, non-overlapping boxes, OLS of log F on log s.
y = cumsum(x - mean(x));  n = numel(y);
if n < 100, error('estimateLocalWhittle:short', 'Series too short for the DFA pre-estimate (N = %d, need >= 100).', n); end
s = unique(round(exp(linspace(log(10), log(floor(n/10)), 15))));
F = zeros(size(s));
for i = 1:numel(s)
    wl = s(i);  nb = floor(n / wl);  Y = reshape(y(1:nb*wl), wl, nb);
    t = (1:wl)' / wl;  [Q, ~] = qr(t .^ (0:order), 0);
    R = Y - Q * (Q' * Y);  F(i) = sqrt(mean(R(:).^2));
end
p = polyfit(log(s), log(F), 1);  a = 2 * p(1) - 1;
end
