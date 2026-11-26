function [alphaEst, kUsed, info] = estimateWhittleAlpha(x, fs, fLow, fHigh, options)
%ESTIMATEWHITTLEALPHA Prewhitened, band-limited local Whittle estimate of a spectral exponent.
%
% Whole-record multitaper spectra, and IRASA computed on them
% (estimateIRASA), are biased at steep alpha when the noise has no low-frequency break: on exact
% k-fold-integrated fGn they read -2 to +4 off at alpha 4-6.5, depending on fs and N. This estimator recovers
% the band-limited exponent within +/-0.11 on exact references and on noiseXu output (fs 133/240 Hz,
% N 1000/3000, alpha 2.5-6.5).
%
% METHOD:
%   1. Pre-estimate alpha with DFA of order DFAOrder (default 3; valid for alpha < 2*DFAOrder + 1, Kiyono 2015).
%   2. Difference x k times, k = floor((alphaDFA + 1)/2) clamped to [0, MaxK] (or the K you pass). This removes
%      the low-frequency power that would otherwise leak into the fitting band ("prewhitening").
%   3. Hann-tapered periodogram of the differenced series; divide by |1 - exp(-i*lambda)|^(2k) to recolour.
%      Why the taper: if step 1 picks k one too low, the differenced series is still non-stationary and an untapered
%      periodogram leaks low-frequency power into the band. On exact references with k forced one too low the median bias is
%      -0.43 to -0.48 without the taper and within 0.03 with it; the price is about 1.4 times the SD when k is right
%      (validation/compareHannTaper_v001.m). This is the textbook trade: the untapered estimate is consistent only for d < 1
%      (alpha < 2) and the cosine bell extends that to d < 1.5 at a variance factor of 35/18 (Velasco 1999, pp. 93-94, 97, 101);
%      tapering and low-frequency trimming of differenced data are both needed under smooth trends (Hurvich, Lang & Soulier 2005).
%   4. Local Whittle objective (Robinson 1995) for S(lambda) ~ G*lambda^(-alpha) over [fLow, fHigh] Hz,
%      minimised with fminbnd over AlphaBounds. A minimum on a bound is an error, not an estimate.
%   The target is the in-band slope of the true spectrum (what IRASA over the same band also targets).
%
% SYNTAX:
%   alphaEst = estimateWhittleAlpha(x, fs, fLow, fHigh)
%   [alphaEst, kUsed, info] = estimateWhittleAlpha(x, fs, fLow, fHigh, K=2)
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
%   PeakCheck   - Name-value. Warn about a narrow peak in the band (default true; see PEAKS).
%   PeakFalseAlarm - Name-value. Chance of a false warning per call on peak-free data (default 1e-4).
%
% OUTPUTS:
%   alphaEst - Spectral exponent (S ~ f^-alpha) over [fLow, fHigh]
%   kUsed    - Differencing order used
%   info     - struct: alphaDFA (pre-estimate), nBins (periodogram bins in the band), objective (at the minimum),
%              peakRatio, peakThreshold, peakHz (the most significant bin of the peak check; it warns when ratio > threshold)
%
% ERRORS (fail loud): non-finite or constant input, fewer than 10 bins in the band, series too short for DFA,
%   optimiser on a bound or not converged.
%
% PEAKS: a rhythmic peak inside the band biases this estimate badly, and the estimate alone does not show it. On exact references
%   with a 10 Hz sine added (amplitude a fraction of the series' SD, band 1-20 Hz), the median moved from 2.52 to 1.01 at 10% for
%   alpha 2.5, and from 4.52 to -0.64 at 1% for alpha 4.5; estimateIRASA with Prewhiten = k stayed within 0.18 of the peak-free
%   value (validation/checkWhittlePeak_v001.m). So each call checks for a narrow peak: a bin far above the mean of its neighbours
%   (up to 12 either side, leaving out +/-2) once the fitted slope is flattened. The threshold comes from theory, not from data:
%   with no peak the ratio is F-distributed, and the threshold sets the chance of any false warning in the band to PeakFalseAlarm.
%   A hit warns "estimateWhittleAlpha:peakInBand" (the estimate is still returned); use estimateIRASA with Prewhiten=k, or raise K
%   if it is a harmonic of the movement left in the residual. Calibration (validation/calibratePeakCheck_v002.m): no warning in
%   16,800 peak-free series (exact references, noiseXu with and without its break, white noise); peaks that bias steep series are
%   caught (alpha 4.5 and 6.5: every 10 and 18 Hz sine from 0.3% of the SD, and from 0.1% at fs 240), but shallow series and low-
%   frequency peaks are caught only when larger, and some missed peaks still bias the estimate by up to 0.85 (alpha 2.5, 18 Hz,
%   1%; fs 133, N 1000). The check cannot see a broad bump.
%
% SAMPLING ERROR: for an untapered local Whittle estimate from m frequencies near zero, the estimate of H has asymptotic variance
%   1/(4m) (Robinson 1995, Theorem 2), so SD(alpha) = 1/sqrt(m), because alpha = 2H - 1. Tapering inflates the variance (Hurvich & Chen 2000
%   obtain 1.5/(4m) for their own complex-valued taper). For a full cosine-bell (Hann) taper on the undifferenced series, using all Fourier
%   frequencies from the second upward, Velasco (1999, p. 101, remark after Theorem 6 and Equation 10) gives a variance factor of 35/18 (1.94),
%   that is (35/18)/(4m), for d < 3/2 and no drift. That remark does not cover this estimator: here the taper follows k differences and the
%   periodogram is recoloured, k comes from the data, and the band is fixed in Hz and starts above the lowest frequencies. Measured on exact
%   references (fs 120-300 Hz, N 1000 and 3000, band 1-20 Hz), the variance of alpha is 1.92 times 1/sum(nu.^2), nu = log(lambda) - mean(log(lambda))
%   (median of 32 conditions, range 1.73 to 2.07), and the SD is about twice 1/sqrt(m). The factor agrees with 35/18, but it is a measurement,
%   not a result of Velasco's theorem; see validation/reportToleranceTable_v003.py.
%
% Requires base MATLAB only. References: Hurvich, C. M., & Chen, W. W. (2000). An efficient taper for potentially overdifferenced long-memory
% time series. Journal of Time Series Analysis, 21(2), 155-180. Robinson, P. M. (1995). Gaussian semiparametric estimation of long range
% dependence. Annals of Statistics, 23(5), 1630-1661. Velasco, C. (1999). Gaussian semiparametric estimation of non-stationary time series.
% Journal of Time Series Analysis, 20(1), 87-127. Hurvich, C., Lang, G., & Soulier, P. (2005). Estimation of long memory in the presence of
% a smooth nonparametric trend. Journal of the American Statistical Association, 100(471), 853-871. Kiyono, K. (2015). Establishing a direct connection between
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
    options.PeakCheck (1,1) logical = true
    options.PeakFalseAlarm (1,1) double {mustBePositive, mustBeLessThan(options.PeakFalseAlarm, 1)} = 1e-4
end
x = double(x(:));
if ~all(isfinite(x)), error('estimateWhittleAlpha:nonFinite', 'x contains NaN or Inf.'); end
if std(x) == 0, error('estimateWhittleAlpha:constant', 'x is constant (e.g. an empty residual); no exponent exists.'); end
if fHigh <= fLow || fHigh >= fs/2, error('estimateWhittleAlpha:band', 'Need fLow < fHigh < fs/2 (got [%g %g] Hz, fs %g).', fLow, fHigh, fs); end
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
    error('estimateWhittleAlpha:band', 'Only %d periodogram bins in [%g %g] Hz (N = %d, fs = %g).', info.nBins, fLow, fHigh, N, fs);
end
lam = 2*pi*f(m)/fs;  Ib = I(m) ./ (2 - 2*cos(lam)).^kUsed;
obj = @(a) log(mean(lam.^a .* Ib)) - a * mean(log(lam));
lo = options.AlphaBounds(1);  hi = options.AlphaBounds(2);
[alphaEst, info.objective, flag] = fminbnd(obj, lo, hi, optimset('TolX', 1e-6));
% Peak check (a warning; the estimate is unchanged). Flatten the fitted slope, then compare each bin with the mean of up to 12
% bins either side, leaving out the bin and its +/-2 neighbours (a Hann-tapered sine spreads over 3 bins). With no peak the bins
% are exponential, so the ratio is F(2, 2n) with tail (1 + x/n)^(-n); n is the number of neighbours shrunk by 18/35, the Hann
% correlation factor (Velasco 1999, p. 101). The threshold makes the chance of any false alarm in the band PeakFalseAlarm.
[info.peakRatio, info.peakThreshold, info.peakHz] = peakOverNeighbours(lam.^alphaEst .* Ib, f(m), options.PeakFalseAlarm);
if options.PeakCheck && info.peakRatio > info.peakThreshold
    warning('estimateWhittleAlpha:peakInBand', ['Narrow peak in the band at %.3g Hz: %.1f times its neighbours (threshold %.1f, ' ...
        'false-alarm rate %g per call). A peak biases this exponent. Use estimateIRASA with Prewhiten=k, or raise K if it is a ' ...
        'harmonic of the movement left in the residual (see harmonicContent).'], info.peakHz, info.peakRatio, info.peakThreshold, options.PeakFalseAlarm);
end
if flag ~= 1 || alphaEst - lo < 1e-3 || hi - alphaEst < 1e-3
    error('estimateWhittleAlpha:bound', 'Optimiser on its bound or not converged (alpha = %.4f, bounds [%g %g], flag %d).', alphaEst, lo, hi, flag);
end
end

function [r, thr, fPk] = peakOverNeighbours(z, f, falseAlarm)
% Most significant bin over its neighbours' mean (see the peak check above). Returns its ratio, its threshold and its frequency;
% NaN when no bin has at least 8 neighbours.
nb = numel(z);  W = 12;  gap = 2;  ratio = NaN(nb, 1);  nUsed = zeros(nb, 1);
for j = 1:nb
    idx = [max(1, j - W):j - gap - 1, j + gap + 1:min(nb, j + W)];
    if numel(idx) >= 8, ratio(j) = z(j) / mean(z(idx));  nUsed(j) = numel(idx); end
end
tested = ~isnan(ratio);
if ~any(tested), r = NaN;  thr = NaN;  fPk = NaN;  return, end
nEff = nUsed * 18 / 35;  pBin = falseAlarm / nnz(tested);
x = nEff .* (pBin .^ (-1 ./ nEff) - 1);              % P(F(2, 2 nEff) > x) = pBin
[~, k] = max(ratio ./ x);                              % the most significant bin
r = ratio(k);  thr = x(k);  fPk = f(k);
end

function a = dfaAlphaLocal(x, order)
% DFA of the given order as a spectral exponent alpha = 2*Hfluct - 1:
% 15 log-spaced box sizes 10..N/10, non-overlapping boxes, OLS of log F on log s.
y = cumsum(x - mean(x));  n = numel(y);
if n < 100, error('estimateWhittleAlpha:short', 'Series too short for the DFA pre-estimate (N = %d, need >= 100).', n); end
s = unique(round(exp(linspace(log(10), log(floor(n/10)), 15))));
F = zeros(size(s));
for i = 1:numel(s)
    wl = s(i);  nb = floor(n / wl);  Y = reshape(y(1:nb*wl), wl, nb);
    t = (1:wl)' / wl;  [Q, ~] = qr(t .^ (0:order), 0);
    R = Y - Q * (Q' * Y);  F(i) = sqrt(mean(R(:).^2));
end
p = polyfit(log(s), log(F), 1);  a = 2 * p(1) - 1;
end
