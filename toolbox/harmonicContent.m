function [K, info] = harmonicContent(x, fs, f0, options)
%HARMONICCONTENT Diagnostic: how many harmonics of f0 stand above the noise between them?
%
%   [K, info] = harmonicContent(x, fs, f0)
%   [K, info] = harmonicContent(x, fs, [], F0Band=[fLow fHigh])      % f0 from the record's spectral peak
%   [K, info] = harmonicContent(x, fs, f0, Cut=2, Kmax=20)
%
%   For each harmonic k = 1, 2, ... of f0 it fits one cosine and sine at k*f0 (with a constant and a linear trend) and
%   at the half-multiples (k - 0.5)*f0 and (k + 0.5)*f0, which fall between harmonics. The harmonic's amplitude divided by
%   the geometric mean of the two between-harmonic amplitudes is its ratio: how far it stands above the noise next to it.
%   The content count K is the highest k whose ratio is at least Cut, counting up from k = 1 and stopping at the first
%   harmonic below Cut; it is never less than 1.
%
%   THIS IS A DIAGNOSTIC, NOT A SELECTOR OF RECORD. Use info.ratio to see where the harmonics stop standing out, and report a
%   ladder of K (harmonicTemplate) if anything built on the residual moves with K. A harmonic that is not above the
%   noise is noise, and removing it carves a notch in the noise; a harmonic that is in the movement and is left behind stays
%   in the "noise". The rule is the content rule of the PowerLawSimulationPreReg v028 amendment 01; an information-criterion rule
%   was tried there and failed its acceptance test, so it is not included. Selection by information criterion and a stacked-spectrum
%   plot are planned for v2.1.0.
%
% INPUTS:
%   x    - Real N-by-C matrix (one column per channel) or a vector; finite. Amplitudes are averaged over channels.
%   fs   - Sampling rate in cycles per unit time (any consistent unit).
%   f0   - Fundamental frequency (> 0), or [] to estimate it (F0Band required).
%   F0Band - Name-value. [fLow fHigh] for the estimate, as in harmonicTemplate.
%   Cut    - Name-value. Ratio threshold (default 2).
%   Kmax   - Name-value. Highest harmonic examined (default 20; reduced so that (k + 0.5)*f0 < fs/2).
%
% OUTPUTS:
%   K    - Content count (integer, at least 1).
%   info - Struct: f0, ratio (1-by-Kmax, NaN where not computable), amp, ampNull (1-by-Kmax), Cut, Kmax, firstBelow (first harmonic
%          below Cut, or NaN if none), noHarmonic (true if even harmonic 1 is below Cut: K is then 1 by the floor, not by evidence).
%
% ERRORS (fail loud): non-finite or non-real x, f0 [] without F0Band, a record shorter than 3 cycles of f0, f0 too high for any harmonic below fs/2.
%
% EXAMPLE:
%   t = (0:1/100:30)';  x = 2*sin(2*pi*1*t) + 0.5*sin(2*pi*2*t+0.7) + cumsum(randn(size(t)))*0.01;
%   [K, info] = harmonicContent(x, 100, 1);    % K is 2
%
% See also: harmonicTemplate, shapeNoiseWrapper

arguments
    x double {mustBeReal, mustBeNonempty}
    fs (1, 1) double {mustBePositive}
    f0 double {mustBeScalarOrEmpty}
    options.F0Band (1, 2) double = [NaN NaN]
    options.Cut (1, 1) double {mustBePositive} = 2
    options.Kmax (1, 1) double {mustBeInteger, mustBePositive} = 20
end
if isvector(x), x = x(:); end
[N, C] = size(x);
if any(~isfinite(x), "all"), error("harmonicContent:nonFinite", "x contains NaN or Inf."); end
if isempty(f0)
    if any(isnan(options.F0Band)) || options.F0Band(1) <= 0 || options.F0Band(2) <= options.F0Band(1) || options.F0Band(2) >= fs/2
        error("harmonicContent:f0Band", "f0 is [] so F0Band = [fLow fHigh] (0 < fLow < fHigh < fs/2) is required: there is no default band.");
    end
    f0 = estimateFundamental(x, fs, options.F0Band);
end
if ~(f0 > 0), error("harmonicContent:f0", "f0 must be positive."); end
if N / fs * f0 < 3, error("harmonicContent:tooShort", "The record spans %.2g cycles of f0; at least 3 are needed.", N / fs * f0); end
Kmax = min(options.Kmax, floor(fs / (2 * f0) - 0.5 - 1e-9));
if Kmax < 1, error("harmonicContent:nyquist", "f0 = %g leaves no harmonic with (k + 0.5)*f0 below fs/2 = %g.", f0, fs / 2); end

t = (0:N - 1)' / fs;  D0 = [ones(N, 1), t];
% each fit below includes the constant and linear trend, so no separate detrending is needed
nf = 2 * Kmax + 1;  A = zeros(C, nf);
for j = 1:nf
    fr = 0.5 * j * f0;
    D = [D0, cos(2 * pi * fr * t), sin(2 * pi * fr * t)];
    b = D \ x;
    A(:, j) = hypot(b(3, :), b(4, :))';
end
Am = mean(A, 1);
amp = NaN(1, Kmax);  ampNull = NaN(1, Kmax);
for k = 1:Kmax
    amp(k) = Am(2 * k);  ampNull(k) = sqrt(Am(2 * k - 1) * Am(2 * k + 1));
end
ratio = amp ./ ampNull;
below = find(~(ratio >= options.Cut), 1);
if isempty(below), K = Kmax;  fb = NaN; else, K = max(1, below - 1);  fb = below; end
info = struct("f0", f0, "ratio", ratio, "amp", amp, "ampNull", ampNull, "Cut", options.Cut, "Kmax", Kmax, ...
    "firstBelow", fb, "noHarmonic", ~(ratio(1) >= options.Cut));
end

function f0 = estimateFundamental(x, fs, band)
% highest peak of the detrended record's zero-padded FFT magnitude in the band, per channel (as harmonicTemplate)
[N, C] = size(x);  nfft = 2^nextpow2(4 * N);  fa = (0:nfft - 1)' * fs / nfft;
idx = find(fa >= band(1) & fa <= band(2));
if isempty(idx), error("harmonicContent:f0Band", "F0Band [%g %g] contains no FFT bin of this record.", band(1), band(2)); end
pk = zeros(1, C);  fp = zeros(1, C);
for c = 1:C
    s = x(:, c) - mean(x(:, c));  t = (0:N - 1)' / fs;  s = s - [ones(N, 1), t] * ([ones(N, 1), t] \ s);
    m = abs(fft(s, nfft));  [pk(c), j] = max(m(idx));  fp(c) = fa(idx(j));
end
if C == 1 || all(abs(fp - fp(1)) / max(fp(1), eps) < 0.2), f0 = mean(fp);
else, [~, c] = max(pk);  f0 = fp(c); end
end
