function [resid, fit, info] = harmonicTemplate(x, fs, f0, K, options)
%HARMONICTEMPLATE Windowed harmonic-regression template: strip a periodic component, keep the residual.
%
%   [resid, fit, info] = harmonicTemplate(x, fs, f0, K)
%   [resid, fit, info] = harmonicTemplate(x, fs, [], K, F0Band=[fLow fHigh])   % f0 from the record's spectral peak
%   [___] = harmonicTemplate(___, Window=4, EndAnchored=true, Trend=true, Crop=true)
%
%   Fits, in overlapping Hann-weighted windows of Window cycles of f0 (hop half a window), a constant, an optional linear
%   trend and K harmonics of f0 (a cosine and sine at each of f0, 2 f0, ..., K f0) by weighted least squares, subtracts the
%   fit and overlap-adds the pieces. With EndAnchored an extra final window ends on the last sample, so every sample is
%   covered. The residual is what shapeNoise takes: pass it (Crop=true) straight in.
%
%   THIS IS A TEMPLATE, NOT A SELECTOR. K and Window are yours to choose. K = 1 removes only the fundamental. A harmonic whose
%   amplitude is no larger than the noise next to it is noise, and removing it carves a notch in the noise and flattens its
%   fitted exponent; a harmonic that really is in the movement and is left behind stays in the "noise". The count K and the window
%   are the decisions that matter, so report them, and report a ladder of K if the exponent of the residual or anything
%   built on it moves with K. Selection rules (a content rule and an information-criterion rule) and a stacked-spectrum
%   diagnostic are planned for v2.1.0.
%
% INPUTS:
%   x    - Real N-by-C matrix (one column per channel) or a vector; finite.
%   fs   - Sampling rate, in cycles per unit of time (any consistent unit; f0 and F0Band use the same unit).
%   f0   - Fundamental frequency (> 0, below fs/(2 K)), one value for all channels, or [] to estimate it.
%   K    - Number of harmonics (positive integer; K f0 must be below fs/2).
%   F0Band - Name-value. [fLow fHigh] in the units of fs for the estimate: the highest peak of the detrended record's
%          zero-padded FFT magnitude, per channel; if the channel estimates agree within 20% their mean is used, else the
%          estimate of the channel with the larger peak. REQUIRED when f0 is []: there is no default band.
%          Keep the band narrow around the expected tempo: on steep noise the highest peak can sit at the band's low edge, and the
%          estimate is then wrong without any error.
%   Window - Name-value. Window length in cycles of f0 (default 4, at least 2). The window is at least 2 cycles and at most N samples.
%   EndAnchored - Name-value. Add a final window that ends on the last sample (default true).
%   Trend  - Name-value. Include a linear trend in each window (default true).
%   Crop   - Name-value. Drop floor(window/4) samples at each end (default true), where the overlap-add weights are small
%          and the residual is least reliable; info.index gives the kept samples. false returns the full length.
%
% OUTPUTS:
%   resid - The residual: (N - 2 floor(window/4))-by-C if Crop, else N-by-C.
%   fit   - The fitted template over the same samples as resid.
%   info  - Struct: f0, K, WindowSamples, Hop, nWindows, index = [first last] of the returned samples in x, Covered
%           (true when every returned sample received a window), Crop.
%
% ERRORS (fail loud): non-finite or non-real x, a record shorter than one window ("harmonicTemplate:tooShort": no window fits,
%   so no template exists; the function does not return zeros), f0 [] without F0Band, K f0 at or above fs/2.
%
% EXAMPLE:
%   t = (0:1/100:30)';  x = 3*sin(2*pi*0.5*t) + cumsum(randn(size(t)))*0.05;
%   [r, ~, info] = harmonicTemplate(x, 100, 0.5, 1);   % the 0.5 Hz fundamental removed
%
% Reference for the construction: the windowed harmonic regression of the PowerLawSimulationPreReg runner
% (runLoopClosureFftnoise_v028, templateSubtract_local), generalised from two movement axes to any channels, with the
% record-too-short case an error instead of a zero residual. Related methods: dynamic harmonic regression (Young, Pedregal
% & Tych 1999, J Forecasting 18, 369-394), and time-varying harmonic removal in overlap-add frames (Zivanovic et al. 2023,
% Mech Syst Signal Process 189, 110043).
%
% See also: shapeNoise, estimateWhittleAlpha

arguments
    x double {mustBeReal, mustBeNonempty}
    fs (1, 1) double {mustBePositive}
    f0 double {mustBeScalarOrEmpty}
    K (1, 1) double {mustBeInteger, mustBePositive}
    options.F0Band (1, 2) double = [NaN NaN]
    options.Window (1, 1) double {mustBeGreaterThanOrEqual(options.Window, 2)} = 4
    options.EndAnchored (1, 1) logical = true
    options.Trend (1, 1) logical = true
    options.Crop (1, 1) logical = true
end
if isvector(x), x = x(:); end
[N, C] = size(x);
if any(~isfinite(x), "all"), error("harmonicTemplate:nonFinite", "x contains NaN or Inf."); end
if isempty(f0)
    if any(isnan(options.F0Band)) || options.F0Band(1) <= 0 || options.F0Band(2) <= options.F0Band(1) || options.F0Band(2) >= fs/2
        error("harmonicTemplate:f0Band", "f0 is [] so F0Band = [fLow fHigh] (0 < fLow < fHigh < fs/2) is required: there is no default band.");
    end
    f0 = estimateFundamental(x, fs, options.F0Band);
end
if ~(f0 > 0), error("harmonicTemplate:f0", "f0 must be positive."); end
if K * f0 >= fs / 2, error("harmonicTemplate:nyquist", "K*f0 = %g is not below fs/2 = %g.", K * f0, fs / 2); end

win = round(options.Window / f0 * fs);  win = min(win, N);  win = max(win, round(2 / f0 * fs));
hop = round(win / 2);
if win > N || win < 2 * (2 * K + 2)
    error("harmonicTemplate:tooShort", ...
        "The record (%d samples) is shorter than one template window (%d samples = %g cycles of f0 = %g, at least 2 cycles): no template can be fitted.", ...
        N, max(win, round(2 / f0 * fs)), options.Window, f0);
end
starts = 1:hop:(N - win + 1);
if options.EndAnchored && starts(end) + win - 1 < N, starts(end + 1) = N - win + 1; end

res = zeros(N, C);  ww = zeros(N, 1);
tw = (0:win - 1)' / fs;  han = 0.5 * (1 - cos(2 * pi * (0:win - 1)' / (win - 1)));
D = ones(win, 1);  if options.Trend, D = [D, tw]; end
for h = 1:K, D = [D, cos(2 * pi * h * f0 * tw), sin(2 * pi * h * f0 * tw)]; end %#ok<AGROW>
Dw = D .* han;
for s = starts
    e = s + win - 1;
    b = Dw \ (x(s:e, :) .* han);
    res(s:e, :) = res(s:e, :) + (x(s:e, :) - D * b) .* han;
    ww(s:e) = ww(s:e) + han;
end
covered = ww > 0;
res(covered, :) = res(covered, :) ./ ww(covered);

if options.Crop
    cl = max(round(win / 4), 1);  first = cl;  last = min(max(N - cl, first + 100), N);
else
    first = 1;  last = N;
end
resid = res(first:last, :);  fit = x(first:last, :) - resid;
info = struct("f0", f0, "K", K, "WindowSamples", win, "Hop", hop, "nWindows", numel(starts), "index", [first last], ...
    "Covered", all(covered(first:last)), "Crop", options.Crop);
end

function f0 = estimateFundamental(x, fs, band)
% highest peak of the detrended record's zero-padded FFT magnitude in the band, per channel
[N, C] = size(x);  nfft = 2^nextpow2(4 * N);  fa = (0:nfft - 1)' * fs / nfft;
idx = find(fa >= band(1) & fa <= band(2));
if isempty(idx), error("harmonicTemplate:f0Band", "F0Band [%g %g] contains no FFT bin of this record.", band(1), band(2)); end
pk = zeros(1, C);  fp = zeros(1, C);
for c = 1:C
    s = x(:, c) - mean(x(:, c));  t = (0:N - 1)' / fs;  s = s - [ones(N, 1), t] * ([ones(N, 1), t] \ s);
    m = abs(fft(s, nfft));  [pk(c), j] = max(m(idx));  fp(c) = fa(idx(j));
end
if C == 1 || all(abs(fp - fp(1)) / max(fp(1), eps) < 0.2), f0 = mean(fp);
else, [~, c] = max(pk);  f0 = fp(c); end
end
