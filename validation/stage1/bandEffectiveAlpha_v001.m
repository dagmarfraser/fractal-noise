function aEff = bandEffectiveAlpha_v001(S, fs, band, N)
%BANDEFFECTIVEALPHA_V001 In-band log-log slope of an analytic spectrum (the target a slope estimator can hit).
%   aEff = bandEffectiveAlpha_v001(S, fs, band, N): S is info.S from exactRefSeries_v001 (lambda in rad/sample).
%   Evaluated on the Welch grid the stage-1 runners use (nperseg = min(1024, N/2)), OLS of log10 S on log10 f.
%   Needed because discrete (anti)differencing is not a pure f^-2 per order: the gain (2 sin(pi f/fs))^2 adds a
%   local slope 2x cot(x), x = pi f/fs (1.73 at f = fs/5), and the fGn spectrum itself flattens near Nyquist.
nps = min(1024, floor(N/2));
f = (1:floor(nps/2)) * fs / nps;  f = f(f >= band(1) & f <= band(2));
if numel(f) < 3; error('bandEffectiveAlpha:band', 'Only %d Welch bins inside [%g %g] Hz at N = %d.', numel(f), band, N); end
p = polyfit(log10(f), log10(S(2*pi*f(:)/fs)).', 1);
aEff = -p(1);
end
