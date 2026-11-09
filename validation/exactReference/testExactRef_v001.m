function results = testExactRef_v001()
%TESTEXACTREF_V001 Self-tests for the V0 exact reference set. Errors on any failure.
%   T1 circulant eigenvalues >= 0 over H in [0.01, 0.99], N up to 16384
%   T2 analytic fGn spectrum integrates to the unit variance
%   T3 Monte Carlo: sample autocovariance of Davies-Harte output matches theory (lags 0-10)
%   T4 Davies-Harte vs Cholesky (N = 512): same covariance target, same sample statistics
%   T5 odd-integer alpha errors loudly
%   T6 alpha -> (H, k) mapping round-trips over the plan grid
rs = RandStream('threefry', Seed=20261107);
fail = strings(0,1);  results = struct();

% T1
Hs = [0.01 0.05:0.05:0.95 0.99];  Ns = [256 1024 4096 16384];  worst = inf;
for H = Hs
    for N = Ns
        g = fgnAutocov_v001(0:N-1, H);  lam = real(fft([g, g(end-1:-1:2)]));
        worst = min(worst, min(lam) / max(lam));
    end
end
results.T1_minRelEig = worst;
if worst < -1e-10; fail(end+1) = sprintf("T1 negative circulant eigenvalue (rel %.3g)", worst); end

% T2
[~, info] = exactRefSeries_v001(64, 0.4, 1, Stream=rs);       % H = 0.7, k = 0
v = 2 * integral(@(l) info.S(l(:)).', 1e-9, pi);
results.T2_specIntegral = v;
if abs(v - 1) > 1e-3; fail(end+1) = sprintf("T2 spectrum integrates to %.5f, not 1", v); end

% T3
nRep = 20000;  lags = 0:10;  Htest = [0.2 0.85];
for i = 1:2
    X = fgnDaviesHarte_v001(128, Htest(i), nRep, Stream=rs);
    est = zeros(nRep, numel(lags));                             % per-replicate autocov (known zero mean)
    for L = lags; est(:, L+1) = sum(X(1:end-L,:) .* X(1+L:end,:), 1)' / (128 - L); end
    z = (mean(est) - fgnAutocov_v001(lags, Htest(i))) ./ (std(est) / sqrt(nRep));   % empirical MC SE
    results.(sprintf('T3_maxAbsZ_H%02d', round(100*Htest(i)))) = max(abs(z));
    if max(abs(z)) > 4.5; fail(end+1) = sprintf("T3 H=%.2f autocov |z| = %.2f (> 4.5)", Htest(i), max(abs(z))); end %#ok<AGROW>
end

% T4
N = 512;  H = 0.75;  nRep = 4000;
C = toeplitz(fgnAutocov_v001(0:N-1, H));  Lc = chol(C, 'lower');
Xc = Lc * randn(rs, N, nRep);  Xd = fgnDaviesHarte_v001(N, H, nRep, Stream=rs);
vc = var(mean(Xc, 1));  vd = var(mean(Xd, 1));  vt = sum(C(:)) / N^2;   % variance of the sample mean
results.T4_varMean = [vt vc vd];
if any(abs([vc vd] / vt - 1) > 0.1); fail(end+1) = sprintf("T4 var(mean): theory %.4g, chol %.4g, DH %.4g", vt, vc, vd); end

% T5
try
    exactRefSeries_v001(64, 3, 1); fail(end+1) = "T5 alpha=3 did not error";
catch e
    if ~strcmp(e.identifier, 'exactRefSeries:oddAlpha'); fail(end+1) = "T5 wrong error: " + e.identifier; end
end

% T6
grid = setdiff(-0.75:0.25:6.75, [1 3 5]);
for a = grid
    [~, I] = exactRefSeries_v001(32, a, 1, Stream=rs);
    if abs((2*I.H - 1) + 2*I.k - a) > 1e-12 || I.H <= 0 || I.H >= 1; fail(end+1) = sprintf("T6 alpha %.2f maps to H %.3f k %d", a, I.H, I.k); end %#ok<AGROW>
end

results.failures = fail;
disp(results)
if ~isempty(fail); error('testExactRef:fail', 'V0 self-tests FAILED (%d):\n  - %s', numel(fail), strjoin(fail, newline + "  - ")); end
fprintf('V0 self-tests: all pass.\n');
end
