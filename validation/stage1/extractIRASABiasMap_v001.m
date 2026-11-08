function T = extractIRASABiasMap_v001(runDir, varargin)
% EXTRACTIRASABIASMAP_V001  Bias map of the MATLAB-arm estimators from a stage-1 run (knife item K-IRASA1).
%   T = extractIRASABiasMap_v001(runDir)  reads runDir/estimates_matlab.csv and runDir/manifest.csv and
%   returns one row per estimator x band x alpha x N: n, nFail, median, IQR, bias vs nominal alpha,
%   biasEff vs band-effective alpha (manifest alpha_eff_<band>), and dMed (step in median from the
%   previous alpha at the same N; <= 0 flags a non-monotone response). Writes biasMap_matlab_v001.csv
%   into runDir and prints the operating-point slice (Alphas, Ns).
%   CAVEAT: stage-1 runs at fs = 100 Hz; the pre-reg IRASA (iraAlphaSigma_v001) runs at 240 Hz (or
%   FS_DS). Band-effective alpha differs with fs, so this is indicative, not the pre-reg operating point.
p = inputParser;
addParameter(p, 'Alphas', [2.5 6.75]);
addParameter(p, 'Ns', [2048 4096]);
addParameter(p, 'Estimator', "matlab_HB_irasa");
parse(p, varargin{:}); o = p.Results;

E = readtable(fullfile(runDir, 'estimates_matlab.csv'), 'TextType', 'string', 'TreatAsMissing', {'NA','nan','NaN'});
M = readtable(fullfile(runDir, 'manifest.csv'), 'TextType', 'string');
[tf, loc] = ismember(E.cell_id, M.cell_id);
assert(all(tf), 'extractIRASABiasMap:manifest', '%d estimate rows have no manifest cell', nnz(~tf));
E.alpha = M.alpha(loc); E.N = M.N(loc);
E.alphaEff = NaN(height(E), 1);
for b = unique(E.band)'
    col = "alpha_eff_" + b;
    assert(ismember(col, M.Properties.VariableNames), 'extractIRASABiasMap:band', 'manifest lacks %s', col);
    r = E.band == b; E.alphaEff(r) = M.(col)(loc(r));
end
ok = E.status == "ok" & isfinite(E.alpha_hat);

[G, est, band, alpha, N] = findgroups(E.estimator, E.band, E.alpha, E.N);
n      = splitapply(@sum, ok, G);
nFail  = splitapply(@(x) sum(~x), ok, G);
med    = splitapply(@(a, k) median(a(k)), E.alpha_hat, ok, G);
iqrv   = splitapply(@(a, k) iqr(a(k)), E.alpha_hat, ok, G);
aEff   = splitapply(@(a) a(1), E.alphaEff, G);
T = table(est, band, alpha, N, n, nFail, med, iqrv, med - alpha, med - aEff, aEff, ...
    'VariableNames', {'estimator','band','alpha','N','n','nFail','median','IQR','bias','biasEff','alphaEff'});
T = sortrows(T, {'estimator','band','N','alpha'});
T.dMed = [NaN; diff(T.median)];
newRun = [true; T.N(2:end) ~= T.N(1:end-1) | T.band(2:end) ~= T.band(1:end-1) | T.estimator(2:end) ~= T.estimator(1:end-1)];
T.dMed(newRun) = NaN;
writetable(T, fullfile(runDir, 'biasMap_matlab_v001.csv'));

S = T(T.estimator == o.Estimator & T.alpha >= o.Alphas(1) & T.alpha <= o.Alphas(2) & ismember(T.N, o.Ns), :);
fprintf('%s, alpha %.2f-%.2f, N %s (fs 100 Hz: indicative only)\n', o.Estimator, o.Alphas, mat2str(o.Ns));
disp(S(:, {'band','N','alpha','n','nFail','median','IQR','bias','biasEff','dMed'}));
nm = S(S.dMed <= 0, :);
if ~isempty(nm), warning('extractIRASABiasMap:nonMonotone', '%d non-monotone steps in the slice', height(nm)); end
end
