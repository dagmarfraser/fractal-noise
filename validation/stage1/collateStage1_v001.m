function S = collateStage1_v001(out)
%COLLATESTAGE1_V001 Collate stage-1 long-format estimates (any subset of languages) against the manifest.
%   S = collateStage1_v001(outDir) writes estimates_all.csv, summary.csv (per estimator x cell: failures,
%   median, bias/RMSE vs nominal alpha AND vs band-effective alpha_eff, SD, s/series) and kselect.csv
%   (share of series where the data-chosen k equals the generator k). FAIL rows are counted, never dropped.
M = readtable(fullfile(out, "manifest.csv"), TextType="string");
f = dir(fullfile(out, "estimates_*.csv"));  f = f(~strcmp({f.name}, 'estimates_all.csv'));
if isempty(f); error('collate:none', 'No estimates_*.csv in %s', out); end
vn = ["cell_id","rep","estimator","language","alpha_hat","seconds","status","k_used"];
A = table();
for i = 1:numel(f)
    T = readtable(fullfile(out, f(i).name), TextType="string", TreatAsMissing=["NA","nan","NaN"], Delimiter=",");
    miss = setdiff(vn, T.Properties.VariableNames);
    if ~isempty(miss); error('collate:cols', '%s lacks column(s): %s', f(i).name, strjoin(miss, ", ")); end
    for v = ["alpha_hat","seconds","k_used"]; if ~isnumeric(T.(v)); T.(v) = double(T.(v)); end; end
    A = [A; T(:, vn)]; %#ok<AGROW>
end
writetable(A, fullfile(out, "estimates_all.csv"));
A = join(A, M(:, ["cell_id","alpha","alpha_eff","N","k"]));
isK = A.estimator == "kselect_whittle";
K = A(isK, :);  A = A(~isK, :);
ok = A.status == "ok";
[g, est, cid] = findgroups(A.estimator, A.cell_id);
first = @(v) splitapply(@(a) a(1), v, g);
S = table(est, cid, first(A.alpha), first(A.alpha_eff), first(A.N), splitapply(@sum, ok, g), splitapply(@sum, ~ok, g), ...
    VariableNames=["estimator","cell_id","alpha","alpha_eff","N","nOK","nFail"]);
S.median  = splitapply(@(a, o) median(a(o)), A.alpha_hat, ok, g);
S.bias    = S.median - S.alpha;
S.biasEff = S.median - S.alpha_eff;
S.sd      = splitapply(@(a, o) std(a(o)), A.alpha_hat, ok, g);
S.rmse    = splitapply(@(a, t, o) sqrt(mean((a(o) - t(o)).^2)), A.alpha_hat, A.alpha, ok, g);
S.rmseEff = splitapply(@(a, t, o) sqrt(mean((a(o) - t(o)).^2)), A.alpha_hat, A.alpha_eff, ok, g);
S.secPerSeries = splitapply(@(t) median(t, "omitnan"), A.seconds, g);
S = sortrows(S, ["N","alpha","estimator"]);
writetable(S, fullfile(out, "summary.csv"));
if ~isempty(K)
    [gk, kc] = findgroups(K.cell_id);
    KS = table(kc, splitapply(@(a) a(1), K.alpha, gk), splitapply(@(a) a(1), K.N, gk), splitapply(@(a) a(1), K.k, gk), ...
        splitapply(@(u, k) mean(u == k), K.k_used, K.k, gk), splitapply(@(u, k) mean(u < k), K.k_used, K.k, gk), ...
        splitapply(@(u, k) mean(u > k), K.k_used, K.k, gk), splitapply(@(s) sum(s ~= "ok"), K.status, gk), ...
        VariableNames=["cell_id","alpha","N","k","pCorrect","pUnder","pOver","nFail"]);
    writetable(sortrows(KS, ["N","alpha"]), fullfile(out, "kselect.csv"));
end
fprintf('Collated %d estimate rows (%d FAIL) from %d file(s) -> %s\n', height(A) + height(K), sum(~ok) + sum(K.status ~= "ok"), numel(f), out);
end
