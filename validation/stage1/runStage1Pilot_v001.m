function S = runStage1Pilot_v001(opts)
%RUNSTAGE1PILOT_V001 Stage-1 timing pilot: all estimators on V0 exact series, 3 languages.
%   Writes a hub-neutral manifest + series CSVs, runs the MATLAB estimators
%   (HB/FT/PM via estimateIRASA), then the Python and R runners on the same
%   files, collates long-format estimates and prints per-estimator bias,
%   SD, RMSE, failures and seconds/series, plus a projected full-grid budget.
%   Pilot only: band and differencing choices are provisional (plan open decisions).
arguments
    opts.Alphas (1,:) double = [0.5 4]
    opts.N (1,1) double = 4096
    opts.NRep (1,1) double = 10
    opts.MasterSeed (1,1) double = 20261108
    opts.Fs (1,1) double = 100
    opts.Band (1,2) double = [1 20]
    opts.PythonExe (1,1) string = "~/venvs/fnval/bin/python"
    opts.RscriptExe (1,1) string = "Rscript"
    opts.OutDir (1,1) string = ""
    opts.FieldTripDir (1,1) string = fullfile(getenv("HOME"), "Library/CloudStorage/Dropbox/MATLAB/fieldtrip-20251201")
    opts.Cores (1,1) double = 10
end
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "exactReference"));
if isempty(which('ft_defaults')); addpath(opts.FieldTripDir); ft_defaults; end
if isempty(which('estimateIRASA')); error('pilot:path', 'fractalnoise toolbox (estimateIRASA) not on path.'); end
out = opts.OutDir;  if out == ""; out = fullfile(here, "pilot_out_" + string(datetime("now", Format="yyyyMMdd_HHmmss"))); end
mkdir(out);  sdir = fullfile(out, "series");  mkdir(sdir);
py = regexprep(opts.PythonExe, "^~", getenv("HOME"));

% Manifest + series
ms = RandStream('threefry', Seed=opts.MasterSeed);  nC = numel(opts.Alphas);
M = table('Size', [nC 8], 'VariableTypes', ["string","double","double","double","double","double","double","string"], ...
    'VariableNames', ["cell_id","alpha","N","k","H","nrep","seed","series_file"]);
Xs = cell(nC, 1);
for c = 1:nC
    seed = randi(ms, 2^31 - 1);
    [X, info] = exactRefSeries_v001(opts.N, opts.Alphas(c), opts.NRep, Stream=RandStream('threefry', Seed=seed));
    id = sprintf("a%+.2f_N%d", opts.Alphas(c), opts.N);
    M(c, :) = {id, opts.Alphas(c), opts.N, info.k, info.H, opts.NRep, seed, id + ".csv"};
    writematrix(X, fullfile(sdir, id + ".csv"));  Xs{c} = X;
end
mf = fullfile(out, "manifest.csv");  writetable(M, mf);

% MATLAB arm
rows = {};
for c = 1:nC
    for r = 1:opts.NRep
        x = Xs{c}(:, r);
        t0 = tic;  aHB = estimateIRASA(x, opts.Fs, opts.Band(1), opts.Band(2));  tHB = toc(t0);
        t0 = tic;  [~, ~, ~, ~, ftI, pmI] = estimateIRASA(x, opts.Fs, opts.Band(1), opts.Band(2), CompareFieldTrip=true, ComparePMTM=true);  tAll = toc(t0);
        rows(end+1, :) = {M.cell_id(c), r, "matlab_HB_irasa", "matlab", aHB, tHB, st(aHB, "")}; %#ok<AGROW>
        rows(end+1, :) = {M.cell_id(c), r, "fieldtrip_irasa", "matlab", ftI.alphaFT, max(tAll - tHB, 0), st(ftI.alphaFT, ftI.message)}; %#ok<AGROW>
        rows(end+1, :) = {M.cell_id(c), r, "matlab_PM_pmtm", "matlab", pmI.alphaPM, NaN, st(pmI.alphaPM, pmI.message)}; %#ok<AGROW>
    end
end
E = cell2table(rows, VariableNames=["cell_id","rep","estimator","language","alpha_hat","seconds","status"]);
writetable(E, fullfile(out, "estimates_matlab.csv"));

% Python and R arms (same files; no MATLAB needed for these two)
nyq = opts.Fs;  % convert band to cycles/sample for R multitaper
cmdP = sprintf('"%s" "%s" "%s" "%s" "%s" %g %g %g 2>&1', py, fullfile(here, "run_stage1_v001.py"), mf, sdir, fullfile(out, "estimates_python.csv"), opts.Fs, opts.Band);
cmdR = sprintf('"%s" "%s" "%s" "%s" "%s" %g %g 2>&1', opts.RscriptExe, fullfile(here, "run_stage1_v001.R"), mf, sdir, fullfile(out, "estimates_R.csv"), opts.Band / nyq);
t0 = tic; [sP, oP] = system(cmdP); wallP = toc(t0);
t0 = tic; [sR, oR] = system(cmdR); wallR = toc(t0);
if sP ~= 0; error('pilot:python', 'Python runner failed (exit %d):\n%s', sP, oP); end
if sR ~= 0; error('pilot:R', 'R runner failed (exit %d):\n%s', sR, oR); end
fprintf('Python arm wall %.1f s, R arm wall %.1f s\n', wallP, wallR);

% Collate
P = readtable(fullfile(out, "estimates_python.csv"), TextType="string");
R = readtable(fullfile(out, "estimates_R.csv"), TextType="string");
E.cell_id = string(E.cell_id); E.estimator = string(E.estimator); E.language = string(E.language); E.status = string(E.status);
A = [E; P(:, E.Properties.VariableNames); R(:, E.Properties.VariableNames)];
writetable(A, fullfile(out, "estimates_all.csv"));
A = join(A, M(:, ["cell_id","alpha"]));
[g, est, cid] = findgroups(A.estimator, A.cell_id);
ok = A.status == "ok";
S = table(est, cid, splitapply(@(a) a(1), A.alpha, g), splitapply(@sum, ok, g), splitapply(@(s) sum(s ~= "ok"), A.status, g), ...
    VariableNames=["estimator","cell_id","alpha","nOK","nFail"]);
S.median = splitapply(@(a, o) median(a(o)), A.alpha_hat, ok, g);
S.bias   = S.median - S.alpha;
S.sd     = splitapply(@(a, o) std(a(o)), A.alpha_hat, ok, g);
S.rmse   = splitapply(@(a, t, o) sqrt(mean((a(o) - t(o)).^2)), A.alpha_hat, A.alpha, ok, g);
S.secPerSeries = splitapply(@(t) median(t, "omitnan"), A.seconds, g);
S = sortrows(S, ["alpha","estimator"]);
writetable(S, fullfile(out, "summary.csv"));  disp(S)
% Budget: 37 alphas x N in {256..4096} (sum N/4096 = 1.94) x mean reps ~273 (200; 500 at alpha >= 5)
perEst = groupsummary(S, "estimator", "max", "secPerSeries");
proj = sum(perEst.max_secPerSeries, "omitnan") * 37 * 273 * 1.94;
fprintf('Projected full-grid compute (serial, all estimators): %.1f h; on %d cores: %.1f h\n', proj/3600, opts.Cores, proj/3600/opts.Cores);
fprintf('Outputs in %s\n', out);
end

function s = st(a, msg)
if isfinite(a); s = "ok"; else; s = "FAIL " + string(msg); end
end
