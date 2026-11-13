function S = runStage1_v001(opts)
%RUNSTAGE1_V001 Stage-1 estimator calibration on V0 exact references, 3 languages (supersedes runStage1Pilot_v001).
%   Writes a hub-neutral manifest (with per-cell seeds and band-effective alpha) + series CSVs, runs the
%   MATLAB estimators in parfor (HB/FT/PM via estimateIRASA), then run_stage1_v002.py (joblib) and
%   run_stage1_v002.R (mclapply) on the same files, one language at a time. Collates long-format estimates;
%   summary scores each estimator against nominal alpha AND band-effective alpha, plus k-selection accuracy.
%   Default grid = full stage 1. Pilot: runStage1_v001(Alphas=[0.5 4], Ns=4096, NRep=10, NRepSteep=10).
%   Odd-integer alphas have no exact reference: they are listed as excluded, never silently nudged.
arguments
    opts.Alphas (1,:) double = -0.75:0.25:6.75
    opts.Ns (1,:) double = [256 512 1024 2048 4096]
    opts.NRep (1,1) double = 200
    opts.NRepSteep (1,1) double = 500          % used where alpha >= opts.SteepFrom
    opts.SteepFrom (1,1) double = 5
    opts.MasterSeed (1,1) double = 20261108
    opts.Fs (1,1) double = 100
    opts.Band (1,2) double = [1 20]
    opts.PythonExe (1,1) string = "~/venvs/fnval/bin/python"
    opts.RscriptExe (1,1) string = "/usr/local/bin/Rscript"
    opts.OutDir (1,1) string = ""
    opts.FieldTripDir (1,1) string = fullfile(getenv("HOME"), "Library/CloudStorage/Dropbox/MATLAB/fieldtrip-20251201")
    opts.Cores (1,1) double = 10
    opts.Languages (1,:) string = ["matlab" "python" "R"]
end
% fractal-noise v2.0.0: these results were produced with the v1.x estimateIRASA default,
% which is now Resampling="legacy" (a multitaper/Welch slope fit that does not remove
% peaks). Calls are pinned to "legacy" so a rerun reproduces them; see plan v010 F13.
fprintf('Note: estimateIRASA pinned to Resampling="legacy" (v1.x behaviour); legacy warnings suppressed in this session.\n');
warnState = warning("off", "estimateIRASA:legacyResampling");
cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>
here = fileparts(mfilename("fullpath"));
addpath(fullfile(here, "..", "exactReference"));
if isempty(which('ft_defaults')); addpath(opts.FieldTripDir); ft_defaults; end
if isempty(which('estimateIRASA')); error('stage1:path', 'fractalnoise toolbox (estimateIRASA) not on path.'); end
out = opts.OutDir;  if out == ""; out = fullfile(here, "runs", "stage1_" + string(datetime("now", Format="yyyyMMdd_HHmmss"))); end
mkdir(out);  sdir = fullfile(out, "series");  mkdir(sdir);
logf = fullfile(out, "progress.log");  plog = @(varargin) logmsg(logf, varargin{:});
py = regexprep(opts.PythonExe, "^~", getenv("HOME"));

% ---- Grid, exclusions, manifest + series
odd = opts.Alphas(abs(opts.Alphas - round(opts.Alphas)) < 1e-12 & mod(round(opts.Alphas), 2) == 1);
alphas = setdiff(opts.Alphas, odd);
if ~isempty(odd); plog('Excluded (odd-integer alpha, no exact reference): %s', mat2str(odd)); end
[A, Nn] = ndgrid(alphas, opts.Ns);  A = A(:);  Nn = Nn(:);  nC = numel(A);
ms = RandStream('threefry', Seed=opts.MasterSeed);
M = table('Size', [nC 9], 'VariableTypes', ["string","double","double","double","double","double","double","string","double"], ...
    'VariableNames', ["cell_id","alpha","N","k","H","nrep","seed","series_file","alpha_eff"]);
for c = 1:nC
    nrep = opts.NRep + (A(c) >= opts.SteepFrom) * (opts.NRepSteep - opts.NRep);
    seed = randi(ms, 2^31 - 1);
    [X, info] = exactRefSeries_v001(Nn(c), A(c), nrep, Stream=RandStream('threefry', Seed=seed));
    id = sprintf("a%+.2f_N%d", A(c), Nn(c));
    M(c, :) = {id, A(c), Nn(c), info.k, info.H, nrep, seed, id + ".csv", bandEffectiveAlpha_v001(info.S, opts.Fs, opts.Band, Nn(c))};
    writematrix(X, fullfile(sdir, id + ".csv"));
end
mf = fullfile(out, "manifest.csv");  writetable(M, mf);
plog('Manifest: %d cells, %d series, master seed %d', nC, sum(M.nrep), opts.MasterSeed);

% ---- MATLAB arm (parfor over series; reuse any existing pool, never restart it)
if any(opts.Languages == "matlab")
    pool = gcp('nocreate');  if isempty(pool); pool = parpool(opts.Cores); end
    parfevalOnAll(pool, @addpath, 0, opts.FieldTripDir, fullfile(here, "..", "exactReference"));
    wait(parfevalOnAll(pool, @ft_defaults, 0));
    R = cell(nC, 1);  t0 = tic;  % parfor over cells: each worker reads its series file once
    fs = opts.Fs;  b = opts.Band;  cid = M.cell_id;  sf = fullfile(sdir, M.series_file);
    parfor c = 1:nC
        X = readmatrix(sf(c));  Tc = cell(size(X, 2), 1);
        for r = 1:size(X, 2); Tc{r} = matlabArm(X(:, r), fs, b, cid(c), r); end
        R{c} = vertcat(Tc{:});
    end
    E = vertcat(R{:});  writetable(E, fullfile(out, "estimates_matlab.csv"));
    plog('MATLAB arm: %d series in %.1f min, %d FAIL rows', sum(M.nrep), toc(t0)/60, sum(E.status ~= "ok"));
end

% ---- Python and R arms (same files; neither needs MATLAB)
cmd.python = sprintf('"%s" "%s" "%s" "%s" "%s" %g %g %g %d 2>&1', py, fullfile(here, "run_stage1_v002.py"), mf, sdir, ...
    fullfile(out, "estimates_python.csv"), opts.Fs, opts.Band, opts.Cores);
cmd.R = sprintf('"%s" "%s" "%s" "%s" "%s" %g %g %d 2>&1', opts.RscriptExe, fullfile(here, "run_stage1_v002.R"), mf, sdir, ...
    fullfile(out, "estimates_R.csv"), opts.Band / opts.Fs, opts.Cores);
for L = intersect(["python" "R"], opts.Languages, "stable")
    plog('%s arm starting', L);  t0 = tic;  [s, o] = system(cmd.(L));
    if s ~= 0; plog('%s arm FAILED (exit %d)', L, s); error('stage1:%s', '%s runner failed (exit %d):\n%s', L, L, s, o); end
    plog('%s arm: %.1f min', L, toc(t0)/60);
end

% ---- Collate + summarise
S = collateStage1_v001(out);
plog('Done. Outputs in %s', out);
end

function T = matlabArm(x, fs, b, cid, r)
rows = cell(3, 8);
try
    t0 = tic;  aHB = estimateIRASA(x, fs, b(1), b(2), Resampling="legacy");  tHB = toc(t0);
    t0 = tic;  [~, ~, ~, ~, ftI, pmI] = estimateIRASA(x, fs, b(1), b(2), CompareFieldTrip=true, ComparePMTM=true, Resampling="legacy");  tAll = toc(t0);
    rows(1, :) = {cid, r, "matlab_HB_irasa", "matlab", aHB, tHB, st(aHB, ""), NaN};
    rows(2, :) = {cid, r, "fieldtrip_irasa", "matlab", ftI.alphaFT, max(tAll - tHB, 0), st(ftI.alphaFT, ftI.message), NaN};
    rows(3, :) = {cid, r, "matlab_PM_pmtm", "matlab", pmI.alphaPM, NaN, st(pmI.alphaPM, pmI.message), NaN};
catch ME
    nm = ["matlab_HB_irasa"; "fieldtrip_irasa"; "matlab_PM_pmtm"];
    for i = 1:3; rows(i, :) = {cid, r, nm(i), "matlab", NaN, NaN, "FAIL " + string(ME.identifier) + ": " + string(ME.message), NaN}; end
end
T = cell2table(rows, VariableNames=["cell_id","rep","estimator","language","alpha_hat","seconds","status","k_used"]);
end

function s = st(a, msg)
if isfinite(a); s = "ok"; else; s = "FAIL " + string(msg); end
end

function logmsg(f, varargin)
msg = sprintf("[%s] %s", string(datetime("now", Format="HH:mm:ss")), sprintf(varargin{:}));
fprintf('%s\n', msg);  fid = fopen(f, 'a');  fprintf(fid, '%s\n', msg);  fclose(fid);
end
