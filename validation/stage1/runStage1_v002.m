function S = runStage1_v002(opts)
%RUNSTAGE1_V002 Stage-1 estimator calibration on V0 exact references, 3 languages, two fitting bands.
%   v002 vs v001: every spectral estimator is run in each band of opts.Bands (default "record" then "fixed").
%     record: f_lo = RecordCycles cycles per record (README_IsoSpectralNoise_v037 s9 epoch rule), f_hi = BandFixed(2)
%     fixed : BandFixed (1-20 Hz, as in the BRM work)
%   Band-free estimators (likelihood, DFA, k selection) run once. Manifest carries f_lo_/f_hi_/alpha_eff_<band>.
%   Writes a hub-neutral manifest (with per-cell seeds and band-effective alpha) + series CSVs, runs the
%   MATLAB estimators in parfor (HB/FT/PM via estimateIRASA), then run_stage1_v002.py (joblib) and
%   run_stage1_v002.R (mclapply) on the same files, one language at a time. Collates long-format estimates;
%   summary scores each estimator against nominal alpha AND band-effective alpha, plus k-selection accuracy.
%   Default grid = full stage 1. Pilot: runStage1_v002(Alphas=[0.5 4], Ns=4096, NRep=10, NRepSteep=10).
%   Odd-integer alphas have no exact reference: they are listed as excluded, never silently nudged.
arguments
    opts.Alphas (1,:) double = -0.75:0.25:6.75
    opts.Ns (1,:) double = [256 512 1024 2048 4096]
    opts.NRep (1,1) double = 200
    opts.NRepSteep (1,1) double = 500          % used where alpha >= opts.SteepFrom
    opts.SteepFrom (1,1) double = 5
    opts.MasterSeed (1,1) double = 20261108
    opts.Fs (1,1) double = 100
    opts.BandFixed (1,2) double = [1 20]
    opts.RecordCycles (1,1) double = 20
    opts.Bands (1,:) string {mustBeMember(opts.Bands, ["record" "fixed"])} = ["record" "fixed"]
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
M = table('Size', [nC 8], 'VariableTypes', ["string","double","double","double","double","double","double","string"], ...
    'VariableNames', ["cell_id","alpha","N","k","H","nrep","seed","series_file"]);
B = struct('fixed', repmat(opts.BandFixed, nC, 1), 'record', [opts.RecordCycles * opts.Fs ./ Nn, repmat(opts.BandFixed(2), nC, 1)]);
for L = opts.Bands
    if any(B.(L)(:, 1) >= B.(L)(:, 2)); error('stage1:band', '%s band empty for N = %s', L, mat2str(Nn(B.(L)(:, 1) >= B.(L)(:, 2)).')); end
    M.("f_lo_" + L) = B.(L)(:, 1);  M.("f_hi_" + L) = B.(L)(:, 2);  M.("alpha_eff_" + L) = NaN(nC, 1);
end
for c = 1:nC
    nrep = opts.NRep + (A(c) >= opts.SteepFrom) * (opts.NRepSteep - opts.NRep);
    seed = randi(ms, 2^31 - 1);
    [X, info] = exactRefSeries_v001(Nn(c), A(c), nrep, Stream=RandStream('threefry', Seed=seed));
    id = sprintf("a%+.2f_N%d", A(c), Nn(c));
    M(c, 1:8) = {id, A(c), Nn(c), info.k, info.H, nrep, seed, id + ".csv"};
    for L = opts.Bands; M.("alpha_eff_" + L)(c) = bandEffectiveAlpha_v001(info.S, opts.Fs, B.(L)(c, :), Nn(c)); end
    writematrix(X, fullfile(sdir, id + ".csv"));
end
mf = fullfile(out, "manifest.csv");  writetable(M, mf);
plog('Manifest: %d cells, %d series, master seed %d', nC, sum(M.nrep), opts.MasterSeed);

% ---- MATLAB arm (parfor over cells; reuse any existing pool, never restart it)
if any(opts.Languages == "matlab")
    pool = gcp('nocreate');  if isempty(pool); pool = parpool(opts.Cores); end
    wait(parfevalOnAll(pool, @addpath, 0, opts.FieldTripDir, fullfile(here, "..", "exactReference"), fileparts(which("estimateIRASA"))));
    wait(parfevalOnAll(pool, @ft_defaults, 0));
    R = cell(nC, 1);  t0 = tic;  % parfor over cells: each worker reads its series file once
    fs = opts.Fs;  cid = M.cell_id;  sf = fullfile(sdir, M.series_file);  bl = opts.Bands;
    bb = cell2mat(arrayfun(@(L) [M.("f_lo_" + L), M.("f_hi_" + L)], bl, UniformOutput=false));  % nC x 2*nBands
    parfor c = 1:nC
        X = readmatrix(sf(c));  Tc = cell(size(X, 2), numel(bl));
        for r = 1:size(X, 2); for j = 1:numel(bl); Tc{r, j} = matlabArm(X(:, r), fs, bb(c, 2*j-1:2*j), cid(c), r, bl(j)); end; end
        R{c} = vertcat(Tc{:});
    end
    E = vertcat(R{:});  writetable(E, fullfile(out, "estimates_matlab.csv"));
    plog('MATLAB arm: %d series in %.1f min, %d FAIL rows', sum(M.nrep), toc(t0)/60, sum(E.status ~= "ok"));
end

% ---- Python and R arms (same files; neither needs MATLAB)
exe = struct('python', sprintf('"%s" "%s"', py, fullfile(here, "run_stage1_v003.py")), ...
             'R', sprintf('"%s" "%s"', opts.RscriptExe, fullfile(here, "run_stage1_v003.R")));
for L = intersect(["python" "R"], opts.Languages, "stable")
    for j = 1:numel(opts.Bands)
        subset = "spectral";  if j == 1; subset = "all"; end     % band-free estimators only once
        cmd = sprintf('%s "%s" "%s" "%s" %g %s %d %s 2>&1', exe.(L), mf, sdir, ...
            fullfile(out, sprintf("estimates_%s_%s.csv", L, opts.Bands(j))), opts.Fs, opts.Bands(j), opts.Cores, subset);
        plog('%s arm, band %s (%s) starting', L, opts.Bands(j), subset);  t0 = tic;  [s, o] = system(cmd);
        if s ~= 0; plog('%s arm FAILED (exit %d)', L, s); error('stage1:%s', '%s runner failed (exit %d):\n%s', L, L, s, o); end
        plog('%s arm, band %s: %.1f min', L, opts.Bands(j), toc(t0)/60);
    end
end

% ---- Collate + summarise
S = collateStage1_v002(out);
plog('Done. Outputs in %s', out);
end

function T = matlabArm(x, fs, b, cid, r, band)
rows = cell(3, 9);
try
    t0 = tic;  aHB = estimateIRASA(x, fs, b(1), b(2), Resampling="legacy");  tHB = toc(t0);
    t0 = tic;  [~, ~, ~, ~, ftI, pmI] = estimateIRASA(x, fs, b(1), b(2), CompareFieldTrip=true, ComparePMTM=true, Resampling="legacy");  tAll = toc(t0);
    rows(1, :) = {cid, r, "matlab_HB_irasa", "matlab", aHB, tHB, st(aHB, ""), NaN, band};
    rows(2, :) = {cid, r, "fieldtrip_irasa", "matlab", ftI.alphaFT, max(tAll - tHB, 0), st(ftI.alphaFT, ftI.message), NaN, band};
    rows(3, :) = {cid, r, "matlab_PM_pmtm", "matlab", pmI.alphaPM, NaN, st(pmI.alphaPM, pmI.message), NaN, band};
catch ME
    nm = ["matlab_HB_irasa"; "fieldtrip_irasa"; "matlab_PM_pmtm"];
    for i = 1:3; rows(i, :) = {cid, r, nm(i), "matlab", NaN, NaN, "FAIL " + string(ME.identifier) + ": " + string(ME.message), NaN, band}; end
end
T = cell2table(rows, VariableNames=["cell_id","rep","estimator","language","alpha_hat","seconds","status","k_used","band"]);
end

function s = st(a, msg)
if isfinite(a); s = "ok"; else; s = "FAIL " + string(msg); end
end

function logmsg(f, varargin)
msg = sprintf("[%s] %s", string(datetime("now", Format="HH:mm:ss")), sprintf(varargin{:}));
fprintf('%s\n', msg);  fid = fopen(f, 'a');  fprintf(fid, '%s\n', msg);  fclose(fid);
end
