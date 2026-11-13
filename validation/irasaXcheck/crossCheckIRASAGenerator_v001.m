function T = crossCheckIRASAGenerator_v001(opts)
%CROSSCHECKIRASAGENERATOR_V001 Does the IRASA turn-over depend on the noise generator? (knife item K-IRASA1)
%   T = crossCheckIRASAGenerator_v001() crosses generator x alpha x fs x N and estimates alpha with the
%   pre-reg production estimator iraAlphaSigma_v001 (pmtm, whole record) and the toolbox estimateIRASA.
%   Generators:
%     xu099 : noiseXu default (GGM, Phi = 0.99, break f_b = (1-Phi) fs / 2pi) - what benchmarkIRASA_v001-v004 used
%     xu1   : Xu (2019) Appendix A recursion with phi = 1 (pure fractional integration, truncated); local function,
%             because noiseXu forbids Phi = 1
%     exact : exactRefSeries_v001 (k-fold integrated Davies-Harte fGn) - the stage-1 references
%   Each row is scored against nominal alpha and against the band-effective alpha of that generator's analytic
%   spectrum (bandEffectiveAlpha_v001; GGM S = |1 - phi e^{-i lam}|^{-alpha}).
%   Writes <OutDir>/xcheck_v001.csv (one row per series x estimator) and xcheck_summary_v001.csv; logs progress.
arguments
    opts.Alphas double = [2.5 3.5 4 4.5 5.5 6 6.5]
    opts.Fs double = [100 120 240]
    opts.Ns double = [1024 3000]
    opts.Generators string = ["xu099" "xu1" "exact"]
    opts.NRep (1,1) double = 50
    opts.Band (1,2) double = [1 20]
    opts.Hset double = 1.1:0.05:1.9
    opts.MasterSeed (1,1) double = 20261108
    opts.PreRegFunctions (1,1) string = "~/Library/CloudStorage/Dropbox/Brain2Bee/PowerLawSimulationPreReg/src/functions"
    opts.OutDir (1,1) string = ""
end
% fractal-noise v2.0.0: these results were produced with the v1.x estimateIRASA default,
% which is now Resampling="legacy" (a multitaper/Welch slope fit that does not remove
% peaks). Calls are pinned to "legacy" so a rerun reproduces them; see plan v010 F13.
fprintf('Note: estimateIRASA pinned to Resampling="legacy" (v1.x behaviour); legacy warnings suppressed in this session.\n');
warnState = warning("off", "estimateIRASA:legacyResampling");
cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>
here = fileparts(mfilename('fullpath'));
stage1 = fullfile(here, "..", "stage1");  exactDir = fullfile(here, "..", "exactReference");
pre = char(opts.PreRegFunctions); if startsWith(pre, '~'); pre = [char(java.lang.System.getProperty('user.home')) pre(2:end)]; end
assert(isfolder(pre), 'xcheck:path', 'PreRegFunctions not found: %s', pre);
addpath(stage1, exactDir, pre);
assert(~isempty(which('noiseXu')) && ~isempty(which('estimateIRASA')), 'xcheck:path', 'fractalnoise toolbox not on path');
assert(~isempty(which('iraAlphaSigma_v001')), 'xcheck:path', 'iraAlphaSigma_v001 not on path');
if opts.OutDir == ""; opts.OutDir = fullfile(here, "runs", "xcheck_seed" + opts.MasterSeed + "_v001"); end
if ~isfolder(opts.OutDir); mkdir(opts.OutDir); end
logf = fullfile(opts.OutDir, 'progress.log');
plog = @(varargin) appendLog(logf, varargin{:});

[G, A, F, NN] = ndgrid(1:numel(opts.Generators), opts.Alphas, opts.Fs, opts.Ns);
cells = table(opts.Generators(G(:)).', A(:), F(:), NN(:), 'VariableNames', {'gen', 'alpha', 'fs', 'N'});
nC = height(cells);  cells.seed = opts.MasterSeed + (1:nC).';
cells.alphaEff = NaN(nC, 1);
for c = 1:nC
    cells.alphaEff(c) = bandEffectiveAlpha_v001(analyticS(cells.gen(c), cells.alpha(c)), cells.fs(c), opts.Band, cells.N(c));
end
plog('%d cells x %d reps, generators %s, band [%g %g] Hz', nC, opts.NRep, strjoin(opts.Generators, '/'), opts.Band);

pool = gcp('nocreate'); if isempty(pool); pool = parpool; end      % reuse an existing pool; never restart it
wait(parfevalOnAll(pool, @addpath, 0, stage1, exactDir, pre, fileparts(which('noiseXu'))));
out = cell(nC, 1);  t0 = tic;
parfor c = 1:nC
    out{c} = runCell(cells(c, :), opts);
end
R = vertcat(out{:});
plog('done: %d rows in %.1f min, %d FAIL', height(R), toc(t0)/60, nnz(R.status ~= "ok"));
R = join(R, cells(:, {'gen','alpha','fs','N','alphaEff'}));
writetable(R, fullfile(opts.OutDir, 'xcheck_v001.csv'));
ok = R.status == "ok";
[g, gen, est, alpha, fs, N] = findgroups(R.gen(ok), R.estimator(ok), R.alpha(ok), R.fs(ok), R.N(ok));
med = splitapply(@median, R.alpha_hat(ok), g);  sdv = splitapply(@std, R.alpha_hat(ok), g);
aE = splitapply(@(a) a(1), R.alphaEff(ok), g);  n = splitapply(@numel, R.alpha_hat(ok), g);
T = table(gen, est, alpha, fs, N, n, med, sdv, med - alpha, med - aE, aE, ...
    'VariableNames', {'gen','estimator','alpha','fs','N','n','median','sd','bias','biasEff','alphaEff'});
T = sortrows(T, {'estimator','gen','fs','N','alpha'});
writetable(T, fullfile(opts.OutDir, 'xcheck_summary_v001.csv'));
plog('summary written');
end

function appendLog(logf, fmt, varargin)
fid = fopen(logf, 'a');  assert(fid > 0, 'xcheck:log', 'Cannot open %s', logf);
fprintf(fid, ['[%s] ' fmt '\n'], string(datetime('now', 'Format', 'HH:mm:ss')), varargin{:});  fclose(fid);
end

function R = runCell(cl, opts)
rng(cl.seed, 'twister');
n = opts.NRep;  est = ["iraAlphaSigma_v001" "estimateIRASA"];
rep = repelem((1:n).', 2);  estimator = repmat(est.', n, 1);
alpha_hat = NaN(2*n, 1);  status = repmat("ok", 2*n, 1);
for r = 1:n
    try
        x = generate(cl.gen, cl.N, cl.alpha, cl.fs);
    catch ME
        status(2*r-1:2*r) = "FAIL gen: " + ME.message;  continue
    end
    try, alpha_hat(2*r-1) = iraAlphaSigma_v001(x, cl.fs, opts.Band(1), opts.Band(2), opts.Hset);
    catch ME, status(2*r-1) = "FAIL: " + ME.message; end
    try, alpha_hat(2*r) = estimateIRASA(x, cl.fs, opts.Band(1), opts.Band(2), Resampling="legacy");
    catch ME, status(2*r) = "FAIL: " + ME.message; end
end
nr = 2*n;
R = table(repmat(cl.gen, nr, 1), repmat(cl.alpha, nr, 1), repmat(cl.fs, nr, 1), repmat(cl.N, nr, 1), rep, estimator, alpha_hat, status, ...
    'VariableNames', {'gen','alpha','fs','N','rep','estimator','alpha_hat','status'});
end

function x = generate(gen, N, alpha, fs)
switch gen
    case "xu099", x = noiseXu(N, alpha, 1, fs);
    case "xu1",   x = xuPhi1(N, alpha);
    case "exact", x = exactRefSeries_v001(N, alpha, 1);
    otherwise, error('xcheck:gen', 'Unknown generator %s', gen);
end
x = x(:);
end

function dx = xuPhi1(N, alpha)
% Xu (2019) Appendix A GGM_FD with phi = 1 (pure fractional integration of order d = alpha/2), unit variance.
d = alpha / 2;  w = randn(N, 1);
k = (1:N-1).';  h = [1; cumprod((k + d - 1) ./ k)];
np = 2^(fix(log2(2*N)) + 1);
dx = real(ifft(fft(w, np) .* fft(h, np)));  dx = dx(1:N);  dx = dx / std(dx);
end

function S = analyticS(gen, alpha)
switch gen
    case "xu099", S = @(lam) abs(1 - 0.99 * exp(-1i * lam)).^(-alpha);
    case "xu1",   S = @(lam) (2 - 2 * cos(lam)).^(-alpha / 2);
    case "exact", [~, info] = exactRefSeries_v001(8, alpha, 1); S = info.S;
end
end
