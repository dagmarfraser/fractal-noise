function report = checkEnvFN_v001(opts)
%CHECKENVFN_V001 Minimal 3-language smoke test for the fractal-noise validation.
%   report = checkEnvFN_v001  checks MATLAB release/toolboxes, Python and R
%   interpreters + pinned packages, then runs one shared-input round trip
%   (periodogram slope of cumsum(white noise), N=4096) in all three languages
%   and requires agreement to 1e-9, then one HB/FT/PM IRASA call on noiseXu
%   (alpha=2). Errors at the end if anything failed.
%   Pattern: docs/README_RValidationMethodology_v001.md (system() + CSV).
arguments
    opts.PythonExe (1,1) string = "python3"   % e.g. "~/venvs/fnval/bin/python"
    opts.RscriptExe (1,1) string = "Rscript"
    opts.Tol (1,1) double = 1e-9
    opts.FieldTripDir (1,1) string = fullfile(getenv("HOME"), "Library/CloudStorage/Dropbox/MATLAB/fieldtrip-20251201")
end
% fractal-noise v2.0.0: these results were produced with the v1.x estimateIRASA default,
% which is now Resampling="legacy" (a multitaper/Welch slope fit that does not remove
% peaks). Calls are pinned to "legacy" so a rerun reproduces them; see plan v010 F13.
fprintf("Note: estimateIRASA pinned to Resampling=\"legacy\" (v1.x behaviour); legacy warnings suppressed in this session.\n");
warnState = warning("off", "estimateIRASA:legacyResampling");
cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>
here = fileparts(mfilename("fullpath"));
% Expand a leading ~ ourselves: system() quotes the path, so the shell never does.
opts.PythonExe = regexprep(opts.PythonExe, "^~", getenv("HOME"));
opts.RscriptExe = regexprep(opts.RscriptExe, "^~", getenv("HOME"));
fail = strings(0,1);  report = struct();

%% 1. MATLAB
report.matlab = version;
needTbx = ["Signal Processing Toolbox","DSP System Toolbox","Sensor Fusion and Tracking Toolbox","Parallel Computing Toolbox"];
have = string({ver().Name});
missingTbx = needTbx(~ismember(needTbx, have));
fprintf('MATLAB %s\n', report.matlab);
if isMATLABReleaseOlderThan("R2025a"); fail(end+1) = "MATLAB older than R2025a"; end
if ~isempty(missingTbx); fail(end+1) = "Missing toolboxes: " + strjoin(missingTbx, ", "); end
if isempty(which('noiseXu')); fail(end+1) = "fractalnoise toolbox (noiseXu) not on path"; end
if isempty(which('ft_defaults')) && isfile(fullfile(opts.FieldTripDir, "ft_defaults.m")); addpath(opts.FieldTripDir); ft_defaults; end
if isempty(which('ft_defaults')); fail(end+1) = "FieldTrip not found (FieldTripDir = " + opts.FieldTripDir + ")";
else; report.fieldtrip = string(fileparts(which('ft_defaults'))); fprintf('FieldTrip %s\n', report.fieldtrip); end

%% 2. Python: interpreter + pinned packages
req = readlines(fullfile(here, "requirements_v001.txt"));
req = strtrim(extractBefore(req + "#", "#"));  req = req(contains(req, "=="));
pyName = extractBefore(req, "==");  pyPin = strtrim(extractAfter(req, "=="));
pyAll = [pyName; "nolitsa"];
py = "import json, importlib.metadata as m" + newline + "out={}" + newline + ...
     "for p in " + jsonencode(cellstr(pyAll)) + ":" + newline + ...
     "    try: out[p]=m.version(p); __import__(p)" + newline + ...
     "    except Exception as e: out[p]='MISSING: '+type(e).__name__" + newline + ...
     "import sys; out['_python']=sys.version.split()[0]; print('JSON|' + json.dumps(out))";
[ok, out] = runScript(opts.PythonExe, py, ".py");
if ~ok; fail(end+1) = "Python not runnable: " + out;
else
    j = splitlines(out);  j = extractAfter(j(startsWith(j, "JSON|")), "JSON|");  % ignore import warnings
    if isempty(j); error('checkEnvFN:pyParse', 'No JSON line from Python. Raw output:\n%s', out); end
    v = jsondecode(j(end));  report.python = v;  fprintf('Python %s\n', v.x_python);
    fail = [fail; pinCheck(pyAll, [pyPin; ""], @(p) getf(v, p), "py")];
end

%% 3. R: interpreter + locked packages
lockF = fullfile(here, "r_lock_v001.csv");
if ~isfile(lockF); fail(end+1) = "r_lock_v001.csv missing: run Rscript installR_fnval_v001.R first";
else
    lk = readtable(lockF, TextType="string");
    r = "pk <- c(" + strjoin('"' + lk.package + '"', ",") + ")" + newline + ...
        "v <- vapply(pk, function(p) if (suppressMessages(requireNamespace(p, quietly=TRUE))) as.character(packageVersion(p)) else 'MISSING', '')" + newline + ...
        "cat(paste0('PKG|', c(as.character(getRversion()), v)), sep='\n')";
    [ok, out] = runScript(opts.RscriptExe, r, ".R");
    if ~ok; fail(end+1) = "R not runnable: " + out;
    else
        L = strtrim(splitlines(string(out)));  L = extractAfter(L(startsWith(L, "PKG|")), "PKG|");  % ignore R startup chatter
        report.R = L(1);  fprintf('R %s\n', L(1));
        got = L(2:end);
        fail = [fail; pinCheck(lk.package, lk.version, @(p) got(lk.package == p), "R")];
    end
end

%% 4. Shared-input round trip
rng(1);  x = randn(4096, 1);
dataF = fullfile(tempdir, "fn_smoke_x.csv");  writematrix(x, dataF);
c = onCleanup(@() delete(dataF));
aM = slopeM(x);
pD = strrep(dataF, '\', '/');
pyRT = "import numpy as np" + newline + "x=np.loadtxt('" + pD + "')" + newline + ...
       "y=np.cumsum(x); N=len(y); k=np.arange(1,N//2); P=np.abs(np.fft.fft(y))[k]**2" + newline + ...
       "print('%.17g' % -np.polyfit(np.log(k/N), np.log(P), 1)[0])";
rRT  = "x <- scan('" + pD + "', quiet=TRUE); y <- cumsum(x); N <- length(y); k <- 1:(N%/%2-1)" + newline + ...
       "P <- Mod(fft(y))[k+1]^2; cat(sprintf('%.17g', -coef(lm(log(P) ~ log(k/N)))[2]))";
[okP, oP] = runScript(opts.PythonExe, pyRT, ".py");
[okR, oR] = runScript(opts.RscriptExe, rRT, ".R");
lastLine = @(o) subsref(splitlines(o), substruct('()', {numel(splitlines(o))}));
aP = str2double(lastLine(oP));  aR = str2double(lastLine(oR));  % last line only: skips any interpreter chatter
report.roundTrip = table(aM, aP, aR, VariableNames=["MATLAB","Python","R"]);
disp(report.roundTrip)
if ~okP || ~isfinite(aP) || abs(aP - aM) > opts.Tol; fail(end+1) = "Python round trip disagrees/failed: " + oP; end
if ~okR || ~isfinite(aR) || abs(aR - aM) > opts.Tol; fail(end+1) = "R round trip disagrees/failed: " + oR; end

%% 5. Minimal three-way IRASA: homebrew (HB) vs FieldTrip (FT) vs direct pmtm (PM)
% Existence/finiteness only; agreement is the validation's job (V1e), not this smoke test's.
if ~isempty(which('noiseXu')) && ~isempty(which('ft_defaults'))
    rng(2);  fs = 100;  xn = noiseXu(4096, 2, 1, fs);
    [aHB, ~, ~, ~, ftI, pmI] = estimateIRASA(xn(:), fs, 1, 20, CompareFieldTrip=true, ComparePMTM=true, Resampling="legacy");
    report.irasa = table(2, aHB, ftI.alphaFT, pmI.alphaPM, VariableNames=["target","HB","FT","PM"]);
    disp(report.irasa)
    if any(~isfinite([aHB ftI.alphaFT pmI.alphaPM])); fail(end+1) = "IRASA three-way not finite: " + ftI.message + " " + pmI.message; end
end

%% Verdict: fail loud
report.failures = fail;
if isempty(fail); fprintf('\nENV OK: MATLAB, Python and R all pass.\n');
else; error('checkEnvFN:fail', 'Environment check FAILED (%d):\n  - %s', numel(fail), strjoin(fail, newline + "  - "));
end
end

function a = slopeM(x)
y = cumsum(x);  N = numel(y);  k = (1:N/2-1)';
P = abs(fft(y)).^2;  P = P(k+1);
b = polyfit(log(k/N), log(P), 1);  a = -b(1);
end

function [ok, out] = runScript(exe, code, ext)
f = string(tempname) + ext;  fid = fopen(f, 'w');  fprintf(fid, '%s', code);  fclose(fid);
[s, out] = system(sprintf('"%s" "%s" 2>&1', exe, f));  delete(f);
ok = s == 0;  out = strtrim(string(out));
end

function f = pinCheck(names, pins, getter, tag)
f = strings(0,1);
for i = 1:numel(names)
    g = string(getter(names(i)));
    if isempty(g) || startsWith(g, "MISSING"); f(end+1,1) = sprintf("%s %s: MISSING", tag, names(i)); %#ok<AGROW>
    elseif pins(i) ~= "" && g ~= pins(i); f(end+1,1) = sprintf("%s %s: installed %s, pinned %s", tag, names(i), g, pins(i)); %#ok<AGROW>
    end
end
end

function v = getf(s, p)
fn = matlab.lang.makeValidName(p);
if isfield(s, fn); v = s.(fn); else; v = "MISSING"; end
end
