function runDFA45Pass_v001(runDir, opts)
%RUNDFA45PASS_V001 Supplementary DFA-4/5 pass on an existing stage-1 run (Kiyono 2015: DFA-m reaches alpha < 2m+1).
%   runDFA45Pass_v001(runDir) launches run_stage1_v004.py and run_stage1_v004.R with subset 'dfa45' IN THE
%   BACKGROUND (system ... &), writing estimates_python_dfa45.csv / estimates_R_dfa45.csv and logs
%   dfa45_python.log / dfa45_R.log in runDir. MATLAB returns at once. When both CSVs exist, re-run
%   collateStage1_v002(runDir) to fold them into summary.csv. Fails loudly if inputs are missing.
arguments
    runDir (1,1) string
    opts.PythonExe (1,1) string = "~/venvs/fnval/bin/python"
    opts.RscriptExe (1,1) string = "/usr/local/bin/Rscript"
    opts.Cores (1,1) double = 10
    opts.Fs (1,1) double = 100
end
here = fileparts(mfilename('fullpath'));
man = fullfile(runDir, "manifest.csv");  sdir = fullfile(runDir, "series");
assert(isfile(man) && isfolder(sdir), 'dfa45:inputs', 'manifest.csv or series/ missing in %s', runDir);
py = regexprep(opts.PythonExe, "^~", getenv("HOME"));
assert(isfile(py) && isfile(opts.RscriptExe), 'dfa45:exe', 'Python or Rscript not found');
q = @(s) """" + s + """";
cmdPy = sprintf('%s %s %s %s %s %g record %d dfa45 > %s 2>&1 &', q(py), q(fullfile(here, "run_stage1_v004.py")), q(man), q(sdir), ...
    q(fullfile(runDir, "estimates_python_dfa45.csv")), opts.Fs, opts.Cores, q(fullfile(runDir, "dfa45_python.log")));
cmdR = sprintf('%s %s %s %s %s %g record %d dfa45 > %s 2>&1 &', q(opts.RscriptExe), q(fullfile(here, "run_stage1_v004.R")), q(man), q(sdir), ...
    q(fullfile(runDir, "estimates_R_dfa45.csv")), opts.Fs, opts.Cores, q(fullfile(runDir, "dfa45_R.log")));
[s1, o1] = system(cmdPy);  [s2, o2] = system(cmdR);
if s1 ~= 0 || s2 ~= 0, error('dfa45:launch', 'Launch failed: py %d (%s), R %d (%s)', s1, o1, s2, o2); end
fprintf('DFA-4/5 pass launched in background for %s; watch dfa45_python.log / dfa45_R.log.\n', runDir);
end
