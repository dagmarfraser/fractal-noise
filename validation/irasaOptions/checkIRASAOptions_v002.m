function T = checkIRASAOptions_v002(outDir, nReps)
%CHECKIRASAOPTIONS_V002 checkIRASAOptions_v001 rerun under the fixed default, Resampling="irasa", with "legacy" beside it (v001 kept, unedited).
%   Same generators, grid, methods, seed and reps as v001: exact integrated fGn (no break), noiseXu Phi 0.99 (break), noiseXu Phi 1;
%   alpha 1.5:1:5.5 x fs {100, 240} x N 3000, band 1-20 Hz, nReps per cell (v001 used 10).
%   Every IRASA method is run twice on the same series: Resampling="irasa" (the v2.0.0 default) and Resampling="legacy" (v1.x, what the
%   estimateIRASA help RANGE figures were measured with). The legacy rows let this run reproduce the help figures and pair them with the new ones.
%   Writes <outDir>/checkIRASAOptions_v002.csv (one row per series x method x resampling) and _summary.csv. Stops rather than overwrites.
%   Usage: T = checkIRASAOptions_v002(fullfile(repoRoot, 'validation', 'irasaOptions', 'runs', 'v002'), 10)
arguments
    outDir (1, :) char
    nReps (1, 1) double {mustBeInteger, mustBePositive} = 10
end
if isfile(fullfile(outDir, 'checkIRASAOptions_v002.csv')), error('checkIRASAOptions:exists', 'Output exists in %s; use a new version.', outDir); end
warnState = warning("off", "estimateIRASA:legacyResampling");   % expected on every legacy call; the legacy rows are labelled
cleanupWarn = onCleanup(@() warning(warnState));
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'toolbox'), fullfile(here, '..', 'exactReference'));
fprintf('estimateIRASA: %s\n', which('estimateIRASA'));
if ~isfolder(outDir), mkdir(outDir); end
logF = fullfile(outDir, 'progress.log');
rng(20261108, 'twister');                                        % as v001, so the series are the same
N = 3000; band = [1 20]; alphas = 1.5:1:5.5; fsSet = [100 240];
gens = ["exact", "xu099", "xu1"];
meths = ["irasa_pmtm", "irasa_welch", "irasa_pmtm_pw", "irasa_welch_pw", "whittle"];
rows = {};
for g = gens
    for a = alphas
        switch g
            case "exact", X = exactRefSeries_v001(N, a, nReps);
            case "xu099", X = cell2mat(arrayfun(@(~) noiseXu(N, a, 1, 1, Phi=0.99), 1:nReps, 'UniformOutput', false));
            case "xu1",   X = cell2mat(arrayfun(@(~) noiseXu(N, a, 1, 1, Phi=1), 1:nReps, 'UniformOutput', false));
        end
        k = min(3, floor((a + 1) / 2));
        for fs = fsSet
            for r = 1:nReps
                x = X(:, r);
                for m = meths
                    if m == "whittle", res = "none"; else, res = ["irasa", "legacy"]; end
                    for rs = res
                        est = NaN;  msg = "";
                        try
                            switch m
                                case "irasa_pmtm",     est = estimateIRASA(x, fs, band(1), band(2), Resampling=rs);
                                case "irasa_welch",    est = estimateIRASA(x, fs, band(1), band(2), PSD="welch", Resampling=rs);
                                case "irasa_pmtm_pw",  est = estimateIRASA(x, fs, band(1), band(2), Prewhiten=k, Resampling=rs);
                                case "irasa_welch_pw", est = estimateIRASA(x, fs, band(1), band(2), PSD="welch", Prewhiten=k, Resampling=rs);
                                case "whittle",        est = estimateWhittleAlpha(x, fs, band(1), band(2));
                            end
                        catch ME
                            msg = string(ME.identifier);   % counted as a failure in the summary, never dropped
                        end
                        rows(end+1, :) = {g, a, fs, r, m, rs, k, est, msg}; %#ok<AGROW>
                    end
                end
            end
            appendLog(logF, sprintf('%s %s alpha %.1f fs %d done\n', char(datetime('now')), g, a, fs));
        end
    end
end
T = cell2table(rows, 'VariableNames', {'generator', 'alpha', 'fs', 'rep', 'method', 'resampling', 'kOracle', 'alphaHat', 'error'});
writetable(T, fullfile(outDir, 'checkIRASAOptions_v002.csv'));
S = groupsummary(T, {'generator', 'alpha', 'fs', 'method', 'resampling'}, {@(v) median(v, 'omitnan'), @(v) sum(isnan(v))}, 'alphaHat');
S.Properties.VariableNames(end-1:end) = {'medianAlphaHat', 'nFail'};
S.bias = S.medianAlphaHat - S.alpha;
writetable(S, fullfile(outDir, 'checkIRASAOptions_v002_summary.csv'));
appendLog(logF, sprintf('%s ALL DONE\n', char(datetime('now'))));
end

function appendLog(f, s)
fid = fopen(f, 'a'); c = onCleanup(@() fclose(fid)); fprintf(fid, '%s', s); %#ok<NASGU>
end
