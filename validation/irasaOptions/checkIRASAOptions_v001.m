function T = checkIRASAOptions_v001(outDir, nReps)
%CHECKIRASAOPTIONS_V001 Bias of estimateIRASA's PSD / Prewhiten options and estimateWhittleAlpha.
%   Generators: exact integrated fGn (no break), noiseXu Phi 0.99 (GGM break), noiseXu Phi 1.
%   Grid: alpha 1.5:1:5.5 x fs {100, 240} x N 3000, band 1-20 Hz, nReps per cell.
%   Methods: IRASA pmtm (toolbox default), IRASA welch, IRASA pmtm + oracle prewhiten,
%   IRASA welch + oracle prewhiten (k = floor((alpha+1)/2), max 3), estimateWhittleAlpha (auto k).
%   Writes <outDir>/checkIRASAOptions_v001.csv (one row per series x method) and _summary.csv.
%   Usage: T = checkIRASAOptions_v001(fullfile(repoRoot,'validation','irasaOptions','runs','v001'), 10)
arguments
    outDir (1, :) char
    nReps (1, 1) double {mustBeInteger, mustBePositive} = 10
end
% fractal-noise v2.0.0: these results were produced with the v1.x estimateIRASA default,
% which is now Resampling="legacy" (a multitaper/Welch slope fit that does not remove
% peaks). Calls are pinned to "legacy" so a rerun reproduces them; see plan v010 F13.
fprintf("Note: estimateIRASA pinned to Resampling=\"legacy\" (v1.x behaviour); legacy warnings suppressed in this session.\n");
warnState = warning("off", "estimateIRASA:legacyResampling");
cleanupWarn = onCleanup(@() warning(warnState)); %#ok<NASGU>
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'toolbox'), fullfile(here, '..', 'exactReference'));
if ~isfolder(outDir), mkdir(outDir); end
logF = fullfile(outDir, 'progress.log');
rng(20261108, 'twister');
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
                est = NaN(1, numel(meths)); msg = strings(1, numel(meths));
                for m = 1:numel(meths)
                    try
                        switch meths(m)
                            case "irasa_pmtm",     est(m) = estimateIRASA(x, fs, band(1), band(2), Resampling="legacy");
                            case "irasa_welch",    est(m) = estimateIRASA(x, fs, band(1), band(2), PSD="welch", Resampling="legacy");
                            case "irasa_pmtm_pw",  est(m) = estimateIRASA(x, fs, band(1), band(2), Prewhiten=k, Resampling="legacy");
                            case "irasa_welch_pw", est(m) = estimateIRASA(x, fs, band(1), band(2), PSD="welch", Prewhiten=k, Resampling="legacy");
                            case "whittle",        est(m) = estimateWhittleAlpha(x, fs, band(1), band(2));
                        end
                    catch ME
                        msg(m) = string(ME.identifier);   % counted as FAIL in the summary, never dropped
                    end
                    rows(end+1, :) = {g, a, fs, r, meths(m), k, est(m), msg(m)}; %#ok<AGROW>
                end
            end
            appendLog(logF, sprintf('%s %s alpha %.1f fs %d done\n', char(datetime('now')), g, a, fs));
        end
    end
end
T = cell2table(rows, 'VariableNames', {'generator', 'alpha', 'fs', 'rep', 'method', 'kOracle', 'alphaHat', 'error'});
writetable(T, fullfile(outDir, 'checkIRASAOptions_v001.csv'));
S = groupsummary(T, {'generator', 'alpha', 'fs', 'method'}, {@(v) median(v, 'omitnan') , @(v) sum(isnan(v))}, 'alphaHat');
S.Properties.VariableNames(end-1:end) = {'medianAlphaHat', 'nFail'};
S.bias = S.medianAlphaHat - S.alpha;
writetable(S, fullfile(outDir, 'checkIRASAOptions_v001_summary.csv'));
appendLog(logF, sprintf('%s ALL DONE\n', char(datetime('now'))));
end

function appendLog(f, s)
fid = fopen(f, 'a'); c = onCleanup(@() fclose(fid)); fprintf(fid, '%s', s);
end
