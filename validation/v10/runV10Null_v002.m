function S = runV10Null_v002(outDir, nRep)
%RUNV10NULL_V002 = v001 plus dhPos = h(1) - h(5) (q > 0 width) for the original and its IAAFT, and the IAAFT
%   95th percentile of each (for v002 diagnostics). Same seed, so dh columns reproduce v001.
%RUNV10NULL_V001 V10 null: MFDFA width dh on exact (monofractal) references, per alpha x N.
%   Prewhitened as on real data: K from estimateWhittleAlpha (toolbox), MFDFA order 2, q = -5:5.
%   Alpha 2.5:1:5.5 (plus 4.0), N = 1000/2000/4000 (Fraser trials are ~1100 samples at 240 Hz),
%   fs-free (MFDFA is scale-based; K uses a 0.01-0.2 cycles/sample band). IAAFT is applied to the
%   DIFFERENCED series (raw steep series break phase surrogates, V3), as on real data; on a monofractal
%   series the IAAFT dh distribution should match the original's.
%   Writes <outDir>/v10null_v002.csv (one row per series) and v10null_summary_v002.csv (median, 95th pct).
arguments
    outDir (1,:) char
    nRep (1,1) double {mustBeInteger, mustBePositive} = 200
end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'toolbox'), fullfile(here, '..', 'exactReference'), here);
if ~isfolder(outDir), mkdir(outDir); end
rng(20261109, 'twister');
alphas = [2.5 3.5 4.0 4.5 5.5];  Ns = [1000 2000 4000];
rows = {};
for a = alphas
    for N = Ns
        X = exactRefSeries_v001(N, a, nRep);
        for r = 1:nRep
            x = X(:, r);  st = "ok";  k = NaN; dh = NaN; h2 = NaN; dhS = NaN; dhP = NaN; dhPS = NaN;
            try
                [~, k] = estimateWhittleAlpha(x, 1, 0.01, 0.2);
                xd = x;  for j = 1:k, xd = diff(xd); end
                [h, dh, info] = mfdfa_v001(xd);  h2 = info.h2 + k;   % h2 reported for the undifferenced series
                q = info.q;  dhP = h(q == 1) - h(q == 5);
                hs = mfdfa_v001(iaaft_v001(xd));  dhS = hs(1) - hs(end);  dhPS = hs(q == 1) - hs(q == 5);
            catch ME
                st = string(ME.identifier) + ": " + string(ME.message);
            end
            rows(end+1, :) = {a, N, r, k, h2, dh, dhS, dhP, dhPS, st}; %#ok<AGROW>
        end
        fprintf('[%s] alpha %.1f N %d done\n', char(datetime('now')), a, N);
    end
end
T = cell2table(rows, 'VariableNames', {'alpha','N','rep','K','h2','dh','dhIAAFT','dhPos','dhPosIAAFT','status'});
writetable(T, fullfile(outDir, 'v10null_v002.csv'));
S = groupsummary(T, {'alpha','N'}, {@(v) median(v,'omitnan'), @(v) prctile(v,95), @(v) sum(isnan(v))}, {'dh','dhIAAFT','dhPos','dhPosIAAFT'});
writetable(S, fullfile(outDir, 'v10null_summary_v002.csv'));
nF = nnz(T.status ~= "ok");
if nF > 0, warning('v10:fails', '%d series failed; see status column.', nF); end
end
