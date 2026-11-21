function T = smokeMFDFA_v001()
%SMOKEMFDFA_V001 Smoke test for mfdfa_v001 / iaaft_v001 before the V10 runs.
%   Monofractal: exact fGn (alpha 0.4) and exact integrated (alpha 4.5, prewhitened K=2) -> dh small.
%   Multifractal: binomial p-model cascade (p = 0.3, 2^12 points; analytic dh = log2(0.7/0.3) = 1.22
%   for q -> +-inf, smaller over q = -5:5) -> dh large; its IAAFT surrogate -> dh small.
here = fileparts(mfilename('fullpath')); addpath(here, fullfile(here, '..', 'exactReference'));
rng(1);
x1 = exactRefSeries_v001(4096, 0.4, 1);
x2 = exactRefSeries_v001(4096, 4.5, 1);
p = 0.3; c = 1; for lev = 1:12, c = reshape([c * p; c * (1 - p)], 1, []); end

x3 = c(:);
[~, d1] = mfdfa_v001(x1);  [~, d2] = mfdfa_v001(x2, K=2);  [~, d3] = mfdfa_v001(x3);  [~, d4] = mfdfa_v001(iaaft_v001(x3));
T = table(["fGn a0.4"; "exact a4.5 (K=2)"; "p-model 0.3"; "IAAFT of p-model"], [d1; d2; d3; d4], 'VariableNames', {'series', 'dh'});
disp(T)
end
