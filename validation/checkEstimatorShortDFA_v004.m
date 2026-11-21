% checkEstimatorShortDFA_v004.m  At which series lengths does the DFA pre-estimate in estimateWhittleAlpha rest on a single box size?
% Run from a shell:  matlab -batch "run('<this file>')"     Prints one line per length; writes nothing. Seed fixed.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'toolbox'));
rng(20261118, 'twister');
fs = 120;   % at fs 300 these lengths have too few band bins, so the estimator stops with estimateWhittleAlpha:band
for n = [99 100 105 109 110 111 119 120]
    x = cumsum(randn(n, 1));
    s = unique(round(exp(linspace(log(10), log(max(floor(n / 10), 1)), 15))));
    lastwarn('');  warning('off', 'MATLAB:polyfit:PolyNotUnique');
    try
        [a, k, info] = estimateWhittleAlpha(x, fs, 1, 20);
        [~, id] = lastwarn;
        fprintf('n %3d: box sizes %d; alpha %.3f, k %d, alphaDFA %.3f, warning "%s"\n', n, numel(s), a, k, info.alphaDFA, id);
    catch err
        fprintf('n %3d: box sizes %d; error %s\n', n, numel(s), err.identifier);
    end
end
