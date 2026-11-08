function P = pivotIRASAOptions_v001(runDir)
%PIVOTIRASAOPTIONS_V001 Bias table (median alphaHat - alpha) from checkIRASAOptions_v001 output.
%   Rows generator x fs x alpha, one column per method; failures (thrown errors) listed separately.
%   Writes <runDir>/pivot_v001.txt and returns the table.
S = readtable(fullfile(runDir, 'checkIRASAOptions_v001_summary.csv'), 'TextType', 'string');
P = unstack(S(:, {'generator', 'fs', 'alpha', 'method', 'bias'}), 'bias', 'method');
P = sortrows(P, {'generator', 'fs', 'alpha'});
F = S(S.nFail > 0, {'generator', 'fs', 'alpha', 'method', 'nFail'});
txt = [evalc('disp(P)'), newline, 'Cells with failures:', newline, evalc('disp(F)')];
fid = fopen(fullfile(runDir, 'pivot_v001.txt'), 'w'); c = onCleanup(@() fclose(fid)); fprintf(fid, '%s', txt);
fprintf('%s', txt);
end
