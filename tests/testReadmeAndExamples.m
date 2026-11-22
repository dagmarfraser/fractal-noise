classdef testReadmeAndExamples < matlab.unittest.TestCase
%TESTREADMEANDEXAMPLES Every README code block, example and documentation page runs as written.
%
% Toolbox: fractalnoise
%
% The README code blocks run top to bottom in one workspace, from the repository folder, as a reader would paste them.
% The block that calls runtests is skipped (it would run this suite again). The quick start and the worked story must run
% without a warning, and the stated results are checked: K = 2 and a knee below the band in the story, K = 2, alpha 0.39 and H 0.70 for HadCET.
% Examples and documentation pages run in their own workspaces with figures hidden; the knee example's three cases are
% checked against what its text says.

    properties
        Root
    end

    methods (TestClassSetup)
        function setup(testCase)
            testCase.Root = fullfile(fileparts(mfilename("fullpath")), "..");
            testCase.assumeTrue(isfile(fullfile(testCase.Root, "README.md")), "README.md not beside tests/ (e.g. a temporary copy); skipping.");
            oldPath = path;  testCase.addTeardown(@() path(oldPath));
            oldDir = pwd;  testCase.addTeardown(@() cd(oldDir));
            addpath(fullfile(testCase.Root, "toolbox"));
            oldVis = get(0, "DefaultFigureVisible");  set(0, "DefaultFigureVisible", "off");
            testCase.addTeardown(@() set(0, "DefaultFigureVisible", oldVis));
            testCase.addTeardown(@() close("all"));
        end
    end

    methods (Test)
        function testReadmeBlocksRunAsWritten(testCase)
            blocks = readmeBlocks(fullfile(testCase.Root, "README.md"));
            testCase.verifyGreaterThanOrEqual(numel(blocks), 5);
            cd(testCase.Root);
            code = strjoin(blocks(~contains(blocks, "runtests(")), newline);
            ws = testCase.verifyWarningFree(@() runIsolated(code), "README code raised a warning.");
            testCase.verifyEqual(ws.K, 2, "Worked story: K stated as 2.");
            testCase.verifyFalse(ws.knee.kneeInBand, "Worked story: knee stated to lie below the band.");
            testCase.verifyGreaterThan(ws.knee.kneeHz, 0.15);
            testCase.verifyLessThan(ws.knee.kneeHz, 0.5, "Worked story: knee stated near 0.3 Hz.");
            testCase.verifyEqual(ws.S.nChannels, 100);
            testCase.verifyEqual(ws.Kc, 2, "HadCET example: Kc stated as 2.");
            testCase.verifyEqual(ws.alphaCET, 0.39, "AbsTol", 0.005, "HadCET anomalies: alpha stated as 0.39.");
            testCase.verifyEqual(ws.H, 0.70, "AbsTol", 0.005, "HadCET anomalies: H stated as 0.70.");
        end

        function testExamplesRun(testCase)
            files = dir(fullfile(testCase.Root, "toolbox", "examples", "*.m"));
            testCase.verifyNotEmpty(files);
            for f = files'
                try
                    runIsolated(fileread(fullfile(f.folder, f.name)), fullfile(f.folder, f.name));
                catch err
                    testCase.verifyFail(f.name + " failed: " + err.identifier + ": " + err.message);
                end
            end
        end

        function testDocPagesRun(testCase)
            files = dir(fullfile(testCase.Root, "toolbox", "doc", "*.m"));
            testCase.verifyNotEmpty(files);
            for f = files'
                try
                    runIsolated(fileread(fullfile(f.folder, f.name)), fullfile(f.folder, f.name));
                catch err
                    testCase.verifyFail(f.name + " failed: " + err.identifier + ": " + err.message);
                end
            end
        end

        function testKneeExampleMatchesItsText(testCase)
            f = fullfile(testCase.Root, "toolbox", "examples", "KneeCheckWorkedExample.m");
            ws = runIsolated(fileread(f), f);
            testCase.verifyFalse(ws.resultSubFloor.kneeInBand, "Case 1: knee below the band.");
            testCase.verifyEqual(ws.resultSubFloor.kneeHz, 0.28, "AbsTol", 0.1, "Case 1: knee recovered near 0.28 Hz.");
            testCase.verifyTrue(ws.resultInBand.kneeInBand, "Case 2: knee inside the band.");
            testCase.verifyLessThan(ws.resultInBand.alphaStd, 4.5 - 0.5, "Case 2: the fixed-band alpha is biased low.");
            testCase.verifyLessThan(ws.resultInBand.alphaWide, 4.5 - 0.5, "Case 2: the wide-band fixed alpha is biased low too.");
            testCase.verifyEqual(string(ws.resultFlat.verdict), "none", "Case 3: break-free noise has no knee.");
        end
    end
end

function blocks = readmeBlocks(file)
% The fenced matlab code blocks of a Markdown file, in order.
txt = fileread(file);
tok = regexp(txt, '```matlab\n(.*?)```', 'tokens');
blocks = strings(1, numel(tok));
for i = 1:numel(tok), blocks(i) = string(tok{i}{1}); end
end

function zzWs = runIsolated(zzCode, zzFile)
% Runs code in this function's own workspace and returns its variables. Example and doc files may define local
% functions, so they are run as scripts (run) rather than eval'ed.
if nargin > 1
    evalc("run(zzFile)");
else
    evalc(zzCode);
end
zzWs = struct();
for zzV = string(who)'
    if ~startsWith(zzV, "zz"), zzWs.(zzV) = eval(zzV); end
end
end
