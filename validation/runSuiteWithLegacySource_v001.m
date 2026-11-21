function R = runSuiteWithLegacySource_v001(legacySrc)
%RUNSUITEWITHLEGACYSOURCE_V001 Run the full test suite with the research code on the path, so the source-equivalence tests run.
%   R = runSuiteWithLegacySource_v001("/path/to/PowerLawSimulationPreReg/src/functions")
%   Sets FRACTALNOISE_LEGACY_SRC to that folder and FRACTALNOISE_REQUIRE_LEGACY to 1 (a missing folder is then an error, not a
%   skipped test), runs tests/, prints the counts and lists every test that did not pass. The environment is restored afterwards.
arguments
    legacySrc (1, 1) string {mustBeFolder}
end
old = [string(getenv("FRACTALNOISE_LEGACY_SRC")), string(getenv("FRACTALNOISE_REQUIRE_LEGACY"))];
cleanup = onCleanup(@() [setenv("FRACTALNOISE_LEGACY_SRC", old(1)), setenv("FRACTALNOISE_REQUIRE_LEGACY", old(2))]);  %#ok<NASGU>
setenv("FRACTALNOISE_LEGACY_SRC", legacySrc);  setenv("FRACTALNOISE_REQUIRE_LEGACY", "1");
root = fileparts(fileparts(mfilename("fullpath")));
R = runtests(fullfile(root, "tests"));
fprintf("\nSUITE WITH LEGACY SOURCE (%s): total %d, passed %d, failed %d, incomplete %d\n", legacySrc, numel(R), nnz([R.Passed]), nnz([R.Failed]), nnz([R.Incomplete]));
for r = R(~[R.Passed]), fprintf("  not passed: %s (failed %d, incomplete %d)\n", r.Name, r.Failed, r.Incomplete); end
end
