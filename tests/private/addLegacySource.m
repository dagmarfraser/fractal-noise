function addLegacySource(testCase)
%ADDLEGACYSOURCE Put the originating research code on the path for the
% source-equivalence tests, from the FRACTALNOISE_LEGACY_SRC environment
% variable (the folder holding iraAlphaSigma_v003.m, XuNoise_v002.m,
% generateCustomNoise_v004.m, alphaPWLW_v001.m).
%
% Fail loud, never silent:
%   - variable set but folder missing -> fatal error;
%   - variable unset -> a logged notice that the source-equivalence tests
%     will be filtered (they then report as Incomplete, naming the variable);
%   - FRACTALNOISE_REQUIRE_LEGACY=1 and variable unset -> fatal error
%     (use this on a machine that is meant to check bit-identity).
d = string(getenv("FRACTALNOISE_LEGACY_SRC"));
if d == ""
    testCase.fatalAssertNotEqual(string(getenv("FRACTALNOISE_REQUIRE_LEGACY")), "1", ...
        "FRACTALNOISE_REQUIRE_LEGACY=1 but FRACTALNOISE_LEGACY_SRC is not set.");
    testCase.log(1, "FRACTALNOISE_LEGACY_SRC not set: source-equivalence tests will be filtered.");
    return
end
testCase.fatalAssertTrue(isfolder(d), "FRACTALNOISE_LEGACY_SRC is set but is not a folder: " + d);
addpath(d);
testCase.addTeardown(@() rmpath(d));
end
