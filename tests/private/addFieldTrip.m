function addFieldTrip(testCase)
%ADDFIELDTRIP Put FieldTrip on the path from FRACTALNOISE_FIELDTRIP if it
% is not already there. Set but missing -> fatal error; unset -> logged
% notice (the FieldTrip comparison tests then report as Incomplete).
if ~isempty(which("ft_defaults")), return, end
d = string(getenv("FRACTALNOISE_FIELDTRIP"));
if d == ""
    testCase.log(1, "FieldTrip not on path and FRACTALNOISE_FIELDTRIP not set: FieldTrip comparison tests will be filtered.");
    return
end
testCase.fatalAssertTrue(isfolder(d), "FRACTALNOISE_FIELDTRIP is set but is not a folder: " + d);
addpath(d);
evalc("ft_defaults");
end
