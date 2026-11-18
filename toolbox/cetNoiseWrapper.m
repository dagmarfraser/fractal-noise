function [n, info] = cetNoiseWrapper(options)
%CETNOISEWRAPPER Surrogate noise for the Central England Temperature series: a worked example of a non-kinematic record.
%
%   [n, info] = cetNoiseWrapper()
%   [n, info] = cetNoiseWrapper(FirstYear=1900, LastYear=2020, K=2, AlphaBand=[0.5 5])
%   cet = cetNoiseWrapper(LoadOnly=true)      % just the monthly series
%
%   Loads the monthly HadCET series bundled in examples/data (fs = 12 samples per year, so frequencies are in cycles per
%   year), removes the annual cycle with harmonicTemplate (f0 = 1 per year), and builds the surrogate with shapeNoise
%   over an explicit band. It shows the unit-agnostic path: nothing here is in Hz or tied to hand movement.
%
%   This is an EXAMPLE, outside the known-truth checks in the shapeNoise help (RANGE: fs 120-300 Hz, N 1000 and 3000, band 1-20 Hz):
%   the series is about 4400 months in cycles per year, and the search cost grows with N. Use FirstYear to keep N in the low thousands.
%
% INPUTS (all name-value):
%   FirstYear - first calendar year (default 1659, the start of the series).
%   LastYear  - last year (default [] = the last year with all twelve months present in the bundled file).
%   K         - number of harmonics of the annual cycle stripped (default 2).
%   Window    - template window in years (cycles of f0; default 4).
%   AlphaBand - [fLow fHigh] in cycles per year, 0 < fLow < fHigh < 6 (default [0.5 5]).
%   File      - path to a HadCET monthly file in the Met Office format (e.g. a newer download of meantemp_monthly_totals.txt); default the bundled one.
%   LoadOnly  - true returns the struct of the loaded series and does nothing else.
%   ShapeArgs - cell array of further shapeNoise name-value pairs.
%
% OUTPUTS:
%   n    - Surrogate (info.shape.L months).
%   info - Struct: cet (year, month, value; fs = 12), template, residual, shape.
%
% DATA: Hadley Centre Central England Temperature (HadCET), monthly mean, degrees Celsius, downloaded 2026-11-17 from
% www.metoffice.gov.uk/hadobs/hadcet/data/meantemp_monthly_totals.txt (md5 837fb3a2ec05dc52c744cd5c11374630).
% Contains public sector information licensed under the Open Government Licence v3.0. Source: Met Office Hadley
% Centre, www.metoffice.gov.uk/hadobs. When publishing work that uses it, cite: Parker, D.E., Legg, T.P. and Folland, C.K.
% (1992). A new daily Central England Temperature series, 1772-1991. Int. J. Climatol. 12, 317-342 (the monthly series
% before 1772 is Manley 1974, Q.J.R. Meteorol. Soc. 100, 389-405; see the Met Office page for the full statement).
% The series is revised each year: a newer file may differ from the bundled one.
%
% See also: shapeNoiseWrapper, harmonicTemplate, shapeNoise

arguments
    options.FirstYear (1, 1) double {mustBeInteger} = 1659
    options.LastYear double {mustBeScalarOrEmpty} = []
    options.K (1, 1) double {mustBeInteger, mustBePositive} = 2
    options.Window (1, 1) double = 4
    options.AlphaBand (1, 2) double = [0.5 5]
    options.File (1, 1) string = ""
    options.LoadOnly (1, 1) logical = false
    options.ShapeArgs (1, :) cell = {}
end
cet = loadHadCET(options.File, options.FirstYear, options.LastYear);
if options.LoadOnly
    n = cet;  info = [];  return
end
[n, w] = shapeNoiseWrapper(cet.value, cet.fs, F0=1, K=options.K, Window=options.Window, ...
    AlphaBand=options.AlphaBand, ShapeArgs=options.ShapeArgs);
info = struct("cet", cet, "template", w.template, "residual", w.residual, "shape", w.shape);
end

function cet = loadHadCET(f, firstYear, lastYear)
if f == "", f = fullfile(fileparts(mfilename("fullpath")), "examples", "data", "HadCET_monthly_downloaded_20261117.txt"); end
if ~isfile(f), error("cetNoiseWrapper:noData", "The HadCET file is missing: %s", f); end
txt = splitlines(string(fileread(f)));
rows = zeros(0, 13);
for k = 1:numel(txt)
    v = sscanf(char(txt(k)), "%f")';
    if numel(v) == 14 && v(1) >= 1659 && v(1) == round(v(1)), rows(end + 1, :) = v(1:13); end %#ok<AGROW>
end
if isempty(rows), error("cetNoiseWrapper:parse", "No data rows found in %s.", f); end
complete = all(rows(:, 2:13) > -99, 2);
lastComplete = rows(find(complete, 1, "last"), 1);
if isempty(lastYear), lastYear = lastComplete; end
if lastYear > lastComplete, error("cetNoiseWrapper:lastYear", "LastYear %d is after the last complete year in the file (%d).", lastYear, lastComplete); end
sel = rows(:, 1) >= firstYear & rows(:, 1) <= lastYear;
if ~any(sel), error("cetNoiseWrapper:years", "No years in [%d, %d].", firstYear, lastYear); end
R = rows(sel, :);
if any(R(:, 2:13) <= -99, "all"), error("cetNoiseWrapper:missing", "Missing months inside the chosen years."); end
T = R(:, 2:13)';
cet = struct("year", repelem(R(:, 1), 12), "month", repmat((1:12)', size(R, 1), 1), "value", T(:), "fs", 12, "unit", "degC, cycles per year");
end
