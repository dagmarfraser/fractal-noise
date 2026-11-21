function [n, info] = shapeNoiseWrapper(x, fs, options)
%SHAPENOISEWRAPPER One call from a record to a surrogate noise: optional harmonic template, then shapeNoise.
%
%   [n, info] = shapeNoiseWrapper(x, fs, F0=f0, K=2, AlphaBand=[lo hi])
%   [n, info] = shapeNoiseWrapper(x, fs, F0Band=[fLow fHigh], K=1)   % f0 estimated from the record
%   [n, info] = shapeNoiseWrapper(x, fs, AlphaBand="kinematic")      % no periodic part to remove: shapeNoise only
%
%   Does two things in order, each documented in its own function:
%     1. harmonicTemplate(x, fs, f0, K, Window=...) removes a periodic component, if F0 or F0Band is given;
%     2. shapeNoise(residual, fs, Alpha, AlphaBand=...) builds the end-matched random-phase surrogate.
%   It adds no method of its own. If you need an option that is not here, call the two functions directly.
%
% INPUTS:
%   x    - Real N-by-C matrix (one column per channel) or a vector; finite.
%   fs   - Sampling rate in cycles per unit time (any consistent unit; F0, F0Band and AlphaBand use the same unit).
%   F0         - Name-value. Fundamental frequency, a scalar. Default [] (none).
%   F0Band     - Name-value. [fLow fHigh]; estimate f0 from the record in this band (see harmonicTemplate).
%                With neither F0 nor F0Band the template step is skipped and x is shaped as it is.
%   K          - Name-value. Number of harmonics stripped. REQUIRED whenever the template step runs (F0 or F0Band given); there is
%                no default, and a missing K is an error ("shapeNoiseWrapper:KRequired"). Choose it for your record with harmonicContent
%                (README, "Choosing K for your data"). Not used when there is no template step.
%   Window     - Name-value. Template window in cycles of f0 (default 4).
%   Alpha      - Name-value. Target exponent per channel (default [] = that of the residual; see shapeNoise).
%   AlphaBand  - Name-value. [fLow fHigh] for the exponent, or "kinematic" (default). See shapeNoise.
%   ShapeArgs  - Name-value. Cell array of further shapeNoise name-value pairs, e.g. {"MinLength", 64}.
%
% OUTPUTS:
%   n    - The surrogate (info.shape.L samples per channel, see shapeNoise).
%   info - Struct: template (harmonicTemplate info, or []), residual (what shapeNoise was given), shape (shapeNoise info).
%
% The measured range in the shapeNoise help (RANGE) applies; outside it this is exploratory.
%
% ERRORS (fail loud): "shapeNoiseWrapper:KRequired" when F0 or F0Band is given without K; otherwise those of harmonicTemplate and shapeNoise.
%
% See also: harmonicContent, harmonicTemplate, shapeNoise, cetNoiseWrapper

arguments
    x double {mustBeReal, mustBeNonempty}
    fs (1, 1) double {mustBePositive}
    options.F0 double {mustBeScalarOrEmpty} = []
    options.F0Band (1, 2) double = [NaN NaN]
    options.K double {mustBeScalarOrEmpty, mustBeInteger, mustBePositive} = []
    options.Window (1, 1) double = 4
    options.Alpha double = []
    options.AlphaBand = "kinematic"
    options.ShapeArgs (1, :) cell = {}
end
useTemplate = ~isempty(options.F0) || all(~isnan(options.F0Band));
if useTemplate && isempty(options.K)
    error("shapeNoiseWrapper:KRequired", "K, the number of harmonics to strip, is required when F0 or F0Band is given; there is no default. " + ...
        "Choose it for this record with harmonicContent(x, fs, f0) and pass K=..., e.g. shapeNoiseWrapper(x, fs, F0=f0, K=2).");
end
if useTemplate
    [resid, ~, tinfo] = harmonicTemplate(x, fs, options.F0, options.K, F0Band=options.F0Band, Window=options.Window);
else
    resid = x;  tinfo = [];
end
[n, sinfo] = shapeNoise(resid, fs, options.Alpha, "AlphaBand", options.AlphaBand, options.ShapeArgs{:});
info = struct("template", tinfo, "residual", resid, "shape", sinfo);
end
