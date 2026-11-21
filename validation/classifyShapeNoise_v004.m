function R = classifyShapeNoise_v004(x, fs, band)
%CLASSIFYSHAPENOISE_V004 One default shapeNoise call, with its outcome returned as data.
%   R.cls is "ok", "noSegment", "segmentAlpha", "noAccepted" or "other". R.id is the error identifier ("" when none); for segmentAlpha it
%   also carries the identifier of the cause, as "shapeNoise:segmentAlpha<cause". R.s and R.info are empty when the call failed.
R = struct('cls', "ok", 'id', "", 's', [], 'info', []);
try
    [R.s, R.info] = shapeNoise(x, fs, [], AlphaBand=band);
    if ~all(R.info.accepted), R.cls = "noAccepted";  R.s = []; end
catch err
    R.id = string(err.identifier);
    switch R.id
        case "shapeNoise:noSegment",    R.cls = "noSegment";
        case "shapeNoise:segmentAlpha", R.cls = "segmentAlpha";
            if ~isempty(err.cause), R.id = R.id + "<" + string(err.cause{1}.identifier); end
        otherwise,                      R.cls = "other";
    end
end
end
