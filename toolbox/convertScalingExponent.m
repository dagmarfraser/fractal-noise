function y = convertScalingExponent(x, from, to, options)
%CONVERTSCALINGEXPONENT Convert between the spectral exponent alpha, Hurst H, long-memory d, DFA exponent and fractal dimension D.
%
%   y = convertScalingExponent(x, from, to)
%   y = convertScalingExponent(x, from, to, Class="fGn")      % needed whenever "H" or "D" is involved
%
%   For a Gaussian, self-affine series with one scaling regime, these are one quantity in different units:
%
%     alpha  spectral exponent, S(f) ~ f^(-alpha) at low frequency (what estimateWhittleAlpha and estimateIRASA return)
%     d      long-memory (fractional differencing) parameter: d = alpha / 2 (noiseXu builds its noise from d)
%     DFA    detrended fluctuation exponent a, F(s) ~ s^a: alpha = 2a - 1 (any DFA order that removes the trend)
%     H      Hurst exponent, in (0, 1). Its meaning depends on the class of the series:
%              Class "fGn" (stationary noise, -1 < alpha < 1): alpha = 2H - 1
%              Class "fBm" (random-walk-like, 1 < alpha < 3):  alpha = 2H + 1
%     D      fractal (Hausdorff) dimension of a graph, in (1, 2), D = 2 - H:
%              Class "fBm": D of the series itself, alpha = 5 - 2D
%              Class "fGn": D of the series' cumulative sum (a noise's own graph is space-filling), alpha = 3 - 2D
%
%   alpha, d and DFA convert without a class. H and D need Class, because the same H means alpha = 2H - 1 for a noise and
%   2H + 1 for its running sum; a value outside the class's range is an error, not a clipped answer. Odd-integer alpha
%   (-1, 1, 3) is a boundary of both classes (H = 0 or 1) and has no H or D.
%
%   CAVEATS. These are identities for self-affine Gaussian processes. In data, D is usually measured at small scales (roughness)
%   and H, d, DFA and alpha at large scales (memory); they agree only if the series scales the same way across both, and in
%   general they need not (Gneiting & Schlather 2004). Converting a published value assumes the authors' estimator measured the
%   same thing over a comparable range of scales; R/S, for example, is biased upwards for anti-persistent series. The
%   conversion carries no uncertainty: convert an interval's ends, not its centre alone.
%
% INPUTS:
%   x       - Numeric array of exponents (finite, real). The output has the same size.
%   from    - "alpha", "d", "DFA", "H" or "D".
%   to      - "alpha", "d", "DFA", "H" or "D".
%   Class   - Name-value. "fGn" or "fBm"; required when from or to is "H" or "D".
%
% OUTPUTS:
%   y - Converted exponents, same size as x.
%
% ERRORS (fail loud): "convertScalingExponent:classRequired" (H or D without Class), "convertScalingExponent:outOfRange" (a value
%   outside the class's range, H outside (0, 1), D outside (1, 2), or an odd-integer boundary).
%
% EXAMPLE:
%   convertScalingExponent(0.39, "alpha", "H", Class="fGn")    % 0.695: monthly CET anomalies (README)
%   convertScalingExponent(0.72, "H", "alpha", Class="fGn")    % 0.44
%   convertScalingExponent(0.9, "DFA", "alpha")                % 0.8
%   convertScalingExponent(1.3, "D", "H", Class="fGn")         % 0.7
%
% REFERENCES: Mandelbrot, B. B., & Van Ness, J. W. (1968). Fractional Brownian motions, fractional noises and applications.
%   SIAM Review, 10(4), 422-437. Gneiting, T., & Schlather, M. (2004). Stochastic models that separate fractal dimension and the
%   Hurst effect. SIAM Review, 46(2), 269-282. Kiyono, K. (2015). Establishing a direct connection between detrended
%   fluctuation analysis and Fourier analysis. Physical Review E, 92, 042925.
%
% Toolbox: fractalnoise
%
% See also: estimateWhittleAlpha, estimateIRASA, noiseXu

arguments
    x {mustBeNumeric, mustBeReal, mustBeFinite}
    from (1, 1) string {mustBeMember(from, ["alpha", "d", "DFA", "H", "D"])}
    to (1, 1) string {mustBeMember(to, ["alpha", "d", "DFA", "H", "D"])}
    options.Class (1, 1) string {mustBeMember(options.Class, ["", "fGn", "fBm"])} = ""
end
x = double(x);
needsClass = any(ismember([from, to], ["H", "D"]));
if needsClass && options.Class == ""
    error("convertScalingExponent:classRequired", ['Converting %s to %s needs Class="fGn" (stationary noise, -1 < alpha < 1) or ' ...
        'Class="fBm" (random-walk-like, 1 < alpha < 3): the same H means alpha = 2H - 1 for a noise and 2H + 1 for its running sum.'], from, to);
end
fBm = options.Class == "fBm";

% To alpha
switch from
    case "alpha", a = x;
    case "d",     a = 2 * x;
    case "DFA",   a = 2 * x - 1;
    case "H"
        checkRange(x, 0, 1, "H");
        a = 2 * x - 1 + 2 * fBm;
    case "D"
        checkRange(x, 1, 2, "D");
        a = 3 - 2 * x + 2 * fBm;
end

% From alpha
if needsClass
    lo = -1 + 2 * fBm;  hi = lo + 2;
    if any(a <= lo | a >= hi, "all")
        bad = a(find(a <= lo | a >= hi, 1));
        error("convertScalingExponent:outOfRange", ['alpha = %g is outside the %s range (%g, %g); H and D exist only strictly inside it. ' ...
            'Check the class (alpha above 1 is fBm-like) or the value.'], bad, options.Class, lo, hi);
    end
end
switch to
    case "alpha", y = a;
    case "d",     y = a / 2;
    case "DFA",   y = (a + 1) / 2;
    case "H",     y = (a + 1 - 2 * fBm) / 2;
    case "D",     y = (3 - a + 2 * fBm) / 2;
end
end

function checkRange(v, lo, hi, name)
if any(v <= lo | v >= hi, "all")
    bad = v(find(v <= lo | v >= hi, 1));
    error("convertScalingExponent:outOfRange", "%s = %g is outside (%g, %g).", name, bad, lo, hi);
end
end
