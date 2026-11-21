%[text] # Comparing Noise Generation Options
%[text] This toolbox offers several ways to get 1/f^alpha-ish noise, each suited to a different situation. This example generates the same target alpha via each route and compares the recovered alpha and spectral shape, to help you choose.
%[text] ## The options
%[text] - `noiseXu` -- parametric generator. You specify alpha directly; it synthesises a signal to match. Use this when you want noise with a  **known** ground-truth alpha (simulation, testing, stress-testing a pipeline).
%[text] - `shapeNoise` (default, end-matched) -- takes an  **empirical** residual, finds its longest segment whose ends meet, and puts random phases on that segment's amplitude spectrum, keeping a draw only if its exponent matches. Use this when you have real data and want a surrogate that matches its spectrum and derivative noise without reusing its exact realisation. The output is shorter than the input (`info.L`).
%[text] - `shapeNoise` with `Method="legacy"` -- the 1.x construction: the whole record's amplitude spectrum on a phase carrier (`PhaseSource` `"xu"`, `"white"`, `"pinknoise"`, `"dsp"`, or a handle). Kept to reproduce earlier results; at steep alpha a share of its surrogates come out shallow, and it warns.
%[text] - `external/fftnoise/fftnoise.m` -- a **phase-surrogate** generator (Aslak Grinsted, 2011): randomises the phases of a whole record's FFT. Different from  `noiseXu` in kind -- it needs a real FFT to start from. Applied to a whole steep record it inherits the end-mismatch leakage; see `examples/CompareFftnoiseVsShapedXu.m`. \
%[text] ## Generating the same target alpha several ways
fs = 100;
N = 3000;
targetAlpha = 2.5;
fLow = 1.0;
fHigh = 25.0;

rng(10);
xuDirect = noiseXu(N, targetAlpha, 3.0, fs);

% For shapeNoise, first need an empirical-ish residual with roughly the
% target spectral shape to draw the magnitude from -- use a noiseXu
% draw for this demo (in real use this would be your actual data).
rng(11);
empiricalStandIn = noiseXu(N, targetAlpha, 3.0, fs);

rng(12);
[shapedEM, info] = shapeNoise(empiricalStandIn, fs, []);
fprintf("End-matched segment: L = %d of %d samples\n", info.L, N);
rng(13);
shapedLegacy = shapeNoise(empiricalStandIn, fs, targetAlpha, Method="legacy", PhaseSource="xu");
%%
%[text] ## Recovering alpha from each
methodNames = ["noiseXu (direct)"; "shapeNoise (end-matched)"; "shapeNoise (legacy, xu carrier)"];
signals = {xuDirect, shapedEM, shapedLegacy};
recoveredAlpha = zeros(numel(signals), 1);
for i = 1:numel(signals)
    recoveredAlpha(i) = estimateIRASA(signals{i}, fs, fLow, fHigh);
end

comparison = table(methodNames, recoveredAlpha, ...
    VariableNames=["Method", "alpha_recovered"]);
disp(comparison)
fprintf("Target alpha: %.2f\n", targetAlpha);
%%
%[text] ## Spectral comparison
figure;
tiledlayout(1, 2);

nexttile;
[pxxXu, fXu] = pmtm(xuDirect, 4, [], fs);
[pxxShapedEM, fShapedEM] = pmtm(shapedEM, 4, [], fs);
[pxxShapedLegacy, fShapedLegacy] = pmtm(shapedLegacy, 4, [], fs);
loglog(fXu, pxxXu, DisplayName="noiseXu (direct)");
hold on;
loglog(fShapedEM, pxxShapedEM, DisplayName="shapeNoise (end-matched)");
loglog(fShapedLegacy, pxxShapedLegacy, DisplayName="shapeNoise (legacy)");
hold off;
xlabel("Frequency (Hz)");
ylabel("PSD");
title("Power spectra (log-log)");
legend(Location="southwest");
grid on;

nexttile;
bar(categorical(methodNames), recoveredAlpha);
yline(targetAlpha, "r--", "target alpha", LineWidth=1.5);
ylabel("Recovered alpha");
title("Alpha recovery by method");
xtickangle(20);
%%
%[text] At alpha 2.5 all three cluster near the target: the end-mismatch problem that sinks the legacy construction is rare below alpha 3. Run the comparison at alpha 4 or 5, or over many seeds, to see the legacy tail (`validation/checkShapeNoiseEndMatched_v001.m` does this).
%[text] ## When would the construction matter beyond alpha?
%[text] If your downstream analysis is sensitive to something **other** than the spectral exponent -- derivative noise (velocity, acceleration, curvature), higher-order statistics, waveform shape -- check the surrogate on that quantity, not on alpha alone. On real movement residuals the legacy construction matched alpha but left velocity and curvature noise up to 30% short for some datasets; the end-matched surrogate matched both.
%[text] ## Cross-references
%[text] - `noiseXu.m`, `shapeNoise.m` docstrings
%[text] - `doc/GettingStarted.mlx` for a broader tour of the toolbox \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
