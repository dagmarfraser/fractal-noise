function out = buildToolbox_v001(version)
%BUILDTOOLBOX_V001 Package toolbox/ as release/Fractal Noise Toolbox.mltbx (saved build script).
%   out = buildToolbox_v001("1.2.0") rebuilds the .mltbx with the same GUID as v1.1.0, minimum
%   release R2025a, and the current description. Overwrites release/Fractal Noise Toolbox.mltbx
%   (git tags/releases hold old versions). Attach the output to the GitHub Release so File Exchange
%   picks up this exact package.
arguments
    version (1,1) string = "1.2.0"
end
root = fileparts(fileparts(mfilename('fullpath')));
tb   = fullfile(root, 'toolbox');
o = matlab.addons.toolbox.ToolboxOptions(tb, "f070d0f9-6828-49e2-9178-1861e71f6fe2");
o.ToolboxName    = "Fractal Noise Toolbox";
o.ToolboxVersion = version;
o.Summary        = "Generate, shape, estimate, and diagnose 1/f^alpha (fractal) noise in MATLAB.";
o.Description    = "Small, standalone toolbox for 1/f^alpha (fractal/colored) noise: noiseXu " + ...
    "(Xu 2019 GGM generation; Phi = 1 for a break-free power law), shapeNoise (impose an empirical " + ...
    "spectrum onto a chosen phase carrier), estimateWhittleAlpha (prewhitened local Whittle exponent, " + ...
    "reliable for steep spectra), estimateIRASA (IRASA, Wen & Liu 2016, with optional Welch PSD and " + ...
    "prewhitening), and checkSpectralKnee (Lorentzian-vs-fixed knee diagnostic).";
o.AuthorName         = "Dagmar Scott Fraser";
o.AuthorEmail        = "d.s.fraser@bham.ac.uk";
o.AuthorCompany      = "University of Birmingham";
o.MinimumMatlabRelease = "R2025a";
o.OutputFile         = fullfile(root, 'release', 'Fractal Noise Toolbox.mltbx');
matlab.addons.toolbox.packageToolbox(o);
out = o.OutputFile;
fprintf('Built %s (v%s, min %s, %d files)\n', out, version, o.MinimumMatlabRelease, numel(o.ToolboxFiles));
end
