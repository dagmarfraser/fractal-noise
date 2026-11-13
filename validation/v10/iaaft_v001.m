function y = iaaft_v001(x, nIter)
%IAAFT_V001 Iterative amplitude-adjusted Fourier surrogate (Schreiber & Schmitz 1996).
%   Keeps the amplitude spectrum and the value distribution of x; destroys nonlinear structure
%   (including multifractality). Uses the global RNG.
arguments
    x (:,1) double
    nIter (1,1) double {mustBePositive} = 100
end
A = abs(fft(x));  xs = sort(x);  N = numel(x);
y = x(randperm(N));
for it = 1:nIter
    Y = fft(y);  y = real(ifft(A .* exp(1i * angle(Y))));
    [~, r] = sort(y);  y(r) = xs;
end
end
