"""Python port of V0 (exactReference/*_v001.m): exact fGn by Davies-Harte, k-fold integration, analytic spectrum.
alpha = (2H - 1) + 2k, 0 < H < 1; odd-integer alpha raises (no exact representative), never nudged."""
import numpy as np

def fgn_autocov(k, H):
    k = np.abs(np.asarray(k, float))
    return 0.5 * ((k + 1) ** (2 * H) - 2 * k ** (2 * H) + np.abs(k - 1) ** (2 * H))

def fgn_davies_harte(N, H, nrep, rng, eig_tol=1e-10):
    g = fgn_autocov(np.arange(N), H)
    c = np.concatenate([g, g[-2:0:-1]]); M = c.size
    lam = np.real(np.fft.fft(c))
    if lam.min() < -eig_tol * lam.max():
        raise ValueError(f"fgnDaviesHarte:negEig min eigenvalue {lam.min():.3g} (H={H}, N={N})")
    lam = np.clip(lam, 0, None)
    nc = int(np.ceil(nrep / 2)); out = []
    for _ in range(nc):
        W = np.fft.fft(np.sqrt(lam / M) * (rng.standard_normal(M) + 1j * rng.standard_normal(M)))
        out += [W.real[:N], W.imag[:N]]
    return np.column_stack(out[:nrep])

def alpha_to_kH(alpha):
    k = int(np.floor((alpha + 1) / 2)); H = (alpha + 1 - 2 * k) / 2
    if H <= 1e-12 or H >= 1 - 1e-12:
        raise ValueError(f"exactRefSeries:oddAlpha alpha={alpha} is an odd-integer boundary (H={H})")
    return k, H

def fgn_spectrum(lam, H, J=200):
    lam = np.atleast_1d(lam)[:, None]; j = np.arange(-J, J + 1)[None, :]
    from math import gamma, sin, pi
    cH = gamma(2 * H + 1) * sin(pi * H) / (2 * pi)
    s = np.sum(np.abs(lam + 2 * pi * j) ** (-2 * H - 1), 1)
    l = lam[:, 0]
    tail = ((2 * pi * (J + .5) + l) ** (-2 * H) + (2 * pi * (J + .5) - l) ** (-2 * H)) / (2 * pi * 2 * H)
    return 2 * cH * (1 - np.cos(l)) * (s + tail)

def exact_ref_series(N, alpha, nrep, rng):
    k, H = alpha_to_kH(alpha)
    X = fgn_davies_harte(N, H, nrep, rng)
    for _ in range(k): X = np.cumsum(X, 0)
    S = lambda lam: fgn_spectrum(lam, H) / (2 - 2 * np.cos(lam)) ** k
    return X, dict(k=k, H=H, S=S)
