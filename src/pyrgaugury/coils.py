"""Coiled-coil prediction: a numpy port of nCoils (Lupas et al. 1991).

Reproduces ``ncoils -f`` as RGAugury ran it (window 21, unweighted, MTIDK
"new.mat" matrix, P >= 0.5): every window gets the geometric mean of its
heptad-position weights, each residue takes the best score of any window
covering it, and that score is converted to a coil probability P.
"""

from __future__ import annotations

import numpy as np

# new.mat from the nCoils distribution: per-residue weights at heptad a..g
_MATRIX = {
    "L": (2.998, 0.269, 0.367, 3.852, 0.510, 0.514, 0.562),
    "I": (2.408, 0.261, 0.345, 0.931, 0.402, 0.440, 0.289),
    "V": (1.525, 0.479, 0.350, 0.887, 0.286, 0.350, 0.362),
    "M": (2.161, 0.605, 0.442, 1.441, 0.607, 0.457, 0.570),
    "F": (0.490, 0.075, 0.391, 0.639, 0.125, 0.081, 0.038),
    "Y": (1.319, 0.064, 0.081, 1.526, 0.204, 0.118, 0.096),
    "G": (0.084, 0.215, 0.432, 0.111, 0.153, 0.367, 0.125),
    "A": (1.283, 1.364, 1.077, 2.219, 0.490, 1.265, 0.903),
    "K": (1.233, 2.194, 1.817, 0.611, 2.095, 1.686, 2.027),
    "R": (1.014, 1.476, 1.771, 0.114, 1.667, 2.006, 1.844),
    "H": (0.590, 0.646, 0.584, 0.842, 0.307, 0.611, 0.396),
    "E": (0.281, 3.351, 2.998, 0.789, 4.868, 2.735, 3.812),
    "D": (0.068, 2.103, 1.646, 0.182, 0.664, 1.581, 1.401),
    "Q": (0.311, 2.290, 2.330, 0.811, 2.596, 2.155, 2.585),
    "N": (1.231, 1.683, 2.157, 0.197, 1.653, 2.430, 2.065),
    "S": (0.332, 0.753, 0.930, 0.424, 0.734, 0.801, 0.518),
    "T": (0.197, 0.543, 0.647, 0.680, 0.905, 0.643, 0.808),
    "C": (0.918, 0.002, 0.385, 0.440, 0.138, 0.432, 0.079),
    "W": (0.066, 0.064, 0.065, 0.747, 0.006, 0.115, 0.014),
    "P": (0.004, 0.108, 0.018, 0.006, 0.010, 0.004, 0.007),
}
# unweighted fitting parameters for each window: m_cc, sd_cc, m_g, sd_g, sc
_FIT = {14: (1.82, 0.28, 0.95, 0.26, 20), 21: (1.74, 0.23, 0.86, 0.21, 25),
        28: (1.69, 0.18, 0.80, 0.18, 30)}

_LOGW = np.zeros((256, 7))
_VALID = np.zeros(256, dtype=bool)
for _aa, _w in _MATRIX.items():
    _LOGW[ord(_aa)] = np.log(_w)
    _VALID[ord(_aa)] = True


def coil_probabilities(seq: str, window: int = 21) -> np.ndarray:
    """Per-residue coiled-coil probability P, as computed by nCoils."""
    codes = np.frombuffer(seq.upper().encode("ascii", "replace"), dtype=np.uint8)
    n = len(codes)
    valid = _VALID[codes]
    score = np.zeros(n)
    n_windows = n - window + 1
    if n_windows > 0:
        window_score = np.zeros(n_windows)
        positions = np.arange(n)
        valid_cum = np.concatenate(([0], np.cumsum(valid)))
        counts = valid_cum[window:] - valid_cum[:-window]
        for frame in range(7):
            # for windows starting at i = frame (mod 7), residue r sits at
            # heptad position (r - frame) mod 7
            logw = _LOGW[codes, (positions - frame) % 7] * valid
            cum = np.concatenate(([0.0], np.cumsum(logw)))
            starts = np.arange(frame, n_windows, 7)
            sums = cum[starts + window] - cum[starts]
            c = counts[starts]
            with np.errstate(divide="ignore", invalid="ignore"):
                window_score[starts] = np.where(c > 0, np.exp(sums / np.maximum(c, 1)), 0.0)
        # each residue takes the best score among windows covering it
        padded = np.concatenate((np.full(window - 1, -np.inf), window_score,
                                 np.full(window - 1, -np.inf)))
        best = np.lib.stride_tricks.sliding_window_view(padded, window).max(axis=1)
        score = np.where(valid, np.maximum(best, 0.0), 0.0)
    m_cc, sd_cc, m_g, sd_g, sc = _FIT[window]
    g_cc = np.exp(-0.5 * ((score - m_cc) / sd_cc) ** 2) / sd_cc
    g_g = np.exp(-0.5 * ((score - m_g) / sd_g) ** 2) / sd_g
    return g_cc / (sc * g_g + g_cc)


def coil_segments(seq: str, window: int = 21, min_p: float = 0.5) -> list[tuple[int, int]]:
    """1-based inclusive (start, end) runs of residues with P >= ``min_p``."""
    is_coil = coil_probabilities(seq, window) >= min_p
    edges = np.flatnonzero(np.diff(np.concatenate(([0], is_coil.astype(np.int8), [0]))))
    return [(int(s) + 1, int(e)) for s, e in zip(edges[::2], edges[1::2])]
