import math
import cmath
import numpy as np
from fractions import Fraction
from numpy import linalg as LA

# Function returns the current potential value
def sub_diag(index, val):
    """Return the diagonal value at the given index (0–3), or 0 if out of range."""
    if not 0 <= index < len(val):
        raise IndexError(f"index {index} out of range for val of length {len(val)}")
    return val[index]


def nested_phase_points(res):
    """
    Build a nested set of phase points in [0, pi].

    For a given resolution parameter res we take all rational numbers
    a/q with

        0 <= a <= q
        1 <= q <= res

    and map them to

        theta = pi * a/q

    This produces a nested set of phase points:
    increasing res only adds new points.

    Returns
    -------
    phases : np.ndarray
        Sorted unique phases in [0, pi].
    """

    points = set()

    for q in range(1, res + 1):
        for a in range(q + 1):
            points.add(Fraction(a, q))

    points = sorted(points)

    phases = np.array([np.pi * float(t) for t in points])

    return phases
    
def op_mat_no_phase(patch, val):
    """
    Build the real-space part of the Schroedinger operator matrix for a
    height x wid grid. Sets nearest-neighbour hopping (=1) and diagonal
    values from val.
    Open boundary conditions — edges do not wrap.
    """
    height = len(patch)
    wid = len(patch[0])
    size = wid*height
    mat = np.zeros((size, size), dtype=complex)

    for i in range(height):
        for j in range(wid):
            idx = i * wid + j
            if i < height - 1: mat[idx, idx + wid] = 1       # hop down
            if i > 0:       mat[idx, idx - wid] = 1       # hop up
            if j < wid - 1: mat[idx, idx + 1]  = 1        # hop right
            if j > 0:       mat[idx, idx - 1]  = 1        # hop left
            mat[idx, idx] = sub_diag(patch[i][j], val)

    return mat


def sample_band_edges(res, patch, val):
    """
    Compute sampled band edges (min/max over phases) for all eigenvalues.

    Parameters
    ----------
    res : int
        Resolution parameter for the phase sampling.
    patch : list or np.ndarray
        2D substitution patch.
    val : dict or array-like
        Potential values for symbols in the patch.
    """

    height = len(patch)
    wid = len(patch[0])
    size = height * wid

    phases = nested_phase_points(res)
    n_pts = len(phases)

    # Base matrix without phase-dependent wrap-around terms
    mat1 = op_mat_no_phase(patch, val)

    # Precompute exponentials
    exp_pos = np.exp(1j * phases)
    exp_neg = np.exp(-1j * phases)

    # Precompute boundary index pairs
    vertical_pairs = [(j, (height - 1) * wid + j) for j in range(wid)]
    horizontal_pairs = [(i * wid, i * wid + (wid - 1)) for i in range(height)]

    band_max = None
    band_min = None

    for ki in range(n_pts):
        e1_pos = exp_pos[ki]
        e1_neg = exp_neg[ki]

        for li in range(n_pts):
            e2_pos = exp_pos[li]
            e2_neg = exp_neg[li]

            mat = mat1.copy()

            # Horizontal wrap
            for left, right in horizontal_pairs:
                mat[left, right] += e1_neg
                mat[right, left] += e1_pos

            # Vertical wrap
            for top, bottom in vertical_pairs:
                mat[top, bottom] += e2_neg
                mat[bottom, top] += e2_pos

            vect = LA.eigvalsh(mat)

            if band_max is None:
                band_max = vect.copy()
                band_min = vect.copy()
            else:
                band_max = np.maximum(band_max, vect)
                band_min = np.minimum(band_min, vect)

    return np.vstack((band_max, band_min))
