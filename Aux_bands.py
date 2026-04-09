import numpy as np
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from numpy import linalg as LA
from Aux_matrix import sample_band_edges
import os

#Adds the proper suffix when safing a file
def get_ordinal_suffix(num: int) -> str:
    """Return the ordinal suffix ('st', 'nd', 'rd', 'th') for a non-negative integer."""
    if not isinstance(num, int) or num < 0:
        raise ValueError(f"Expected a non-negative integer, got {num!r}")
    if 10 <= num <= 20:
        return "th"
    last_digit = num % 10
    if last_digit == 1: return "st"
    if last_digit == 2: return "nd"
    if last_digit == 3: return "rd"
    return "th"

def _band_filename(itera, eig_num, val_temp):
    """Construct the canonical filename for a saved band computation."""
    rounded = [round(v, 4) for v in val_temp]
    return f"Finished computations/{itera}{get_ordinal_suffix(itera)}, {eig_num} comput for values {rounded}.npy"

def save_bands(itera, res, eig_num, patch, val_temp):
    """
    Compute and save band data for the given parameters.
    Saves a file .npy binary file with the full array
    """
    save_arr = sample_band_edges(res, patch, val_temp)
    base_path = _band_filename(itera, eig_num, val_temp)  # includes .npy extension

    # Save binary
    np.save(base_path, save_arr)

def read_bands(itera, eig_num, val_temp):
    """Load and return saved band data for the given parameters."""
    path = _band_filename(itera, eig_num, val_temp)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Band file not found: {path}")
    return np.load(path)



##################### plotting function #####################
def print_several_bands(
        init_itera,
        fin_itera,
        val_temp,
        init_size,
        vert_dist=0.5,
        plot_tiny=False,
        tiny_tol=1e-8,
        tiny_height_factor=0.15,
        color='black',
        linewidth=0.9,
        fig_width=8.0,
        save_path=None):
    """
    Plot periodic approximation spectra for several substitution iterations.

    Parameters
    ----------
    init_itera : int
        First iteration to plot.

    fin_itera : int
        Last iteration to plot.

    val_temp : list or dict
        Potential values used to label the data files.

    init_size : int
        Number of eigenvalues at iteration 0.
        At iteration it, the number of bands is
            size = init_size * 4**it.

    vert_dist : float, optional
        Visual vertical spacing between consecutive iterations.
        IMPORTANT:
        This does NOT change the y-data coordinates.
        The iterations are still plotted at y = 0,1,2,...
        Instead, vert_dist controls the figure height per iteration.

    plot_tiny : bool, optional
        If True, numerically tiny bands are shown as vertical markers.
        If False, such bands are omitted.

    tiny_tol : float, optional
        A band is treated as tiny if right - left < tiny_tol.

    tiny_height_factor : float, optional
        Height of the vertical tiny-band marker relative to one iteration step.

    color : str, optional
        Plot color for all bands.

    linewidth : float, optional
        Line thickness for both usual and tiny bands.

    fig_width : float, optional
        Width of the figure in inches.

    save_path : str or None, optional
        If not None, the figure is saved to this path.

    Notes
    -----
    The y-axis always uses the natural iteration numbers.
    The parameter vert_dist changes only the *visual* distance between them
    by changing the figure height.
    """

    n_levels = fin_itera - init_itera + 1

    # Figure height is controlled by vert_dist:
    # smaller vert_dist -> more compact plot
    fig_height = max(2.0, vert_dist * n_levels + 1.0)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    ax.set_xlabel('Spectrum')
    ax.set_ylabel('Iteration number n')
    #ax.set_title(f"Periodic approximation spectra, iterations {init_itera}-{fin_itera}")

    horizontal_segments = []
    vertical_segments = []

    for it in range(init_itera, fin_itera + 1):
        y = it
        size = init_size * (4 ** it)

        mat = read_bands(it, 0, val_temp)

        for i in range(size):
            left = min(mat[1][i], mat[0][i])
            right = max(mat[1][i], mat[0][i])
            width = right - left

            if width >= tiny_tol:
                horizontal_segments.append([(left, y), (right, y)])
            else:
                if plot_tiny:
                    center = 0.5 * (left + right)
                    half_height = tiny_height_factor
                    vertical_segments.append([
                        (center, y - half_height),
                        (center, y + half_height)
                    ])

    if horizontal_segments:
        h_collection = LineCollection(
            horizontal_segments,
            colors=color,
            linewidths=linewidth
        )
        ax.add_collection(h_collection)

    if vertical_segments:
        v_collection = LineCollection(
            vertical_segments,
            colors=color,
            linewidths=linewidth
        )
        ax.add_collection(v_collection)

    # Keep natural iteration labels on the y-axis
    ax.set_yticks(np.arange(init_itera, fin_itera + 1, 1))

    # Add a small vertical margin
    ax.set_ylim(init_itera - 0.5, fin_itera + 0.5)

    ax.autoscale_view(scalex=True, scaley=False)

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, bbox_inches='tight')
        print(f"Figure saved to {save_path}")

    plt.show()
    
