import numpy as np

def rule_defn_read(name):
    """Load a substitution rule from a .npy file and return as a NumPy array."""
    return np.load(name + ".npy")

def sub_on_patch(sub_deg, patch, tile_rule):
    """Apply one step of the substitution rule to a patch."""
    wid = [len(patch), len(patch[0])]
    new_wid = [ w*sub_deg for w in wid ]
    arr = np.full((new_wid[0], new_wid[1]), -1, dtype=int)
    patch = np.array(patch, dtype=int)
    for i in range(wid[0]):
        for j in range(wid[1]):
            for li in range(sub_deg):
                for ki in range(sub_deg):
                    arr[i * sub_deg + li, j * sub_deg + ki] = tile_rule[patch[i, j], li, ki]
    return arr


def iterate_sub(iter_num, sub_deg, patch, tile_rule):
    """Apply the substitution rule iter_num times. Returns original patch if iter_num=0."""
    for _ in range(iter_num):
        patch = sub_on_patch(sub_deg, patch, tile_rule)
    return patch
