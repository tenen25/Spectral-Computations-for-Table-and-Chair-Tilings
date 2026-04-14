# Spectral Computation Code

> Computes sampled spectral bands of discrete Schrödinger operators of periodic approximants associated with substitution-generated patterns.

---

## Table of Contents

- [Requirements](#requirements)
- [Setup](#setup)
- [File Structure](#file-structure)
- [Inputs](#inputs)
- [Workflow](#workflow)
- [Key Concepts](#key-concepts)

---

## Requirements

- Python 3.x
- The following libraries:

```bash
pip install numpy matplotlib
```

> `fractions` is part of the Python standard library — no installation needed.

---

## Setup

Before running, create the output folder manually:

```bash
mkdir "Finished computations"
```

All computation results will be saved there as `.npy` band data files and `.png` plot images, named using ordinal labels (e.g. `1st`, `2nd`, `3rd`).

---

## File Structure

```
.
├── main.py                    # Driver script — runs the full pipeline
├── Sub_Rule.py                # Substitution dynamics
├── Aux_matrix.py              # Matrix construction and eigenvalue routines
├── Aux_bands.py               # Band storage and plotting routines
└── Finished computations/     # Output folder (create manually)
```

### Module Dependencies

```
main.py ──► Aux_bands.py ──► Aux_matrix.py
main.py ──► Sub_Rule.py
```

| File | Key Functions |
|---|---|
| `Sub_Rule.py` | `rule_defn_read`, `sub_on_patch`, `iterate_sub` |
| `Aux_matrix.py` | `sub_diag`, `nested_phase_points`, `op_mat_no_phase`, `sample_band_edges` |
| `Aux_bands.py` | `save_bands`, `read_bands`, `print_several_bands`, `get_ordinal_suffix` |

---

## Inputs

The following parameters are set directly in `main.py` before running:

| Parameter | Type | Description |
|---|---|---|
| `RULE_NAME` | `str` | Name of the substitution rule to load (either `"Table_Rule"` or `"Chair_Rule"`) |
| `SUB_DEG` | `int` | Substitution degree — typically `2`, each tile expands into a `SUB_DEG × SUB_DEG` block |
| `start_tile` | `list[list[int]]` | 2D seed tile for the first substitution iteration |
| `max_itera` | `int` | Number of iteration steps |
| `res` | `int` | Phase resolution — controls density of Bloch phase sampling (`q ≤ res`) |
| `val` | `list[float]` | On-site potential values assigned to each tile symbol |
| `init_itera` | `int` | Plot range start — `0` includes the seed iteration |
| `fin_itera` | `int` | Plot range end — inclusive, must be `≤ max_itera` |

### Example of inputs

```python
RULE_NAME  = "Chair_Rule"
SUB_DEG    = 2                          # each tile expands into a 2×2 block
start_tile = [[1,2,3], [0,1,2]]         # a 2×3 seed patch
max_itera  = 4                          # compute for `0 ≤ it ≤ max_itera`
res        = 9                          # samples phases with denominator up to 9
val        = [0.0, 9.0, 18.0, 27.0]     # one potential per tile symbol
init_itera = 0                          # plot from the seed iteration
fin_itera  = 4                          # plot up to and including max_itera
```

> **Note:** The length of `val` must match the number of distinct tile symbols in the substitution rule.

---

## Workflow

When `main.py` is executed, it runs the following steps:

1. **Read substitution rule**
   ```python
   tile_rule = rule_defn_read(RULE_NAME)
   ```

2. **Generate substituted patch** — for each `it` with `0 ≤ it ≤ max_itera`, produces a patch of size `(len(start_tile) × 2^it) × (len(start_tile[0]) × 2^it)`
   ```python
   patch = iterate_sub(it, SUB_DEG, start_tile, tile_rule)
   ```

3. **Compute and save sampled bands** for each `it` with `0 ≤ it ≤ max_itera`
   ```python
   save_bands(it, res, patch, val)
   ```

4. **Plot the saved bands** over the iteration range `[init_itera, fin_itera]`
   ```python
   print_several_bands(init_itera, fin_itera, ...)
   ```

### Core Call Chain

```
main → save_bands → sample_band_edges(res, patch, val_temp) → op_mat_no_phase → eigvalsh
main → print_several_bands → read_bands
```

---

## Key Concepts

### Patch Numbering
The 2D patch is stored as `patch[i][j]` and flattened **row-major**:

```
(i, j)  →  k = i × wid + j
```

For a 4×4 patch, the enumeration is as follows:

```
 0   1   2   3
 4   5   6   7
 8   9  10  11
12  13  14  15
```

### Phase Sampling
Bloch phase pairs `(θ₁, θ₂)` are sampled, where each component independently
takes values of the form:

    θ = π · a/q,   where gcd(a, q) = 1,   q ≤ res

Duplicates are removed from each component set automatically using Python's
`fractions.Fraction` and `set`.

### Sampled Bands
For each eigenvalue branch `k`, the code records the sampled min and max over all phase pairs `(θ₁, θ₂)`:

```
band_min[k] = min λₖ(θ₁, θ₂)
band_max[k] = max λₖ(θ₁, θ₂)
```

### Ordinal Suffixes
`get_ordinal_suffix` in `Aux_bands.py` generates labels for file naming — `1st`, `2nd`, `3rd`, `4th` — with correct handling of the `11th`/`12th`/`13th` exceptions.

---

## Remarks

> **Generalisability:** Although the code is provided with the `Table_Rule` and `Chair_Rule` substitutions (stored in `Table_Rule.npy` and `Chair_Rule.npy`), it is not limited to these two examples. For any two-dimensional block substitution of degree `SUB_DEG`, the code can be used directly to compute sampled spectral approximations of the associated Schrödinger operator. Provided the corresponding `.npy` file exists, it suffices to set `RULE_NAME` accordingly and to choose appropriate values for `val`.
