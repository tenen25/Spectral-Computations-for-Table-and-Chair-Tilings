# Spectral Computation Code

> Computes sampled spectral bands of periodic approximants associated with substitution-generated patterns.

---

## Table of Contents

- [Requirements](#requirements)
- [Setup](#setup)
- [File Structure](#file-structure)
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

All computation results will be saved there.

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
| `Aux_matrix.py` | `sub_diag`, `nested_phase_points`, `op_mat_no_phase`, `sample_nth_eig` |
| `Aux_bands.py` | `save_bands`, `read_bands`, `print_several_bands`, `get_ordinal_suffix` |

---

## Workflow

When `main.py` is executed, it runs the following steps:

1. **Read substitution rule**
   ```python
   tile_rule = rule_defn_read(RULE_NAME)
   ```

2. **Generate substituted patch** — produces a `2ⁿ × 2ⁿ` patch after `n` iterations
   ```python
   patch = iterate_sub(n, SUB_DEG, start_tile, tile_rule)
   ```

3. **Compute and save sampled bands**
   ```python
   save_bands(n, res, 0, patch, val)
   ```

4. **Plot the saved bands**
   ```python
   print_several_bands(...)
   ```

### Core Call Chain

```
main → save_bands → sample_nth_eig → op_mat_no_phase → eigvalsh
main → print_several_bands → read_bands
```

---

## Key Concepts

### Patch Numbering
The 2D patch is stored as `patch[i][j]` and flattened **row-major**:

```
(i, j)  →  k = i × wid + j
```

For a 4×4 patch:

```
 0   1   2   3
 4   5   6   7
 8   9  10  11
12  13  14  15
```

### Phase Sampling
Bloch phases are sampled over all rational points of the form:

```
θ = π · a/q,   where gcd(a, q) = 1,   q ≤ res
```

Duplicates are removed automatically using Python's `fractions.Fraction` and `set`.

### Sampled Bands
For each eigenvalue branch `k`, the code records the sampled min and max over all phase pairs `(θ₁, θ₂)`:

```
band_min[k] = min λₖ(θ₁, θ₂)
band_max[k] = max λₖ(θ₁, θ₂)
```

### Ordinal Suffixes
`get_ordinal_suffix` in `Aux_bands.py` generates labels for file naming — `1st`, `2nd`, `3rd`, `4th` — with correct handling of the `11th`/`12th`/`13th` exceptions.
