# Spectral Computation Code

This repository provides a Python implementation for computing **sampled spectral bands of periodic approximants** associated with substitution-generated patterns.

---

## Overview

The code constructs substitution patches, builds finite Schrödinger operators, and computes spectral bands via Bloch phase sampling.

### Conceptual workflow

substitution rule→ iterated patch→ finite Hamiltonian→ Bloch phase sampling→ band edges→ saving and plotting


### Main execution chain

rule_defn_read→ iterate_sub→ save_bands→ sample_nth_eig→ op_mat_no_phase→ eigvalsh

### Plotting workflow

read_bands → print_several_bands

---

## Requirements

### Python libraries

- `numpy`
- `matplotlib`
- `fractions`

### Folder setup

Create a subfolder named: Finished computations

This is required for saving output data.

---

## Practical Description

The code performs the following steps:

1. Generate a substitution patch of size `2^n × 2^n`
2. Flatten the patch into a 1D index set (row-major order)
3. Construct a finite Schrödinger matrix:
   - on-site potentials
   - nearest-neighbour hopping
4. Add Bloch boundary conditions (horizontal & vertical)
5. Sample eigenvalues over rational phase pairs:

(π a/q, π a'/q'), with q, q' ≤ res

6. Record min/max values for each eigenvalue branch
7. Save and plot the resulting spectral bands

---

## File Structure

### `main.py`
Driver script:
- reads substitution rules
- generates patches
- computes spectral bands
- produces plots

### `Sub_Rule.py`
Implements substitution dynamics:
- `rule_defn_read`
- `sub_on_patch`
- `iterate_sub`

### Patch Structure and Indexing

- Patches are stored as 2D arrays: `patch[i][j]`  
- Indexing:
  - `i`: row (top to bottom)  
  - `j`: column (left to right)  

### Flattening rule

(i, j) -> k = i * wid + j

This gives row-major ordering:

0, 1, ..., wid-1
wid, ..., 2*wid-1
...
