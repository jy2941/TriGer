# TriGer — Dense Community Detection for Multi-Omics Integration

Python, MATLAB and R implementations of the **TriGer** algorithm accompanying our academic paper.

## Overview

TriGer identifies dense bipartite communities from multi-modal omics data by jointly analysing three correlation matrices:
- **X** (m×m): symmetric correlation matrix for data type 1 (e.g., transcriptomics)
- **Y** (n×n): symmetric correlation matrix for data type 2 (e.g., metabolomics)
- **XY** (m×n): cross-association matrix between the two data types

Community size is controlled by a penalty parameter λ ∈ [1, 2]: higher values favour smaller, denser communities.

## Repository Structure

```
TriGer/
├── python/triger/        # Python implementation (methods and plotting)
├── matlab/
│   ├── methods/          # MATLAB implementation of the methods
│   └── visualization/    # MATLAB plotting functions
├── R/
│   ├── methods/          # R implementation of the methods
│   └── visualization/    # R plotting function
├── simulation/           # Simulation studies in the manuscript (Python + R)
└── Real_data/            # Real-data analyses (MATLAB)
    ├── CRC/              # Colorectal cancer
    └── IBD/              # Inflammatory bowel disease
```

The three implementations use the same function names and arguments.

## Functions

| Function | Description | Python | MATLAB | R |
|----------|-------------|:------:|:------:|:-:|
| `greedy_peeling_XY_one` | Extract one dense block from a cross-omics (X-Y) matrix | ✓ | ✓ | ✓ |
| `greedy_peeling_XY_all` | Extract dense X-Y blocks repeatedly | ✓ | ✓ | ✓ |
| `greedy_peeling_X_one` | Extract one dense community from a within-omics (symmetric) matrix | ✓ | ✓ | ✓ |
| `greedy_peeling_X_all` | Extract within-omics communities repeatedly | ✓ | ✓ | ✓ |
| `greedy_3step` | Full pipeline: X-Y block detection, then within-omics refinement of large blocks | ✓ | | ✓ |
| `sort_result` | Order features within each block by their within-block sums | ✓ | ✓ | ✓ |
| `sort_result_XandY` | Order features using the X-X and Y-Y matrices | | ✓ | |
| `plot3in1_v1` | X-Y heatmap with X-X and Y-Y triangles (top / right) | ✓ | ✓ | ✓ |
| `plot3in1_v2` | Same plot with the triangles bottom / left | | ✓ | |
| `plot_all_result` | Heatmap of the detected blocks within the full or selected matrix | ✓ | ✓ | |

Indices are 0-based in Python and 1-based in MATLAB and R. `greedy_peeling_XY_all`
returns at most `max_round` blocks in Python and `max_round - 1` in MATLAB and R.

## Usage

### Python

```python
import sys; sys.path.append("python")
from triger import greedy_peeling_XY_all, sort_result, greedy_3step, plot3in1_v1

# Dense blocks in a non-negative X-Y association matrix (e.g. -log10 p)
blocks = sort_result(XY, greedy_peeling_XY_all(XY, lam=1.5))   # [(X_idx, Y_idx, density), ...]

# Full pipeline: P and R are the (m+n) x (m+n) association and signed correlation
# matrices over [X | Y]; blocks with more than X_size[1] / Y_size[1] features are
# refined by within-omics peeling on R
blocks = greedy_3step(P, R, m, n, lambda1=1.1, lambda2=1.8, lambda3=1.8,
                      X_size=(10, 15), Y_size=(10, 250))
```

### MATLAB

```matlab
addpath('matlab/methods/'); addpath('matlab/visualization/')
result = sort_result(XY, greedy_peeling_XY_all(XY, 1.5));
plot3in1_v1(X, Y, XY, result, true);
```

### R

```r
source("R/methods/greedy.R"); source("R/visualization/plot3in1_v1.R")
blocks <- sort_result(XY, greedy_peeling_XY_all(XY, 1.5))
```

## Parameters

- **λ (lambda)**: size penalty in [1, 2]. Higher values give smaller, denser blocks.
  In `greedy_3step`, `lambda1` is the cross-omics penalty and `lambda2` / `lambda3`
  are the within-omics penalties for X and Y (λ₀, λ₁ and λ₂ in the manuscript).
  Typical values: 1.1–1.6 for X-Y detection, 1.8–1.9 for within-omics refinement.

## Simulation

Simulation studies from the manuscript: the main nine-method comparison, signal
strength and sample size, penalty sensitivity, and null datasets. See
[`simulation/README.md`](simulation/README.md).

## Real Data

- **CRC**: `CRC_correlation_matrix.mat` (precomputed correlation and -log10 p matrices) and `CRC_Analysis.m`
- **IBD**: `IBD_Analysis.m`, comparing CD, UC and non-IBD controls

## Requirements

- Python 3 with NumPy (matplotlib for plotting; SciPy and pandas for the simulations)
- MATLAB R2023a or later, no extra toolboxes
- R 4.0 or later (base graphics); the simulation comparison methods need mixOmics, PMA, RGCCA, mogsa, IntNMF and r.jive
