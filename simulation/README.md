# Simulation

Code and per-dataset results for the simulation studies in the manuscript.
Each dataset has n samples, 500 Omics I (X) features and 1,000 Omics II (Y)
features, with a planted module of 20 X and 60 Y features (indices 1-20 and
1-60). Singleton features correlate with the module across layers but not
within it.

## Files

| File | Description |
|------|-------------|
| `common.py` | Data generators, association matrices, TriGer / Greedy 1-step / BSP runners, scoring |
| `bsp.py` | Python implementation of BSP |
| `competitors.R` | DIABLO, sCCA, RGCCA, moCluster, IntNMF and JIVE with the manuscript settings |
| `01_main.py` | Main simulation (500 datasets): TriGer, Greedy 1-step, BSP; writes the data to `data/main/` |
| `01_main_competitors.R` | Main simulation: the six R methods on `data/main/` |
| `02_signal_sample_size.py` | Signal strength (s = 0.5, 1.0, 1.5) x sample size (n = 300-700), 100 datasets each |
| `03_penalty_sensitivity.py` | Each penalty varied over 1.0-1.5, 100 datasets |
| `04_null.py` | 50 null datasets (Y rows permuted): TriGer, Greedy 1-step, BSP; writes the data to `data/null/` |
| `04_null_competitors.R` | Null datasets: the six R methods |
| `results/` | Per-dataset results reported in the manuscript |

## Run

From this folder:

```bash
python 01_main.py
Rscript 01_main_competitors.R          # or in chunks: Rscript 01_main_competitors.R 1 4
python 02_signal_sample_size.py
python 03_penalty_sensitivity.py
python 04_null.py
Rscript 04_null_competitors.R
```

Python needs NumPy, SciPy and pandas. R needs mixOmics, PMA, RGCCA, mogsa,
IntNMF and r.jive. `data/` (about 1.5 GB for the main simulation) is not
tracked by git.

## Notes

- Penalties: code `lambda1`, `lambda2`, `lambda3` correspond to the
  manuscript's lambda_0 (cross-omics), lambda_1 and lambda_2 (within-omics).
- Main simulation seeds: of seeds 1-500, the 89 with a < 1.20 and q > 18 are
  replaced by seeds 501-589 drawn from restricted ranges (`main_seeds()` in
  `common.py`).
- `results/null_r.csv` has no true-positive counts for JIVE and IntNMF.
- Per-dataset runtimes were not stored; `results/runtime_mean_sec.json` has the
  mean runtime per method. The scripts record `runtime_sec` when re-run.
