"""Null datasets: seeds 1-50 of sim_v with fixed parameters (q = 15, a = 1.2,
g = 1.8, b = 3.8, sigma = 1.7), Y rows permuted independently per dataset.
Runs TriGer, Greedy 1-step and BSP on the null and the original (planted) data,
and writes the null data to data/null/ for 04_null_competitors.R."""
import os
import numpy as np, pandas as pd
import common as c

OUT = os.path.join(c.HERE, "data", "null"); os.makedirs(OUT, exist_ok=True)
rows = []
for s in range(1, 51):
    X, Y, iaX, iaY, _ = c.sim_v(s)
    Yn = Y[np.random.default_rng(777 + s).permutation(Y.shape[0])]
    np.savetxt(os.path.join(OUT, f"X_{s}.csv"), X, delimiter=",", fmt="%.6g")
    np.savetxt(os.path.join(OUT, f"Y_{s}.csv"), Yn, delimiter=",", fmt="%.6g")
    for setting, YY in (("null", Yn), ("planted", Y)):
        for m, (sx, sy) in c.run_methods(X, YY).items():
            rows.append(dict(setting=setting, method=m, seed=s, nX=len(sx), nY=len(sy),
                             tpX=len(sx & set(iaX.tolist())), tpY=len(sy & set(iaY.tolist())),
                             nonempty=bool(sx) and bool(sy)))
pd.DataFrame(rows).to_csv(os.path.join(c.HERE, "results", "null_python.csv"), index=False)
