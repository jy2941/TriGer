"""Main simulation: TriGer, Greedy 1-step and BSP on the 500 datasets.
Also writes each dataset to data/main/seed<k>/ for 01_main_competitors.R."""
import os, time
import numpy as np, pandas as pd
import common as c

OUT = os.path.join(c.HERE, "data", "main")
rows = []
for s in c.main_seeds():
    X, Y, iaX, iaY, singX = c.sim_main(s)
    d = os.path.join(OUT, f"seed{s}"); os.makedirs(d, exist_ok=True)
    np.savetxt(os.path.join(d, "X.csv"), X, delimiter=",", fmt="%.6g")
    np.savetxt(os.path.join(d, "Y.csv"), Y, delimiter=",", fmt="%.6g")
    for m in ("TriGer", "Greedy 1-step", "BSP"):
        t0 = time.time()
        sel = c.run_methods(X, Y, methods=(m,))[m]
        sc = c.score(sel, iaX, iaY, singX)
        rows.append(dict(method=m, seed=s, F1=sc["F1"], precision=sc["precision"],
                         recall=sc["recall"], runtime_sec=time.time() - t0))
pd.DataFrame(rows).to_csv(os.path.join(c.HERE, "results", "main_python.csv"), index=False)
