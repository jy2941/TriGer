"""Penalty sensitivity: each of lambda1, lambda2, lambda3 (manuscript lambda_0,
lambda_1, lambda_2) varied over 1.0-1.5 with the others at their defaults;
s = 1.0, n = 200, seeds 1-100."""
import os
import pandas as pd
from multiprocessing import Pool
import common as c

TRI = dict(c.TRIGER, Y_size=(10, 60))
GRID = (1.0, 1.1, 1.2, 1.3, 1.4, 1.5)


def task(seed):
    X, Y, iaX, iaY, singX, singY = c.sim_s3(seed, s=1.0, n=200)
    PR = c.build_PR(X, Y)
    return [dict(param=p, value=v, seed=seed,
                 **c.score(c.run_methods(X, Y, triger=dict(TRI, **{p: v}), PR=PR, methods=("TriGer",))["TriGer"],
                           iaX, iaY, singX, singY))
            for p in ("lambda1", "lambda2", "lambda3") for v in GRID]


if __name__ == "__main__":
    with Pool(os.cpu_count()) as p:
        rows = [r for res in p.map(task, range(1, 101), chunksize=2) for r in res]
    pd.DataFrame(rows).to_csv(os.path.join(c.HERE, "results", "penalty_sensitivity.csv"), index=False)
