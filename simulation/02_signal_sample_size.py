"""Signal strength x sample size: s in {0.5, 1.0, 1.5}, n in {300, ..., 700},
100 datasets each, with 10 Y singletons and the Y refinement threshold at 60."""
import os
import pandas as pd
from multiprocessing import Pool
import common as c

SIGNAL = {"weak": 0.5, "moderate": 1.0, "strong": 1.5}
TRI = dict(c.TRIGER, Y_size=(10, 60))


def task(args):
    level, n, seed = args
    X, Y, iaX, iaY, singX, singY = c.sim_s3(seed, s=SIGNAL[level], n=n)
    return [dict(level=level, s=SIGNAL[level], n=n, seed=seed, method=m, **c.score(sel, iaX, iaY, singX, singY))
            for m, sel in c.run_methods(X, Y, triger=TRI).items()]


if __name__ == "__main__":
    tasks = [(lv, n, sd) for lv in SIGNAL for n in (300, 400, 500, 600, 700) for sd in range(1, 101)]
    with Pool(os.cpu_count()) as p:
        rows = [r for res in p.imap_unordered(task, tasks, chunksize=4) for r in res]
    pd.DataFrame(rows).sort_values(["level", "n", "method", "seed"]).to_csv(
        os.path.join(c.HERE, "results", "signal_sample_size.csv"), index=False)
