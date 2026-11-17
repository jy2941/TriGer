"""Shared code for the simulation studies: data generators, association
matrices, TriGer / Greedy 1-step / BSP runners and scoring."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
          "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(v, "1")
import sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from triger import greedy_3step
from bsp import largest_bimodule

# Code lambda1 / lambda2 / lambda3 = manuscript lambda_0 / lambda_1 / lambda_2.
TRIGER = dict(lambda1=1.1, lambda2=1.8, lambda3=1.8,
              X_size=(10, 15), Y_size=(10, 250), min_community_size=2)
GREEDY1 = dict(lambda1=1.1, lambda2=1.1, lambda3=1.1,
               X_size=(10, 10**6), Y_size=(10, 10**6), min_community_size=2)


# ----------------------------------------------------------------------
# Main simulation generator
# ----------------------------------------------------------------------
def sim_v(seed, n_sing=15, true_fs=1.2, gx=1.8, true_noise=1.2,
          fy=3.8, y_noise=1.0, sing_fs=1.4, sing_noise=1.7,
          n_clust=3, clust_fs=2.4,
          n=200, pX=500, pY=1000, nX_true=20, nY_true=60):
    """True module: 20 X features (Fs + Gx) and 60 Y features (Fs).
    n_sing X singletons load on Fs only; the first n_clust also share Fc.
    Parameter names -> manuscript: n_sing = q, true_fs = a, gx = g, fy = b,
    sing_noise = sigma."""
    rng = np.random.default_rng(seed)
    Fs = rng.standard_normal(n); Gx = rng.standard_normal(n); Fc = rng.standard_normal(n)
    X = rng.standard_normal((n, pX)); Y = rng.standard_normal((n, pY))
    iaX = np.arange(nX_true); iaY = np.arange(nY_true)
    singX = np.arange(nX_true, nX_true + n_sing)
    sgn = rng.choice([-1, 1], nX_true)
    for k, i in enumerate(iaX):
        X[:, i] = true_fs * sgn[k] * Fs + gx * sgn[k] * Gx + rng.standard_normal(n) * true_noise
    for j in iaY:
        Y[:, j] = fy * rng.choice([-1, 1]) * Fs + rng.standard_normal(n) * y_noise
    for c, i in enumerate(singX):
        base = sing_fs * rng.choice([-1, 1]) * Fs + rng.standard_normal(n) * sing_noise
        if c < n_clust:
            base = base + clust_fs * Fc
        X[:, i] = base
    return X, Y, iaX, iaY, singX


def main_params(seed):
    """Randomized parameters for seeds 1-500 and the replacement seeds 501-589."""
    if seed <= 500:
        meta = np.random.default_rng(30_000 + seed)
        ns = int(meta.integers(6, 27)); tf = meta.uniform(0.95, 1.45)
    else:
        meta = np.random.default_rng(40_000 + seed)
        ns = int(meta.integers(6, 19)); tf = meta.uniform(1.20, 1.45)
    gx = meta.uniform(1.3, 2.2); fy = meta.uniform(2.8, 4.4); sn = meta.uniform(1.6, 2.1)
    return dict(n_sing=ns, true_fs=tf, gx=gx, fy=fy, sing_noise=sn)


def main_seeds():
    """Seeds 1-500 without the 89 with a < 1.20 and q > 18, plus seeds 501-589."""
    keep = [s for s in range(1, 501)
            if not (main_params(s)["true_fs"] < 1.20 and main_params(s)["n_sing"] > 18)]
    return keep + list(range(501, 590))


def sim_main(seed):
    return sim_v(seed, **main_params(seed))


# ----------------------------------------------------------------------
# Generator for the signal / sample-size grid and penalty sensitivity
# ----------------------------------------------------------------------
def sim_s3(seed, s=1.0, n=200, n_sing=16, gx=1.75, sing_noise=1.85,
           n_sing_y=10, sing_noise_y=1.5, true_noise=1.2, y_noise=1.0,
           n_clust=3, clust_fs=2.4, pX=500, pY=1000, nX_true=20, nY_true=60):
    """sim_v with every loading on Fs scaled by s, plus n_sing_y Y singletons
    (s * (1.4 Gx + 1.0 Fs) + 1.5 e) drawn from a separate random stream."""
    X, Y, iaX, iaY, singX = sim_v(seed, n_sing=n_sing, true_fs=1.2 * s, gx=gx,
                                  true_noise=true_noise, fy=3.6 * s, y_noise=y_noise,
                                  sing_fs=1.4 * s, sing_noise=sing_noise, n_clust=n_clust,
                                  clust_fs=clust_fs, n=n, pX=pX, pY=pY,
                                  nX_true=nX_true, nY_true=nY_true)
    rng = np.random.default_rng(seed)
    Fs = rng.standard_normal(n); Gx = rng.standard_normal(n)
    ry = np.random.default_rng(10_000_000 + seed)
    singY = np.arange(nY_true, nY_true + n_sing_y)
    for j in singY:
        sg, sf = ry.choice([-1, 1]), ry.choice([-1, 1])
        Y[:, j] = s * (1.4 * sg * Gx + 1.0 * sf * Fs) + ry.standard_normal(n) * sing_noise_y
    return X, Y, iaX, iaY, singX, singY


# ----------------------------------------------------------------------
# Association matrices, methods and scoring
# ----------------------------------------------------------------------
def build_PR(X, Y):
    """P = -log10 p of Pearson correlations over [X | Y]; R = signed correlations."""
    n, m = X.shape; nn = Y.shape[1]
    R = np.nan_to_num(np.corrcoef(np.hstack([X, Y]).T))
    r = np.clip(R, -0.999999, 0.999999)
    t = r * np.sqrt((n - 2) / (1 - r ** 2))
    P = -np.log10(np.clip(2 * stats.t.sf(np.abs(t), df=n - 2), 1e-300, 1.0))
    np.fill_diagonal(P, 0.0)
    return P, R, m, nn


def selected(blocks):
    sx, sy = set(), set()
    for b in blocks:
        sx |= set(np.asarray(b[0]).tolist()); sy |= set(np.asarray(b[1]).tolist())
    return sx, sy


def run_methods(X, Y, triger=TRIGER, PR=None, methods=("TriGer", "Greedy 1-step", "BSP")):
    """Selected (X, Y) feature sets for each method."""
    if PR is None and ("TriGer" in methods or "Greedy 1-step" in methods):
        PR = build_PR(X, Y)
    out = {}
    if "TriGer" in methods:
        out["TriGer"] = selected(greedy_3step(*PR, **triger))
    if "Greedy 1-step" in methods:
        out["Greedy 1-step"] = selected(greedy_3step(*PR, **GREEDY1))
    if "BSP" in methods:
        bx, by = largest_bimodule(X, Y, alpha=0.05, frac=0.5)
        out["BSP"] = (set(bx.tolist()), set(by.tolist())) if len(bx) and len(by) else (set(), set())
    return out


def score(sel, iaX, iaY, singX=(), singY=()):
    sx, sy = sel
    tX, tY = set(np.asarray(iaX).tolist()), set(np.asarray(iaY).tolist())
    reX = len(sx & tX) / len(tX); reY = len(sy & tY) / len(tY)
    puX = len(sx & tX) / len(sx) if sx else 0.0; puY = len(sy & tY) / len(sy) if sy else 0.0
    f1 = lambda p, r: 2 * p * r / (p + r) if p + r > 0 else 0.0
    return dict(F1=(f1(puX, reX) + f1(puY, reY)) / 2, precision=(puX + puY) / 2,
                recall=(reX + reY) / 2, nX=len(sx), nY=len(sy),
                singX=len(sx & set(np.asarray(singX).tolist())),
                singY=len(sy & set(np.asarray(singY).tolist())))
