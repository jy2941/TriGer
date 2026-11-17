"""
Python port of the Bimodule Search Procedure (BSP)
(Dewaskar et al., arXiv:2009.05079; reference R package: cbce).

Method: from a seed feature set on one omic, iteratively extract the
associated features on the other omic and alternate back until the
(X, Y) feature sets stabilise -> one "bimodule". Association uses a
Fisher z-transform of the Pearson correlation with BH-FDR control on the
*pairwise* p-values (NOT an aggregate chi-square, which is anti-conservative
when the target features are themselves a correlated block). A target
feature joins the module if at least `frac` of its pairwise correlations to
the current set are individually significant.
"""
import numpy as np
from scipy.stats import norm


def _bh(pvals, alpha):
    p = np.asarray(pvals); m = len(p)
    ranked = np.sort(p)
    thr = alpha * np.arange(1, m + 1) / m
    passed = ranked <= thr
    if not passed.any():
        return np.zeros(m, dtype=bool)
    cut = ranked[np.max(np.where(passed)[0])]
    return p <= cut


def _extract(A_idx, Zsrc, Ztgt, n, alpha=0.05, frac=0.5):
    """Target features whose pairwise correlation to >=frac of set A is
    individually BH-significant. Robust to within-set dependence."""
    if len(A_idx) == 0:
        return np.array([], dtype=int)
    Rc = (Zsrc[:, A_idx].T @ Ztgt) / n
    Rc = np.clip(Rc, -0.9995, 0.9995)
    z = np.arctanh(Rc) * np.sqrt(n - 3)
    pv = 2 * norm.sf(np.abs(z))
    sig = _bh(pv.ravel(), alpha).reshape(pv.shape)
    return np.where(sig.mean(axis=0) >= frac)[0]


def bsp(X, Y, alpha=0.05, frac=0.5, min_size=3, max_seeds=120,
        maxit=25, rng=None):
    """Run BSP. Returns list of (Xset, Yset) index arrays (0-based), one per
    stable bimodule, de-duplicated by X-set."""
    rng = np.random.default_rng(0) if rng is None else rng
    n = X.shape[0]
    Zx = (X - X.mean(0)) / X.std(0)
    Zy = (Y - Y.mean(0)) / Y.std(0)
    seeds = rng.permutation(X.shape[1])[:max_seeds]
    found, visited = [], set()
    for s in seeds:
        if s in visited:
            continue
        Xset = np.array([s]); Yset = np.array([], dtype=int)
        for _ in range(maxit):
            Yset = _extract(Xset, Zx, Zy, n, alpha, frac)
            if len(Yset) < min_size:
                Xset = np.array([], dtype=int); break
            newX = _extract(Yset, Zy, Zx, n, alpha, frac)
            if len(newX) < min_size:
                Xset = np.array([], dtype=int); break
            if set(newX.tolist()) == set(Xset.tolist()):
                Xset = newX; break
            Xset = newX
        if len(Xset) >= min_size and len(Yset) >= min_size:
            visited.update(Xset.tolist())
            found.append((np.sort(Xset), np.sort(Yset)))
    uniq = {tuple(xs.tolist()): (xs, ys) for xs, ys in found}
    return list(uniq.values())


def largest_bimodule(X, Y, **kw):
    """Convenience: return the single largest (Xset, Yset) bimodule, or
    (empty, empty) if none found."""
    bims = bsp(X, Y, **kw)
    if not bims:
        return np.array([], dtype=int), np.array([], dtype=int)
    return sorted(bims, key=lambda b: -(len(b[0]) + len(b[1])))[0]
