"""
Python implementation of TriGer greedy peeling.

Mirrors the MATLAB functions in this folder (greedy_peeling_XY_one/all,
greedy_peeling_X_one/all, sort_result) and adds greedy_3step, the full
cross-omics + within-omic refinement pipeline. All indices are 0-based.
"""
import numpy as np


# ----------------------------------------------------------------------
# Step 1 : cross-omics (X x Y) peeling
# ----------------------------------------------------------------------
def greedy_peeling_XY_one(Wp, lam):
    """Peel one dense submatrix from a non-negative m x n matrix Wp.

    Iteratively removes whichever remaining row or column has the smallest
    marginal sum (row/col balanced by the current aspect ratio), tracking the
    size-penalised density  sum(remaining) / (|X|*|Y|)^(lam/2)  and returning
    the submatrix that maximised it.

    Returns (select_X, select_Y) : 0-based index arrays into Wp's rows/cols.
    """
    Wp = np.asarray(Wp, dtype=np.float64)
    m, n = Wp.shape
    C = Wp.sum(axis=0).astype(float)      # column sums
    R = Wp.sum(axis=1).astype(float)      # row sums
    remove_X, remove_Y = [], []
    remove_X_final, remove_Y_final = [], []
    len_X, len_Y = m, n
    W_density = 0.0

    for _ in range(m + n):
        Cc = C.copy(); Rr = R.copy()
        if remove_Y: Cc[remove_Y] = np.inf
        if remove_X: Rr[remove_X] = np.inf
        IndT = int(np.argmin(Cc)); dT = Cc[IndT]      # smallest column marginal
        IndS = int(np.argmin(Rr)); dS = Rr[IndS]      # smallest row marginal
        if dT == np.inf or dS == np.inf:
            break
        c = len_X / len_Y if len_Y > 0 else np.inf
        if c * dS <= dT:                              # remove a row (X feature)
            remove_X.append(IndS)
            C = C - Wp[IndS, :]
            if remove_Y: C[remove_Y] = 0.0
            Wp_sum = C.sum()
            len_X -= 1
        else:                                         # remove a col (Y feature)
            remove_Y.append(IndT)
            R = R - Wp[:, IndT]
            if remove_X: R[remove_X] = 0.0
            Wp_sum = R.sum()
            len_Y -= 1
        if len_X <= 0 or len_Y <= 0:
            break
        score = Wp_sum / (len_X * len_Y) ** (lam / 2.0)
        if np.isfinite(score) and score > W_density:
            W_density = score
            remove_X_final = remove_X.copy()
            remove_Y_final = remove_Y.copy()

    select_X = np.array(sorted(set(range(m)) - set(remove_X_final)), dtype=int)
    select_Y = np.array(sorted(set(range(n)) - set(remove_Y_final)), dtype=int)
    return W_density, select_X, select_Y


def greedy_peeling_XY_all(XY, lam=1.5, max_round=5):
    """Repeatedly peel dense X x Y blocks, removing each from the matrix.

    Returns a list of (X_idx, Y_idx, density) tuples (0-based indices into the
    original XY), one per block, only while a block beats the whole-matrix
    density and the round budget is not exhausted.
    """
    XY = np.asarray(XY, dtype=np.float64)
    m, n = XY.shape
    org_density = XY.sum() / (m * n) ** (lam / 2.0)
    remain_X = np.arange(m); remain_Y = np.arange(n)
    Wp = XY.copy()
    result = []
    for i in range(max_round):
        if Wp.size == 0 or Wp.shape[0] < 1 or Wp.shape[1] < 1:
            break
        dens, sX, sY = greedy_peeling_XY_one(Wp, lam)
        if dens > org_density and len(sX) > 0 and len(sY) > 0:
            gX = remain_X[sX]; gY = remain_Y[sY]
            result.append((np.sort(gX), np.sort(gY), dens))
            remain_X = np.array(sorted(set(range(m)) - set(np.concatenate([r[0] for r in result]).tolist())), dtype=int)
            remain_Y = np.array(sorted(set(range(n)) - set(np.concatenate([r[1] for r in result]).tolist())), dtype=int)
            if len(remain_X) == 0 or len(remain_Y) == 0:
                break
            Wp = XY[np.ix_(remain_X, remain_Y)]
        else:
            break
    return result


# ----------------------------------------------------------------------
# Step 2 : within-omic (X x X, symmetric) peeling
# ----------------------------------------------------------------------
def greedy_peeling_X_one(Wp, lam):
    """Peel one dense community from a symmetric within-omic weight matrix
    (signed correlation r; negative entries penalise grouping anti-correlated features).

    Removes the smallest-column-sum node iteratively, tracking
    sum(remaining) / (m-i)^lam .  Returns (density, select_X 0-based).
    """
    Wp = np.asarray(Wp, dtype=np.float64)
    m = Wp.shape[0]
    C = Wp.sum(axis=0).astype(float)
    remove_X, remove_X_final = [], []
    W_density = 0.0
    for i in range(1, m + 1):
        Cc = C.copy()
        if remove_X: Cc[remove_X] = np.inf
        Ind = int(np.argmin(Cc))
        remove_X.append(Ind)
        C = C - Wp[Ind, :]
        C[remove_X] = 0.0
        Wp_sum = C.sum()
        denom = (m - i) ** lam
        score = Wp_sum / denom if denom > 0 else np.inf
        if np.isfinite(score) and score > W_density:
            W_density = score
            remove_X_final = remove_X.copy()
    select_X = np.array(sorted(set(range(m)) - set(remove_X_final)), dtype=int)
    return W_density, select_X


def greedy_peeling_X_all(Wp, lam, min_size=1):
    """Repeatedly peel symmetric communities. Returns list of (idx, density).

    Parameters
    ----------
    Wp : symmetric within-omic weight matrix (signed correlation r submatrix).
    lam : size penalty.
    min_size : minimum community size to keep. Communities smaller than this are
        discarded from the result. The default of 1 keeps every community
        (original behaviour). Set min_size=2 to drop degenerate single-feature
        communities: an isolated feature (e.g. a within-omic singleton) is peeled
        into its own size-1 community whose self-density (1.0) clears the
        acceptance threshold, so with min_size=1 it is retained. Dropping size-1
        communities removes such singletons during within-omic refinement.
    """
    Wp = np.asarray(Wp, dtype=np.float64)
    assert np.allclose(Wp, Wp.T, atol=1e-8), "Matrix must be symmetric."
    m = Wp.shape[0]
    org_density = Wp.sum() / (m ** lam)
    X = Wp.copy()
    remain_X = np.arange(m)
    result = []
    cur = Wp.copy()
    while True:
        dens, sel = greedy_peeling_X_one(cur, lam)
        if dens > org_density and len(sel) > 0:
            result.append((remain_X[sel], dens))
            allsel = np.concatenate([r[0] for r in result])
            remain_X = np.array(sorted(set(range(m)) - set(allsel.tolist())), dtype=int)
            if len(remain_X) == 0:
                break
            cur = X[np.ix_(remain_X, remain_X)]
        else:
            break
    if min_size > 1:
        result = [r for r in result if len(r[0]) >= min_size]
    return result


# ----------------------------------------------------------------------
# sort_result : order features within each block by within-block marginal
# ----------------------------------------------------------------------
def sort_result(Wp, blocks, symmetric=None):
    """Reorder feature indices inside each block by descending within-block
    marginal sum.  `blocks` = list of (X_idx, Y_idx[, ...]); for a symmetric
    (X x X) matrix pass single-index blocks and symmetric=True."""
    Wp = np.asarray(Wp, dtype=np.float64)
    if symmetric is None:
        symmetric = (Wp.shape[0] == Wp.shape[1])
    out = []
    for blk in blocks:
        if symmetric:
            idx = np.asarray(blk[0], dtype=int)
            sub = Wp[np.ix_(idx, idx)]
            order = np.argsort(-sub.sum(axis=1))
            out.append((idx[order],) + tuple(blk[1:]))
        else:
            xi = np.asarray(blk[0], dtype=int); yi = np.asarray(blk[1], dtype=int)
            sub = Wp[np.ix_(xi, yi)]
            ox = np.argsort(-np.abs(sub.sum(axis=1)))
            oy = np.argsort(-np.abs(sub.sum(axis=0)))
            out.append((xi[ox], yi[oy]) + tuple(blk[2:]))
    return out


# ----------------------------------------------------------------------
# Step 3 : the full driver
# ----------------------------------------------------------------------
def greedy_3step(P, R, m, n, lambda1=1.5, lambda2=1.5, lambda3=1.5,
                 X_size=(10, 200), Y_size=(10, 200), min_community_size=2,
                 keep_unrefined_if_empty=True):
    """Full 3-step greedy peeling.

    Parameters
    ----------
    P : (m+n) x (m+n) combined non-negative association matrix (-log10 p or
        correlation); the X x Y sub-block  P[0:m, m:m+n]  drives step 1.
    R : (m+n) x (m+n) combined SIGNED correlation matrix used for within-omic refinement.
    m, n : number of X (type-1) and Y (type-2) features.
    lambda1/2/3 : size penalties in [1, 2] for the XY / X / Y peels.
    X_size, Y_size : (min, max) community-size limits. Blocks smaller than the
        min are dropped; blocks larger than the max are refined by a within-omic
        peel and reordered.
    min_community_size : minimum size for a within-omic community to be retained
        during step-2/step-3 refinement. Default 2 drops degenerate size-1
        communities, i.e. features that are isolated within their own omic
        (within-omic singletons): the within-omic peel separates each such
        feature into its own size-1 community, and keeping it (min=1) lets it
        survive refinement even though it does not cohere with the recovered
        block. Pass min_community_size=1 to restore the previous behaviour of
        keeping every peeled community. Only affects blocks large enough to be
        refined; the 1-step configuration (loose X_size/Y_size caps) never
        triggers refinement and is unaffected by this parameter.
    keep_unrefined_if_empty : fallback for the within-omic refinement (steps 2/3).
        If the within-omic peel accepts no community (no subset of the block is denser
        than the block as a whole, e.g. an uncontaminated homogeneous block, or a lambda
        at the edge of [1, 2]), keep that side of the block unrefined instead of
        discarding it. Default True. Pass False to empty that side instead.

    Returns
    -------
    list of (X_idx, Y_idx) blocks, 0-based indices into the original X / Y feature
    axes (Y indices are returned in 0..n-1, i.e. already offset back from P).
    """
    P = np.asarray(P, dtype=np.float64); R = np.asarray(R, dtype=np.float64)
    Wp = P[0:m, m:m + n]

    # --- step 1: peel X x Y blocks ---
    blocks = greedy_peeling_XY_all(Wp, lambda1)

    # --- filter by minimum size ---
    blocks = [b for b in blocks if len(b[0]) > X_size[0] and len(b[1]) > Y_size[0]]

    refined = []
    for (xi, yi, dens) in blocks:
        xi = np.asarray(xi, dtype=int); yi = np.asarray(yi, dtype=int)
        # --- step 2: refine X side if oversized, using within-omic correlation ---
        if len(xi) > X_size[1]:
            Xsub = R[np.ix_(xi, xi)]; Xsub = (Xsub + Xsub.T) / 2.0
            comm = greedy_peeling_X_all(Xsub, lambda2, min_size=min_community_size)
            comm = sort_result(Xsub, comm, symmetric=True)
            if comm:
                xi = xi[np.concatenate([c[0] for c in comm])]
            elif not keep_unrefined_if_empty:
                xi = xi[np.array([], dtype=int)]
        # --- step 3: refine Y side if oversized ---
        if len(yi) > Y_size[1]:
            idxY = m + yi
            Ysub = R[np.ix_(idxY, idxY)]; Ysub = (Ysub + Ysub.T) / 2.0
            comm = greedy_peeling_X_all(Ysub, lambda3, min_size=min_community_size)
            comm = sort_result(Ysub, comm, symmetric=True)
            if comm:
                yi = yi[np.concatenate([c[0] for c in comm])]
            elif not keep_unrefined_if_empty:
                yi = yi[np.array([], dtype=int)]
        refined.append((xi, yi))
    return refined
