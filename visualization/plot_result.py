"""
Python versions of plot_all_result.m and plot3in1_v1.m (matplotlib).

`result` is a list of (X_idx, Y_idx, ...) tuples, as returned by
methods/greedy_peeling.py; for a symmetric matrix, (idx, ...) tuples.
Indices are 0-based unless one_based=True.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.patches import Rectangle

MATLAB_GREY = (0.15, 0.15, 0.15)


def matlab_jet(m=256):
    """Exact port of MATLAB's jet(m) colormap (R2019b+ default length 256)."""
    n = int(np.ceil(m / 4))
    u = np.concatenate([np.arange(1, n + 1) / n, np.ones(n - 1), np.arange(n, 0, -1) / n])
    g = int(np.ceil(n / 2)) - (m % 4 == 1) + np.arange(1, len(u) + 1)
    r = g + n
    b = g - n
    r = r[r <= m]
    g = g[g <= m]
    b = b[b >= 1]
    J = np.zeros((m, 3))
    J[r - 1, 0] = u[:len(r)]
    J[g - 1, 1] = u[:len(g)]
    J[b - 1, 2] = u[len(u) - len(b):]
    return ListedColormap(J, name="matlab_jet")


def _unique_keep_last(v):
    """MATLAB: [~,idx]=unique(flip(v),'stable'); v=flip(v(end-idx+1))
    -> keep the LAST occurrence of each value, preserving order."""
    rev = v[::-1]
    _, first = np.unique(rev, return_index=True)
    return rev[np.sort(first)][::-1]


def _as_idx(x):
    return np.atleast_1d(np.asarray(x, dtype=int)).ravel() if x is not None else np.array([], int)


def plot3in1_v1(X, Y, XY, result, line=False, step=None, bar=False,
                one_based=False, figsize=(5.6, 4.2), cbar_pos=(0.15, 0.13, 0.013, 0.45),
                clim_tri=None, clim_mid=None, box_lw=2):
    """
    X, Y   : square matrices (X is m x m, Y is n x n)
    XY     : m x n matrix
    result : p x 2 list, each row = (X_indices, Y_indices) -- same layout as
             the MATLAB cell array (the function swaps columns internally,
             exactly as the original does).
    line   : draw red rectangles around blocks
    step   : (xs, ys) tick step; default (100, 10) as in MATLAB
    bar    : add colorbar for the middle panel (MATLAB: nargin == 7)
    one_based : set True if indices come from MATLAB (1-based)
    cbar_pos  : colorbar [left bottom width height] in figure units, or "right" = just right
                of the plot (past the right triangle tip), bottom/height aligned with the middle square
                (MATLAB hard-coded value by default)
    clim_tri  : colour limits of the two triangles (MATLAB: [-1 1])
    box_lw    : line width of the red block boxes (MATLAB: 2)
    clim_mid  : colour limits of the middle panel (MATLAB: [0 10] if max(A)>1 else [-1 1])
    """
    # result = result(:, [2 1])
    result = [(_as_idx(r[1]), _as_idx(r[0])) for r in result]
    off = 1 if one_based else 0
    select_X = np.concatenate([r[0] for r in result]) - off
    select_Y = np.concatenate([r[1] for r in result]) - off
    select_X = _unique_keep_last(select_X)
    select_Y = _unique_keep_last(select_Y)

    X, Y, XY = np.asarray(X), np.asarray(Y), np.asarray(XY)
    C = Y[np.ix_(select_Y, select_Y)]
    B = X[np.ix_(select_X, select_X)]
    A = XY[np.ix_(select_X, select_Y)]

    m, n = A.shape
    ratio = m / n / 1.2

    cmap = matlab_jet(256)
    with mpl.rc_context({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 10,
        "xtick.color": MATLAB_GREY, "ytick.color": MATLAB_GREY,
        "xtick.labelcolor": MATLAB_GREY, "ytick.labelcolor": MATLAB_GREY,
        "xtick.direction": "in", "ytick.direction": "in",
        "axes.linewidth": 0.5, "axes.unicode_minus": False, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    }):
        fig = plt.figure(figsize=figsize, facecolor="white")
        ax = fig.add_axes([0.13, 0.11, 0.775, 0.815])   # MATLAB default axes Position

        # ---------- Right: diamonds for Y (lower triangle) ----------
        vx = np.array([0, 1, 0, -1]); vy = np.array([1, 0, -1, 0])
        polys, cols = [], []
        for i in range(1, n + 1):
            for j in range(i, n + 1):
                x0 = i - 1 + m * 2 / ratio
                y0 = 2 * n - j * 2 + i
                polys.append(np.column_stack([x0 + vx, y0 + vy]))
                cols.append(C[j - 1, j - i])
        norm1 = Normalize(*(clim_tri or (-1, 1)), clip=True)   # caxis(ax1, [-1 1])
        pc_right = PolyCollection(polys, array=np.array(cols), cmap=cmap, norm=norm1,
                                  edgecolors="face", linewidths=0.2, zorder=1)
        ax.add_collection(pc_right)

        # ---------- Top: diamonds for X (lower triangle) ----------
        vx = np.array([0, 1, 0, -1]) / ratio; vy = np.array([1, 0, -1, 0]) / ratio
        polys, cols = [], []
        for i in range(1, m + 1):
            for j in range(i, m + 1):
                x0 = (j * 2 - i) / ratio
                y0 = (i - 1) / ratio + n * 2
                polys.append(np.column_stack([x0 + vx, y0 + vy]))
                cols.append(B[j - 1, j - i])
        pc_top = PolyCollection(polys, array=np.array(cols), cmap=cmap, norm=norm1,
                                edgecolors="face", linewidths=0.2, zorder=1)
        ax.add_collection(pc_top)

        # ---------- Middle square: XY ----------
        vxs = np.array([0, 2 / ratio, 2 / ratio, 0]); vys = np.array([0, 0, 2, 2])
        polys, cols = [], []
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                x0 = (i - 1) * 2 / ratio
                y0 = (j - 1) * 2
                polys.append(np.column_stack([x0 + vxs, y0 + vys]))
                cols.append(A[i - 1, n - j])
        clim2 = clim_mid or ((0, 10) if np.nanmax(A) > 1 else (-1, 1))   # same rule as MATLAB
        pc_mid = PolyCollection(polys, array=np.array(cols), cmap=cmap,
                                norm=Normalize(*clim2, clip=True),
                                edgecolors="face", linewidths=0.2, zorder=2)
        ax.add_collection(pc_mid)

        # ---------- axis equal tight + linkaxes (union of both data ranges) ----------
        ax.set_xlim(0, max(2 * m / ratio, n + 2 * m / ratio))
        ax.set_ylim(0, 2 * n + m / ratio)
        ax.set_aspect("equal", adjustable="box", anchor="C")
        for s in ax.spines.values():                    # Axle 'none' + box off
            s.set_visible(False)
        ax.set_facecolor("none")

        # ---------- Clean tick labels ----------
        xs, ys = (100, 10) if step is None else (step[0], step[1])
        clean = lambda maxv, base: np.arange(base, np.floor(maxv / base) * base + 1, base)
        xTickIdx = clean(m, xs).astype(int)
        yTickIdx = clean(n, ys).astype(int)[::-1]
        ax.set_xticks((xTickIdx - 1) * 2 / ratio, [str(v) for v in xTickIdx])
        ax.set_yticks((n - yTickIdx) * 2 + 1, [str(v) for v in yTickIdx])
        ax.tick_params(top=False, right=False)

        # MATLAB TickLength 0.01 of the longest axes side
        ax.apply_aspect()
        bb = ax.get_position()
        longest_pt = max(bb.width * figsize[0], bb.height * figsize[1]) * 72
        ax.tick_params(length=0.01 * longest_pt)

        if bar:
            if isinstance(cbar_pos, str) and cbar_pos == "right":
                sq_h = bb.height * (2 * n) / (2 * n + m / ratio)      # middle-square height
                cbar_pos = (bb.x1 + 0.02, bb.y0, 0.013, sq_h)
            cax = fig.add_axes(list(cbar_pos))
            cb = fig.colorbar(pc_mid, cax=cax)
            cb.ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{v:g}"))
            cb.outline.set_linewidth(0.5)
            cb.outline.set_edgecolor(MATLAB_GREY)
            cax.tick_params(direction="in", length=0.01 * longest_pt,
                            color=MATLAB_GREY, labelcolor=MATLAB_GREY)

        # ---------- Lines separating blocks ----------
        if line:
            L = np.array([[len(r[0]), len(r[1])] for r in result], dtype=float)
            L[:, 0] *= 2 / ratio
            L[:, 1] *= 2
            x, y = 0.0, L[:, 1].sum()
            i = 0
            while i < len(L) and L[i, 0] != 0:
                ax.add_patch(Rectangle((x, y - L[i, 1]), L[i, 0], L[i, 1], fill=False,
                                       edgecolor=(1, 0, 0), linewidth=box_lw, zorder=3))
                y -= L[i, 1]
                x += L[i, 0]
                i += 1

    return fig, ax


def plot_all_result(Wp, result, full=False, line=True, crange=None, figsize=(8, 6)):
    """Heatmap of Wp reordered by the detected communities.
    full=True shows the whole matrix (selected features first) with one box
    around the selected region; full=False shows only the selected region
    with one box per community."""
    Wp = np.asarray(Wp)
    m, n = Wp.shape
    if m == n:
        result = [(r[0], r[0]) for r in result]
    rows = [_as_idx(r[0]) for r in result]
    cols = [_as_idx(r[1]) for r in result]
    sel_X = np.concatenate(rows) if rows else np.array([], int)
    sel_Y = np.concatenate(cols) if cols else np.array([], int)
    if full:
        M = Wp[np.ix_(np.r_[sel_X, np.setdiff1d(np.arange(m), sel_X)],
                      np.r_[sel_Y, np.setdiff1d(np.arange(n), sel_Y)])]
    else:
        M = Wp[np.ix_(sel_X, sel_Y)]
    if crange is None:
        crange = (np.nanmin(M), np.nanmax(M))
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(M, aspect="auto", cmap=matlab_jet(256), vmin=crange[0], vmax=crange[1],
                   interpolation="nearest")
    fig.colorbar(im, ax=ax)
    if line:
        if full:
            ax.add_patch(Rectangle((-0.5, -0.5), len(sel_Y), len(sel_X), fill=False, ec="r", lw=2))
        else:
            x = y = -0.5
            for r, c in zip(rows, cols):
                if len(r) == 0:
                    break
                ax.add_patch(Rectangle((x, y), len(c), len(r), fill=False, ec="r", lw=2))
                x += len(c); y += len(r)
    return fig, ax
