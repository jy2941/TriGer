# R version of plot_all_result.m (base graphics)
#
# Heatmap of Wp reordered by the detected blocks. full = TRUE shows the whole
# matrix (selected features first) with one box around the selected region;
# full = FALSE shows only the selected region with one box per block.
# result: list of blocks, each list(X, Y, ...) (1-based); for a symmetric
# matrix only X is used.

if (!exists("matlab_jet")) {
  matlab_jet <- function(m = 256) {
    n <- ceiling(m / 4)
    u <- c((1:n) / n, rep(1, n - 1), (n:1) / n)
    g <- ceiling(n / 2) - (m %% 4 == 1) + seq_along(u)
    r <- g + n; b <- g - n
    r <- r[r <= m]; g <- g[g <= m]; b <- b[b >= 1]
    J <- matrix(0, m, 3)
    J[r, 1] <- u[seq_along(r)]; J[g, 2] <- u[seq_along(g)]
    J[b, 3] <- u[(length(u) - length(b) + 1):length(u)]
    rgb(J[, 1], J[, 2], J[, 3])
  }
}

plot_all_result <- function(Wp, result, full = FALSE, line = TRUE, crange = NULL) {
  Wp <- as.matrix(Wp); m <- nrow(Wp); n <- ncol(Wp)
  rows <- lapply(result, function(r) as.integer(r[[1]]))
  cols <- if (m == n) rows else lapply(result, function(r) as.integer(r[[2]]))
  keep <- lengths(rows) > 0; rows <- rows[keep]; cols <- cols[keep]
  sX <- unlist(rows); sY <- unlist(cols)
  M <- if (full) Wp[c(sX, setdiff(seq_len(m), sX)), c(sY, setdiff(seq_len(n), sY)), drop = FALSE]
       else Wp[sX, sY, drop = FALSE]
  if (is.null(crange)) crange <- range(M, na.rm = TRUE)
  M <- pmin(pmax(M, crange[1]), crange[2])
  H <- nrow(M); cmap <- matlab_jet(256)

  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(1:2, 1), widths = c(12, 1))
  par(mar = c(3, 3.5, 1, 0.5))
  image(seq_len(ncol(M)), seq_len(H), t(M[H:1, , drop = FALSE]), col = cmap, zlim = crange,
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  axis(1)
  yt <- pretty(c(1, H)); yt <- yt[yt >= 1 & yt <= H]
  axis(2, at = H + 1 - yt, labels = yt, las = 1)
  if (line) {
    if (full) {
      rect(0.5, H - length(sX) + 0.5, length(sY) + 0.5, H + 0.5, border = "red", lwd = 2)
    } else {
      r0 <- 0; c0 <- 0
      for (k in seq_along(rows)) {
        a <- length(rows[[k]]); b <- length(cols[[k]])
        rect(c0 + 0.5, H - r0 - a + 0.5, c0 + b + 0.5, H - r0 + 0.5, border = "red", lwd = 2)
        r0 <- r0 + a; c0 <- c0 + b
      }
    }
  }
  par(mar = c(3, 0.5, 1, 3))
  z <- seq(crange[1], crange[2], length.out = 256)
  image(1, z, matrix(z, 1), col = cmap, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1); box()
  invisible(M)
}
