## Null datasets: the six R comparison methods on the data written by 04_null.py.
source("competitors.R")
rows <- list()
for (sd in 1:50) {
  X <- as.matrix(read.csv(sprintf("data/null/X_%d.csv", sd), header = FALSE))
  Y <- as.matrix(read.csv(sprintf("data/null/Y_%d.csv", sd), header = FALSE))
  for (meth in METHODS) {
    set.seed(sd)
    s <- run_competitor(meth, X, Y); sX <- unique(s$X); sY <- unique(s$Y)
    rows[[length(rows) + 1]] <- data.frame(setting = "null", method = meth, seed = sd, nX = length(sX), nY = length(sY),
                                           tpX = sum(sX <= 20), tpY = sum(sY <= 60),
                                           nonempty = length(sX) > 0 && length(sY) > 0)
  }
}
write.csv(do.call(rbind, rows), "results/null_r.csv", row.names = FALSE)
