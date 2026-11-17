## Main simulation: the six R comparison methods on the datasets written by 01_main.py.
## Usage: Rscript 01_main_competitors.R [chunk n_chunks]   (JIVE takes ~40 s per dataset)
source("competitors.R")
args <- as.integer(commandArgs(trailingOnly = TRUE))
seeds <- sort(as.integer(sub("seed", "", list.files("data/main"))))
if (length(args) == 2) seeds <- seeds[(seq_along(seeds) - 1) %% args[2] == args[1] - 1]
out <- if (length(args) == 2) sprintf("results/main_r_chunk%d.csv", args[1]) else "results/main_r.csv"

f1 <- function(p, r) ifelse(p + r > 0, 2 * p * r / (p + r), 0)
rows <- list()
for (sd in seeds) {
  X <- as.matrix(read.csv(sprintf("data/main/seed%d/X.csv", sd), header = FALSE))
  Y <- as.matrix(read.csv(sprintf("data/main/seed%d/Y.csv", sd), header = FALSE))
  for (meth in METHODS) {
    set.seed(sd)
    t0 <- Sys.time(); s <- run_competitor(meth, X, Y); rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    tpX <- sum(unique(s$X) <= 20); tpY <- sum(unique(s$Y) <= 60)
    puX <- if (length(s$X)) tpX / length(unique(s$X)) else 0; puY <- if (length(s$Y)) tpY / length(unique(s$Y)) else 0
    reX <- tpX / 20; reY <- tpY / 60
    rows[[length(rows) + 1]] <- data.frame(method = meth, seed = sd, F1 = (f1(puX, reX) + f1(puY, reY)) / 2,
                                           precision = (puX + puY) / 2, recall = (reX + reY) / 2, runtime_sec = rt)
  }
}
write.csv(do.call(rbind, rows), out, row.names = FALSE)
