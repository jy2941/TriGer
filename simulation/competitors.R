## Comparison methods shared by 01_main_competitors.R and 04_null_competitors.R.
## Each returns list(X = selected X indices, Y = selected Y indices), 1-based.
options(cores = 1, mc.cores = 1)
local({ ns <- asNamespace("parallel"); try(unlockBinding("detectCores", ns), silent = TRUE)
  assign("detectCores", function(...) 1L, envir = ns); try(lockBinding("detectCores", ns), silent = TRUE) })
suppressMessages({ library(mixOmics); library(PMA); library(RGCCA); library(mogsa); library(IntNMF); library(r.jive) })

nz  <- function(w, tol = 1e-8) which(abs(w) > tol)
mad3 <- function(w) which(w > median(w) + 3 * mad(w, constant = 1))

run_competitor <- function(meth, X, Y) {
  colnames(X) <- paste0("V", seq_len(ncol(X))); colnames(Y) <- paste0("V", seq_len(ncol(Y)))
  m <- ncol(X)
  if (meth == "DIABLO") {
    sp <- spls(X, Y, ncomp = 2, keepX = c(20, 20), keepY = c(60, 60), mode = "regression")
    g <- function(b) as.integer(gsub("V", "", unique(unlist(lapply(1:2, function(k) selectVar(sp, comp = k)[[b]]$name)))))
    list(X = g("X"), Y = g("Y"))
  } else if (meth == "sCCA") {
    pen <- CCA.permute(X, Y, typex = "standard", typez = "standard", nperms = 5, trace = FALSE)
    cc <- CCA(X, Y, typex = "standard", typez = "standard", penaltyx = pen$bestpenaltyx,
              penaltyz = pen$bestpenaltyz, K = 2, trace = FALSE)
    list(X = unique(c(nz(cc$u[, 1]), nz(cc$u[, 2]))), Y = unique(c(nz(cc$v[, 1]), nz(cc$v[, 2]))))
  } else if (meth == "RGCCA") {
    fit <- rgcca(blocks = list(X = X, Y = Y), method = "sgcca", sparsity = c(0.3, 0.15), ncomp = 2,
                 scheme = "factorial", verbose = FALSE)
    list(X = unique(c(nz(fit$a[[1]][, 1]), nz(fit$a[[1]][, 2]))), Y = unique(c(nz(fit$a[[2]][, 1]), nz(fit$a[[2]][, 2]))))
  } else if (meth == "moCluster") {
    mb <- mbpca(list(t(X), t(Y)), ncomp = 3, method = "globalScore", k = 0.2, option = "lambda1",
                center = TRUE, scale = FALSE, moa = TRUE, verbose = FALSE)
    ld <- mb@loading
    list(X = unique(unlist(lapply(1:ncol(ld), function(k) nz(ld[1:m, k])))),
         Y = unique(unlist(lapply(1:ncol(ld), function(k) nz(ld[-(1:m), k])))))
  } else if (meth == "IntNMF") {
    Xp <- X - min(X); Yp <- Y - min(Y)
    fit <- nmf.mnnals(dat = list(Xp / max(Xp), Yp / max(Yp)), k = 3, maxiter = 200, st.count = 10, n.ini = 3)
    list(X = mad3(apply(fit$H[[1]], 2, max)), Y = mad3(apply(fit$H[[2]], 2, max)))
  } else if (meth == "JIVE") {
    fit <- jive(list(t(X), t(Y)), method = "perm", showProgress = FALSE)
    if (fit$rankJ == 0) return(list(X = integer(0), Y = integer(0)))
    list(X = mad3(sqrt(rowSums(fit$joint[[1]]^2))), Y = mad3(sqrt(rowSums(fit$joint[[2]]^2))))
  }
}

METHODS <- c("DIABLO", "sCCA", "RGCCA", "moCluster", "IntNMF", "JIVE")
