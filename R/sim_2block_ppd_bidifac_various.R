## =========================================================
## Parallel 2-block grid simulation: PPD + BIDIFAC+
## d x angle grid, rep = 100
## save each rep/method/setting separately
## =========================================================

project_dir <- normalizePath("~/R")
setwd(project_dir)

.libPaths(c(
  file.path(project_dir, "multiview_rlibs"),
  file.path(project_dir, "PPD/multiview_rlibs"),
  .libPaths()
))

library(dplyr)
library(foreach)
library(doParallel)
library(Matrix)
library(Ckmeans.1d.dp)
library(PRIMME)
library(denoiseR)
library(RMTstat)

source("PPD/src/sim_3block_helpers_ppd_bidifacplus.R")
source("bidifac/bidifac.plus.R")
source("bidifac/bidifac.plus.given.R")

## =========================================================
## PPD source in isolated env
## =========================================================

env2 <- new.env(parent = .GlobalEnv)

source("PPD/src/utils.R", local = env2)
source("PPD/src/models_2_views.R", local = env2)

est_sigma_fixed <- function(Y) {
  sing.vals <- svd(Y, nu = 0, nv = 0)$d
  med.sing.val <- median(sing.vals)
  
  ndf_use  <- max(nrow(Y), ncol(Y))
  pdim_use <- min(nrow(Y), ncol(Y))
  
  med.mp <- RMTstat::qmp(
    0.5,
    ndf = ndf_use,
    pdim = pdim_use
  )
  
  med.sing.val / sqrt(med.mp * ndf_use)
}

env2$est.sigma <- est_sigma_fixed
ppd_2view_func <- env2$proposed_func

## =========================================================
## helpers
## =========================================================

random_orthonormal <- function(n, r) {
  if (r <= 0) return(matrix(0, n, 0))
  qr.Q(qr(matrix(rnorm(n * r), n, r)))[, seq_len(r), drop = FALSE]
}

random_orthonormal_orthogonal <- function(r_new, V_old, n = nrow(V_old)) {
  if (r_new <= 0) return(matrix(0, n, 0))
  Z <- matrix(rnorm(n * r_new), n, r_new)
  Z <- Z - V_old %*% crossprod(V_old, Z)
  qr.Q(qr(Z))[, seq_len(r_new), drop = FALSE]
}

rel_err <- function(A, B) {
  sum((A - B)^2) / max(sum(A^2), .Machine$double.eps)
}

rank_num <- function(X, tol = 1e-6) {
  if (is.null(X) || length(X) == 0) return(0L)
  dd <- svd(X, nu = 0, nv = 0)$d
  sum(dd > tol)
}

top_svd_recon <- function(X, r) {
  if (r <= 0) return(matrix(0, nrow(X), ncol(X)))
  ss <- svd(X)
  ss$u[, seq_len(r), drop = FALSE] %*%
    diag(ss$d[seq_len(r)], nrow = r) %*%
    t(ss$v[, seq_len(r), drop = FALSE])
}

estimate_rank <- function(X, true_rank = NULL) {
  if (exists("BEMA_combined")) {
    out <- BEMA_combined(X)
    return(as.integer(out[[2]]))
  }
  if (!is.null(true_rank)) return(as.integer(true_rank))
  stop("No BEMA_combined found and no true_rank given.")
}

run_ppd_2block <- function(X1, X2, r1_hat, r2_hat, bootstrap_iters = 50) {
  n <- ncol(X1)
  
  fit <- ppd_2view_func(
    t(X1), t(X2),
    r1_hat, r2_hat,
    bootstrap_iters = bootstrap_iters
  )
  
  Pj <- fit$Pjoint
  if (is.null(Pj)) Pj <- matrix(0, n, n)
  
  Pind1 <- fit$Pindiv1
  Pind2 <- fit$Pindiv2
  
  if (is.null(Pind1)) Pind1 <- matrix(0, n, n)
  if (is.null(Pind2)) Pind2 <- matrix(0, n, n)
  
  rJ <- fit$rj
  if (!is.finite(rJ)) rJ <- 0
  
  rI1 <- if (!is.null(fit$ri1) && is.finite(fit$ri1)) fit$ri1 else max(r1_hat - rJ, 0)
  rI2 <- if (!is.null(fit$ri2) && is.finite(fit$ri2)) fit$ri2 else max(r2_hat - rJ, 0)
  
  list(
    ranks = c(rJ, rI1, rI2),
    mats = list(
      J = list(X1 %*% Pj, X2 %*% Pj),
      I = list(X1 %*% Pind1, X2 %*% Pind2)
    ),
    fit = fit
  )
}

run_bidifacplus_given_2block <- function(X1, X2) {
  p1 <- nrow(X1)
  p2 <- nrow(X2)
  n  <- ncol(X1)
  
  X0 <- rbind(X1, X2)
  
  p.ind <- list(
    1:p1,
    (p1 + 1):(p1 + p2)
  )
  
  n.ind <- list(1:n)
  
  p.ind.list <- list(
    unlist(p.ind[c(1, 2)]),  # joint
    unlist(p.ind[1]),        # individual 1
    unlist(p.ind[2])         # individual 2
  )
  
  n.ind.list <- rep(list(1:n), length(p.ind.list))
  
  fit <- bidifac.plus.given(
    X0 = X0,
    p.ind = p.ind,
    n.ind = n.ind,
    p.ind.list = p.ind.list,
    n.ind.list = n.ind.list,
    max.iter = 500,
    conv.thresh = 1e-3
  )
  
  S <- fit$S
  zero_all <- matrix(0, nrow(X0), n)
  
  getS <- function(i) {
    if (length(S) < i || is.null(S[[i]])) return(zero_all)
    S[[i]]
  }
  
  J_all  <- getS(1)
  I1_all <- getS(2)
  I2_all <- getS(3)
  
  J1 <- J_all[p.ind[[1]], , drop = FALSE]
  J2 <- J_all[p.ind[[2]], , drop = FALSE]
  
  I1 <- I1_all[p.ind[[1]], , drop = FALSE]
  I2 <- I2_all[p.ind[[2]], , drop = FALSE]
  
  rJ  <- round(mean(c(rank_num(J1), rank_num(J2))))
  rI1 <- rank_num(I1)
  rI2 <- rank_num(I2)
  
  list(
    ranks = c(rJ, rI1, rI2),
    mats = list(
      J = list(J1, J2),
      I = list(I1, I2)
    ),
    fit = fit
  )
}

## =========================================================
## simulation setting
## =========================================================

n <- 500
p1 <- 400
p2 <- 300

r1 <- 6
r2 <- 9
r  <- 3

s1 <- 1
s2 <- 1

d_vals <- c(0.1, 0.25, 0.5, 1, 2)
ang_vals <- c(0, pi / 8, pi / 4, 3 * pi / 8)

rep <- 100
ncores <- 25

methods <- c("PPD", "BIDIFAC+")

true_ranks <- c(rJ = r, rI1 = r1 - r, rI2 = r2 - r)

set.seed(0)

u1 <- random_orthonormal(p1, r1)
u2 <- random_orthonormal(p2, r2)

save_dir <- file.path(project_dir, "1 ppd,bidifacplus.simul")
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## job grid
## =========================================================

job_grid <- expand.grid(
  rep_idx = seq_len(rep),
  d = d_vals,
  ang = ang_vals,
  method = methods,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

job_grid$job_id <- seq_len(nrow(job_grid))

cat("Total jobs:", nrow(job_grid), "\n")
cat("Save dir:", save_dir, "\n")

## =========================================================
## cluster
## =========================================================

cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterExport(cl, c(
  "project_dir",
  "save_dir",
  
  "random_orthonormal",
  "random_orthonormal_orthogonal",
  "rel_err",
  "rank_num",
  "top_svd_recon",
  "estimate_rank",
  
  "ppd_2view_func",
  "run_ppd_2block",
  "run_bidifacplus_given_2block",
  "bidifac.plus.given",
  
  "n", "p1", "p2",
  "r1", "r2", "r",
  "s1", "s2",
  "u1", "u2",
  "true_ranks"
))

clusterEvalQ(cl, {
  setwd(project_dir)
  
  .libPaths(c(
    file.path(project_dir, "1 multiview_rlibs"),
    file.path(project_dir, "1 PPD/1 multiview_rlibs"),
    .libPaths()
  ))
  
  library(dplyr)
  library(Matrix)
  library(Ckmeans.1d.dp)
  library(PRIMME)
  library(denoiseR)
  library(RMTstat)
  
  NULL
})

## =========================================================
## run
## =========================================================

iters <- foreach(
  jj = seq_len(nrow(job_grid)),
  .packages = c("dplyr", "Matrix", "Ckmeans.1d.dp", "PRIMME", "denoiseR", "RMTstat")
) %dopar% {
  
  job <- job_grid[jj, ]
  
  rep_idx <- job$rep_idx
  d       <- job$d
  ang     <- job$ang
  method  <- job$method
  
  method_tag <- ifelse(method == "BIDIFAC+", "BIDIFACplus", "PPD")
  
  out_file <- file.path(
    save_dir,
    sprintf(
      "rep_%03d_d%.2f_ang%.3f_%s.rds",
      rep_idx,
      d,
      ang,
      method_tag
    )
  )
  
  if (file.exists(out_file)) {
    return(readRDS(out_file))
  }
  
  t0 <- proc.time()[3]
  
  out <- tryCatch({
    
    set.seed(1000000 + 10000 * rep_idx + 1000 * match(d, c(0.1, 0.25, 0.5, 1, 2)) +
               10 * match(ang, c(0, pi / 8, pi / 4, 3 * pi / 8)))
    
    d1 <- d * 3 * sqrt(n) * c(6,5,4,2,3,7)
    d2 <- d * 2 * sqrt(n) * c(7,6,5,2,3,4,8,9,10)
    
    v0 <- random_orthonormal(n, r)
    
    v1 <- random_orthonormal_orthogonal(
      r1 - r,
      v0,
      n
    )
    
    v_tmp <- random_orthonormal_orthogonal(
      r2 - r,
      cbind(v0, v1),
      n
    )
    
    v2 <- cbind(
      v_tmp[, seq_len(r1 - r), drop = FALSE] * sin(ang) +
        v1 * cos(ang),
      v_tmp[, ((r1 - r) + 1):(r2 - r), drop = FALSE]
    )
    
    rot <- random_orthonormal(r, r)
    
    basis1 <- cbind(v0, v1)
    basis2 <- cbind(v0 %*% rot, v2)
    
    X1 <- u1 %*% diag(d1, nrow = r1) %*% t(basis1) +
      s1 * matrix(rnorm(p1 * n), p1, n)
    
    X2 <- u2 %*% diag(d2, nrow = r2) %*% t(basis2) +
      s2 * matrix(rnorm(p2 * n), p2, n)
    
    # true_j1 <- u1 %*% diag(c(rep(0, r), d1[(r + 1):r1]), nrow = r1) %*% t(basis1)
    # true_j2 <- u2 %*% diag(c(d2[1:r], rep(0, r2 - r)), nrow = r2) %*% t(basis2)
    # 
    # true_i1 <- u1 %*% diag(c(d1[1:r], rep(0, r1 - r)), nrow = r1) %*% t(basis1)
    # true_i2 <- u2 %*% diag(c(rep(0, r), d2[(r + 1):r2]), nrow = r2) %*% t(basis2)
    
    true_j1 <- u1 %*% diag(c(d1[1:r], rep(0, r1 - r)), nrow = r1) %*% t(basis1)
    true_j2 <- u2 %*% diag(c(d2[1:r], rep(0, r2 - r)), nrow = r2) %*% t(basis2)
    
    true_i1 <- u1 %*% diag(c(rep(0, r), d1[(r + 1):r1]), nrow = r1) %*% t(basis1)
    true_i2 <- u2 %*% diag(c(rep(0, r), d2[(r + 1):r2]), nrow = r2) %*% t(basis2)
    
    r1_hat <- estimate_rank(X1, r1)
    r2_hat <- estimate_rank(X2, r2)
    
    if (method == "PPD") {
      fit <- run_ppd_2block(
        X1, X2,
        r1_hat, r2_hat,
        bootstrap_iters = 50
      )
    } else {
      fit <- run_bidifacplus_given_2block(X1, X2)
    }
    
    tibble(
      rep = rep_idx,
      d = d,
      ang = ang,
      ang_deg = ang * 180 / pi,
      method = method,
      
      rJ = fit$ranks[1],
      rI1 = fit$ranks[2],
      rI2 = fit$ranks[3],
      
      eJ1 = rel_err(true_j1, fit$mats$J[[1]]),
      eJ2 = rel_err(true_j2, fit$mats$J[[2]]),
      eI1 = rel_err(true_i1, fit$mats$I[[1]]),
      eI2 = rel_err(true_i2, fit$mats$I[[2]]),
      
      time = proc.time()[3] - t0,
      error = NA_character_
    )
    
  }, error = function(e) {
    
    tibble(
      rep = rep_idx,
      d = d,
      ang = ang,
      ang_deg = ang * 180 / pi,
      method = method,
      
      rJ = NA_real_,
      rI1 = NA_real_,
      rI2 = NA_real_,
      
      eJ1 = NA_real_,
      eJ2 = NA_real_,
      eI1 = NA_real_,
      eI2 = NA_real_,
      
      time = proc.time()[3] - t0,
      error = e$message
    )
  })
  
  saveRDS(out, out_file)
  
  out
}

stopCluster(cl)

raw_df <- dplyr::bind_rows(iters)

write.csv(
  raw_df,
  file.path(save_dir, "raw_results_grid_ppd_bidifacplus.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    raw = raw_df,
    true_ranks = true_ranks,
    d_vals = d_vals,
    ang_vals = ang_vals
  ),
  file.path(save_dir, "summary_grid_ppd_bidifacplus.rds")
)

cat("Saved raw CSV:\n")
cat(file.path(save_dir, "raw_results_grid_ppd_bidifacplus.csv"), "\n")