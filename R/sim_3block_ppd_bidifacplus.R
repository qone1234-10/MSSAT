## =========================================================
## Parallel 3-block simulation: PPD + BIDIFAC+
## rep = 100, 25 cores, save each rep
## =========================================================

project_dir <- normalizePath("~/R")
setwd(project_dir)

.libPaths(c(
  file.path(project_dir, "multiview_rlibs"),
  file.path(project_dir, "PPD/multiview_rlibs"),
  .libPaths()
))

library(dplyr)
library(tidyr)
library(Matrix)
library(foreach)
library(doParallel)
library(Ckmeans.1d.dp)
library(PRIMME)
library(denoiseR)
library(RMTstat)

source("PPD/src/sim_3block_helpers_ppd_bidifacplus.R")
source("bidifac/bidifac.plus.R")
source("bidifac/bidifac.plus.given.R") 
## =========================================================
## PPD source in separated environments
## =========================================================

env2 <- new.env(parent = .GlobalEnv)
env3 <- new.env(parent = .GlobalEnv)

source("PPD/src/utils.R", local = env2)
source("PPD/src/models_2_views.R", local = env2)

source("PPD/src/utils.R", local = env3)
source("PPD/src/models_3_views.R", local = env3)

est_sigma_fixed <- function(Y) {
  sing.vals <- svd(Y, nu = 0, nv = 0)$d
  med.sing.val <- median(sing.vals)
  
  ndf_use  <- max(nrow(Y), ncol(Y))
  pdim_use <- min(nrow(Y), ncol(Y))
  
  med.mp <- RMTstat::qmp(0.5, ndf = ndf_use, pdim = pdim_use)
  
  med.sing.val / sqrt(med.mp * ndf_use)
}

env2$est.sigma <- est_sigma_fixed
env3$est.sigma <- est_sigma_fixed

ppd_2view_func <- env2$proposed_func
ppd_3view_func <- env3$proposed_func

## =========================================================
## Safe PPD 3-block decomposition
## =========================================================

run_ppd_3block_safe <- function(X1, X2, X3, r1_hat, r2_hat, r3_hat,
                                bootstrap_iter_global = 50,
                                bootstrap_iter_pair = 50) {
  n <- ncol(X1)
  I_n <- diag(n)
  
  outG <- ppd_3view_func(
    t(X1), t(X2), t(X3),
    r1_hat, r2_hat, r3_hat,
    bootstrap_iter = bootstrap_iter_global
  )
  
  PG <- outG$Pjoint
  if (is.null(PG)) PG <- matrix(0, n, n)
  
  rG <- outG$rj
  if (!is.finite(rG)) rG <- 0
  
  G1 <- X1 %*% PG
  G2 <- X2 %*% PG
  G3 <- X3 %*% PG
  
  R1 <- X1 %*% (I_n - PG)
  R2 <- X2 %*% (I_n - PG)
  R3 <- X3 %*% (I_n - PG)
  
  ## pair 12
  r1_12 <- max(r1_hat - rG, 0)
  r2_12 <- max(r2_hat - rG, 0)
  
  if (r1_12 > 0 && r2_12 > 0) {
    out12 <- ppd_2view_func(
      t(R1), t(R2),
      r1_12, r2_12,
      bootstrap_iters = bootstrap_iter_pair
    )
    P12 <- out12$Pjoint
    if (is.null(P12)) P12 <- matrix(0, n, n)
    r12_hat <- out12$rj
    if (!is.finite(r12_hat)) r12_hat <- 0
  } else {
    P12 <- matrix(0, n, n)
    r12_hat <- 0
  }
  
  J12_1 <- R1 %*% P12
  J12_2 <- R2 %*% P12
  R1 <- R1 %*% (I_n - P12)
  R2 <- R2 %*% (I_n - P12)
  
  ## pair 13
  r1_13 <- max(r1_hat - rG - r12_hat, 0)
  r3_13 <- max(r3_hat - rG, 0)
  
  if (r1_13 > 0 && r3_13 > 0) {
    out13 <- ppd_2view_func(
      t(R1), t(R3),
      r1_13, r3_13,
      bootstrap_iters = bootstrap_iter_pair
    )
    P13 <- out13$Pjoint
    if (is.null(P13)) P13 <- matrix(0, n, n)
    r13_hat <- out13$rj
    if (!is.finite(r13_hat)) r13_hat <- 0
  } else {
    P13 <- matrix(0, n, n)
    r13_hat <- 0
  }
  
  J13_1 <- R1 %*% P13
  J13_3 <- R3 %*% P13
  R1 <- R1 %*% (I_n - P13)
  R3 <- R3 %*% (I_n - P13)
  
  ## pair 23
  r2_23 <- max(r2_hat - rG - r12_hat, 0)
  r3_23 <- max(r3_hat - rG - r13_hat, 0)
  
  if (r2_23 > 0 && r3_23 > 0) {
    out23 <- ppd_2view_func(
      t(R2), t(R3),
      r2_23, r3_23,
      bootstrap_iters = bootstrap_iter_pair
    )
    P23 <- out23$Pjoint
    if (is.null(P23)) P23 <- matrix(0, n, n)
    r23_hat <- out23$rj
    if (!is.finite(r23_hat)) r23_hat <- 0
  } else {
    P23 <- matrix(0, n, n)
    r23_hat <- 0
  }
  
  J23_2 <- R2 %*% P23
  J23_3 <- R3 %*% P23
  R2 <- R2 %*% (I_n - P23)
  R3 <- R3 %*% (I_n - P23)
  
  ## individual
  rI1 <- max(r1_hat - rG - r12_hat - r13_hat, 0)
  rI2 <- max(r2_hat - rG - r12_hat - r23_hat, 0)
  rI3 <- max(r3_hat - rG - r13_hat - r23_hat, 0)
  
  I1 <- top_svd_recon(R1, rI1)
  I2 <- top_svd_recon(R2, rI2)
  I3 <- top_svd_recon(R3, rI3)
  
  list(
    ranks = c(rG, r12_hat, r13_hat, r23_hat, rI1, rI2, rI3),
    mats = list(
      G   = list(G1, G2, G3),
      J12 = list(J12_1, J12_2),
      J13 = list(J13_1, J13_3),
      J23 = list(J23_2, J23_3),
      I   = list(I1, I2, I3)
    )
  )
}

## =========================================================
## BIDIFAC+ given structure
## =========================================================

run_bidifacplus_given_3block <- function(X1, X2, X3) {
  p1 <- nrow(X1)
  p2 <- nrow(X2)
  p3 <- nrow(X3)
  n  <- ncol(X1)
  
  X0 <- rbind(X1, X2, X3)
  
  p.ind <- list()
  p.ind[[1]] <- 1:p1
  p.ind[[2]] <- (p1 + 1):(p1 + p2)
  p.ind[[3]] <- (p1 + p2 + 1):(p1 + p2 + p3)
  
  n.ind <- list()
  n.ind[[1]] <- 1:n
  
  p.ind.list <- list()
  p.ind.list[[1]] <- unlist(p.ind[c(1, 2, 3)])  # Global
  p.ind.list[[2]] <- unlist(p.ind[c(1, 2)])     # J12
  p.ind.list[[3]] <- unlist(p.ind[c(1, 3)])     # J13
  p.ind.list[[4]] <- unlist(p.ind[c(2, 3)])     # J23
  p.ind.list[[5]] <- unlist(p.ind[c(1)])        # Ind1
  p.ind.list[[6]] <- unlist(p.ind[c(2)])        # Ind2
  p.ind.list[[7]] <- unlist(p.ind[c(3)])        # Ind3
  
  n.ind.list <- vector("list", length(p.ind.list))
  for (k in seq_along(p.ind.list)) {
    n.ind.list[[k]] <- 1:n
  }
  
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
  
  G_all   <- getS(1)
  J12_all <- getS(2)
  J13_all <- getS(3)
  J23_all <- getS(4)
  I1_all  <- getS(5)
  I2_all  <- getS(6)
  I3_all  <- getS(7)
  
  blk <- function(M, k) M[p.ind[[k]], , drop = FALSE]
  
  G1 <- blk(G_all, 1)
  G2 <- blk(G_all, 2)
  G3 <- blk(G_all, 3)
  
  J12_1 <- blk(J12_all, 1)
  J12_2 <- blk(J12_all, 2)
  
  J13_1 <- blk(J13_all, 1)
  J13_3 <- blk(J13_all, 3)
  
  J23_2 <- blk(J23_all, 2)
  J23_3 <- blk(J23_all, 3)
  
  I1 <- blk(I1_all, 1)
  I2 <- blk(I2_all, 2)
  I3 <- blk(I3_all, 3)
  
  rG <- round(mean(c(rank_num(G1), rank_num(G2), rank_num(G3))))
  r12_hat <- round(mean(c(rank_num(J12_1), rank_num(J12_2))))
  r13_hat <- round(mean(c(rank_num(J13_1), rank_num(J13_3))))
  r23_hat <- round(mean(c(rank_num(J23_2), rank_num(J23_3))))
  rI1 <- rank_num(I1)
  rI2 <- rank_num(I2)
  rI3 <- rank_num(I3)
  
  list(
    ranks = c(rG, r12_hat, r13_hat, r23_hat, rI1, rI2, rI3),
    mats = list(
      G   = list(G1, G2, G3),
      J12 = list(J12_1, J12_2),
      J13 = list(J13_1, J13_3),
      J23 = list(J23_2, J23_3),
      I   = list(I1, I2, I3)
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
p3 <- 200

r1 <- 6
r2 <- 9
r3 <- 12

r <- 3
r12 <- 1
r13 <- 1
r23 <- 2

d <- 1
ang <- pi / 8
ang_label <- "pi8"

d1 <- d/2 * 3 * sqrt(n) * c(6,5,4,7,3,2)
d2 <- d/2 * 2 * sqrt(n) * c(7,6,5,8,9,10,2,3,4)
d3 <- d/2 * 4 * sqrt(n) * c(10,9,8,7,5,6,2,3,4,11,12,13)


s1 <- 1
s2 <- 1
s3 <- 1

set.seed(0)

v0 <- random_orthonormal(n, r)

v_temp <- random_orthonormal_orthogonal(r12 + r13 + r23, v0, n)

v12 <- v_temp[, 1:r12, drop = FALSE]
v13 <- v_temp[, (r12 + 1):(r12 + r13), drop = FALSE]
v23 <- v_temp[, (r12 + r13 + 1):(r12 + r13 + r23), drop = FALSE]

v11 <- random_orthonormal_orthogonal(
  r1 - r - r12 - r13,
  cbind(v0, v12, v13),
  n
)

v22 <- random_orthonormal_orthogonal(
  r2 - r - r12 - r23,
  cbind(v0, v12, v23),
  n
)

v33_temp <- random_orthonormal_orthogonal(
  r3 - r - r13 - r23,
  cbind(v0, v13, v23, v22),
  n
)

v33 <- cbind(
  v22 * cos(ang) +
    v33_temp[, 1:(r2 - r - r12 - r23), drop = FALSE] * sin(ang),
  v33_temp[, ((r2 - r - r12 - r23) + 1):(r3 - r - r13 - r23), drop = FALSE]
)

u1 <- random_orthonormal(p1, r1)
u2 <- random_orthonormal(p2, r2)
u3 <- random_orthonormal(p3, r3)

true_ranks <- c(
  r, r12, r13, r23,
  r1 - (r + r12 + r13),
  r2 - (r + r12 + r23),
  r3 - (r + r13 + r23)
)

## =========================================================
## parallel setting
## =========================================================

rep <- 100
ncores <- 25

save_dir <- file.path(project_dir, "ppd,bidifacplus.3block.simul")
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

cl <- makeCluster(ncores)
registerDoParallel(cl)

clusterExport(cl, c(
  "project_dir",
  "save_dir",
  
  "random_orthonormal",
  "random_orthonormal_orthogonal",
  "rel_err",
  "make_comp",
  "rank_num",
  "top_svd_recon",
  "estimate_rank",
  
  "ppd_2view_func",
  "ppd_3view_func",
  "run_ppd_3block_safe",
  "run_bidifacplus_given_3block",
  
  "bidifac.plus.given",
  
  "n", "p1", "p2", "p3",
  "r1", "r2", "r3",
  "r", "r12", "r13", "r23",
  "d", "ang",
  "d1", "d2", "d3",
  "s1", "s2", "s3",
  "v0", "v12", "v13", "v23", "v11", "v22", "v33",
  "u1", "u2", "u3",
  "true_ranks"
))

clusterEvalQ(cl, {
  setwd(project_dir)
  
  .libPaths(c(
    file.path(project_dir, "multiview_rlibs"),
    file.path(project_dir, "PPD/multiview_rlibs"),
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

methods <- c("PPD", "BIDIFAC+")

iters <- foreach(
  rep_idx = 1:rep,
  .packages = c("dplyr", "Matrix", "Ckmeans.1d.dp", "PRIMME", "denoiseR", "RMTstat")
) %dopar% {
  
  out_file <- file.path(save_dir, sprintf("rep_%03d_%s.rds", rep_idx, ang_label))
  
  set.seed(10000 + rep_idx)
  
  ## =====================================================
  ## Generate one dataset per rep
  ## =====================================================
  rotG2  <- random_orthonormal(r, r)
  rotG3  <- random_orthonormal(r, r)
  rot12  <- random_orthonormal(r12, r12)
  rot13  <- random_orthonormal(r13, r13)
  rot23  <- random_orthonormal(r23, r23)
  
  basis1 <- cbind(v0, v12, v13, v11)
  basis2 <- cbind(v0 %*% rotG2, v12 %*% rot12, v23, v22)
  basis3 <- cbind(v0 %*% rotG3, v13 %*% rot13, v23 %*% rot23, v33)
  
  S1 <- u1 %*% diag(d1, nrow = r1) %*% t(basis1)
  S2 <- u2 %*% diag(d2, nrow = r2) %*% t(basis2)
  S3 <- u3 %*% diag(d3, nrow = r3) %*% t(basis3)
  
  X1 <- S1 + s1 * matrix(rnorm(p1 * n), p1, n)
  X2 <- S2 + s2 * matrix(rnorm(p2 * n), p2, n)
  X3 <- S3 + s3 * matrix(rnorm(p3 * n), p3, n)
  
  idx0_1  <- 1:r
  idx12_1 <- (r + 1):(r + r12)
  idx13_1 <- (r + r12 + 1):(r + r12 + r13)
  idxI1_1 <- (r + r12 + r13 + 1):r1
  
  idx0_2  <- 1:r
  idx12_2 <- (r + 1):(r + r12)
  idx23_2 <- (r + r12 + 1):(r + r12 + r23)
  idxI2_2 <- (r + r12 + r23 + 1):r2
  
  idx0_3  <- 1:r
  idx13_3 <- (r + 1):(r + r13)
  idx23_3 <- (r + r13 + 1):(r + r13 + r23)
  idxI3_3 <- (r + r13 + r23 + 1):r3
  
  trueG1 <- make_comp(u1, d1, basis1, idx0_1)
  trueG2 <- make_comp(u2, d2, basis2, idx0_2)
  trueG3 <- make_comp(u3, d3, basis3, idx0_3)
  
  true12_1 <- make_comp(u1, d1, basis1, idx12_1)
  true12_2 <- make_comp(u2, d2, basis2, idx12_2)
  
  true13_1 <- make_comp(u1, d1, basis1, idx13_1)
  true13_3 <- make_comp(u3, d3, basis3, idx13_3)
  
  true23_2 <- make_comp(u2, d2, basis2, idx23_2)
  true23_3 <- make_comp(u3, d3, basis3, idx23_3)
  
  trueI1 <- make_comp(u1, d1, basis1, idxI1_1)
  trueI2 <- make_comp(u2, d2, basis2, idxI2_2)
  trueI3 <- make_comp(u3, d3, basis3, idxI3_3)
  
  r1_hat <- estimate_rank(X1, r1)
  r2_hat <- estimate_rank(X2, r2)
  r3_hat <- estimate_rank(X3, r3)
  
  run_one_method <- function(method) {
    t0 <- proc.time()[3]
    
    tryCatch({
      if (method == "PPD") {
        fit <- run_ppd_3block_safe(
          X1, X2, X3,
          r1_hat, r2_hat, r3_hat,
          bootstrap_iter_global = 50,
          bootstrap_iter_pair = 50
        )
      } else {
        fit <- run_bidifacplus_given_3block(X1, X2, X3)
      }
      
      eG <- mean(c(
        rel_err(trueG1, fit$mats$G[[1]]),
        rel_err(trueG2, fit$mats$G[[2]]),
        rel_err(trueG3, fit$mats$G[[3]])
      ))
      
      e12 <- mean(c(
        rel_err(true12_1, fit$mats$J12[[1]]),
        rel_err(true12_2, fit$mats$J12[[2]])
      ))
      
      e13 <- mean(c(
        rel_err(true13_1, fit$mats$J13[[1]]),
        rel_err(true13_3, fit$mats$J13[[2]])
      ))
      
      e23 <- mean(c(
        rel_err(true23_2, fit$mats$J23[[1]]),
        rel_err(true23_3, fit$mats$J23[[2]])
      ))
      
      eI1 <- rel_err(trueI1, fit$mats$I[[1]])
      eI2 <- rel_err(trueI2, fit$mats$I[[2]])
      eI3 <- rel_err(trueI3, fit$mats$I[[3]])
      
      tibble(
        rep = rep_idx,
        method = method,
        rG = fit$ranks[1],
        r12 = fit$ranks[2],
        r13 = fit$ranks[3],
        r23 = fit$ranks[4],
        rI1 = fit$ranks[5],
        rI2 = fit$ranks[6],
        rI3 = fit$ranks[7],
        eG = eG,
        e12 = e12,
        e13 = e13,
        e23 = e23,
        eI1 = eI1,
        eI2 = eI2,
        eI3 = eI3,
        time = proc.time()[3] - t0,
        error = NA_character_
      )
    }, error = function(e) {
      tibble(
        rep = rep_idx,
        method = method,
        rG = NA_real_, r12 = NA_real_, r13 = NA_real_, r23 = NA_real_,
        rI1 = NA_real_, rI2 = NA_real_, rI3 = NA_real_,
        eG = NA_real_, e12 = NA_real_, e13 = NA_real_, e23 = NA_real_,
        eI1 = NA_real_, eI2 = NA_real_, eI3 = NA_real_,
        time = proc.time()[3] - t0,
        error = e$message
      )
    })
  }
  
  rep_df <- dplyr::bind_rows(
    run_one_method("PPD"),
    run_one_method("BIDIFAC+")
  )
  
  saveRDS(rep_df, out_file)
  rep_df
}

stopCluster(cl)

raw_df <- dplyr::bind_rows(iters)

write.csv(
  raw_df,
  file.path(save_dir, "raw_results_3block_ppd_bidifacplus.csv"),
  row.names = FALSE
)

cat("Saved raw results:\n")
cat(file.path(save_dir, "raw_results_3block_ppd_bidifacplus.csv"), "\n")