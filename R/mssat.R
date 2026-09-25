# Multi-block joint-direction scan, R translation of the uploaded MATLAB mssat.m.
# Input convention: each X_list[[k]] is a p_k x n matrix, with common column size n.
# If rank_method = "bema" or sigma_method = "bema", a function BEMA(data, alpha = 0.2)
# must already exist in the R environment. BEMA is called as BEMA(t(X), alpha = 0.2).

mssat <- function(X_list,
                  center = TRUE,
                  alpha = 0.05,
                  verbose = TRUE,
                  tol_proj = 1e-10,
                  rank_method = c("bema", "mp", "given"),
                  sigma_method = c("bema", "mp", "given"),
                  r_init = NULL,
                  sigma_init = NULL,
                  test_mode = c("alpha", "alpha_seq", "altmean"),
                  alt_angle_deg = 10,
                  mode = c("test", "all")) {

  rank_method <- match.arg(rank_method)
  sigma_method <- match.arg(sigma_method)
  test_mode <- match.arg(test_mode)
  mode <- match.arg(mode)

  if (!is.list(X_list)) stop("X_list must be a list of matrices.")
  X_list <- lapply(X_list, as.matrix)
  K <- length(X_list)
  if (K < 2) stop("K >= 2 required.")

  n <- ncol(X_list[[1]])
  if (is.null(n) || n <= 0) stop("Each block must be a p_k x n matrix.")
  for (k in seq_len(K)) {
    if (ncol(X_list[[k]]) != n) stop("All blocks must share the same n, i.e. the number of columns.")
  }

  if (!is.logical(center) && !(center %in% c(0, 1))) stop("center must be TRUE or FALSE.")
  center <- as.logical(center)
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a scalar in (0, 1).")
  }
  if (!is.numeric(tol_proj) || length(tol_proj) != 1 || tol_proj < 0) {
    stop("tol_proj must be a nonnegative scalar.")
  }
  if (!is.numeric(alt_angle_deg) || length(alt_angle_deg) != 1 || alt_angle_deg < 0 || alt_angle_deg >= 180) {
    stop("alt_angle_deg must be in [0, 180).")
  }
  if (rank_method == "given" && (is.null(r_init) || length(r_init) != K)) {
    stop("r_init must have length K when rank_method = 'given'.")
  }
  if (sigma_method == "given" && (is.null(sigma_init) || length(sigma_init) != K)) {
    stop("sigma_init must have length K when sigma_method = 'given'.")
  }
  if ((rank_method == "bema" || sigma_method == "bema") && !exists("BEMA", mode = "function")) {
    stop("BEMA(data, alpha) must exist when rank_method or sigma_method is 'bema'.")
  }

  zero_mat <- function(nr, nc) matrix(0, nrow = nr, ncol = nc)

  if (center) {
    X_list <- lapply(X_list, function(X) sweep(X, 1, rowMeans(X), FUN = "-"))
  }

  mp_sigma_rank <- function(X) {
    nk <- ncol(X)
    pk <- nrow(X)
    S <- crossprod(X) / nk
    S <- (S + t(S)) / 2
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    ev <- pmax(ev, 0)
    y <- pk / nk
    sigma2_hat <- median(ev)
    edge <- sigma2_hat * (1 + sqrt(y))^2
    r_mp <- sum(ev > edge)
    list(sigma2 = sigma2_hat, r = r_mp)
  }

  bema_sigma_rank <- function(X) {
    nk <- ncol(X)
    pk <- nrow(X)
    sigma2_hat <- BEMA(t(X), alpha = 0.2)
    sigma2_hat <- max(as.numeric(sigma2_hat), 0)
    S <- crossprod(X) / nk
    S <- (S + t(S)) / 2
    ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    ev <- pmax(ev, 0)
    y <- pk / nk
    edge <- sigma2_hat * (1 + sqrt(y))^2
    r_hat <- sum(ev > edge)
    list(sigma2 = sigma2_hat, r = r_hat)
  }

  est_block_init <- function(X, kidx) {
    sigma2_hat <- NULL
    r_hat <- NULL

    if (sigma_method == "mp" && rank_method == "mp") {
      tmp <- mp_sigma_rank(X)
      sigma2_hat <- tmp$sigma2
      r_hat <- tmp$r
    } else if (sigma_method == "mp" && rank_method != "mp") {
      tmp <- mp_sigma_rank(X)
      sigma2_hat <- tmp$sigma2
    } else if (sigma_method != "mp" && rank_method == "mp") {
      tmp <- mp_sigma_rank(X)
      r_hat <- tmp$r
    } else if (sigma_method == "bema" && rank_method == "bema") {
      tmp <- bema_sigma_rank(X)
      sigma2_hat <- tmp$sigma2
      r_hat <- tmp$r
    } else if (sigma_method == "bema" && rank_method != "bema") {
      tmp <- bema_sigma_rank(X)
      sigma2_hat <- tmp$sigma2
    } else if (sigma_method != "bema" && rank_method == "bema") {
      tmp <- bema_sigma_rank(X)
      r_hat <- tmp$r
    }

    if (sigma_method == "given") {
      sigma2_hat <- max(as.numeric(sigma_init[kidx]), 0)
    } else {
      sigma2_hat <- max(as.numeric(sigma2_hat), 0)
    }

    if (rank_method == "given") {
      r_hat <- max(as.numeric(r_init[kidx]), 0)
    } else {
      r_hat <- max(as.numeric(r_hat), 0)
    }

    if (abs(r_hat - round(r_hat)) > .Machine$double.eps^0.5) {
      warning(sprintf("Initial rank for block %d is non-integer; using floor(r).", kidx))
    }
    r_hat <- as.integer(floor(r_hat))

    if (r_hat > min(dim(X))) {
      stop(sprintf("Initial rank for block %d exceeds min(dim(X)).", kidx))
    }

    list(sigma2 = sigma2_hat, r = r_hat, y = nrow(X) / ncol(X))
  }

  blks <- lapply(seq_len(K), function(k) est_block_init(X_list[[k]], k))
  rhat_vec <- vapply(blks, function(b) b$r, numeric(1))

  if (verbose) {
    cat(sprintf("Initial ranks (%s): %s\n", rank_method,
                paste(sprintf("r%d=%.0f", seq_len(K), rhat_vec), collapse = ", ")))
    cat(sprintf("Initial sigma^2 (%s): %s\n", sigma_method,
                paste(sprintf("%.4f", vapply(blks, function(b) b$sigma2, numeric(1))), collapse = ", ")))
  }

  VA_list <- vector("list", K)
  for (k in seq_len(K)) {
    rk <- rhat_vec[k]
    if (rk > 0) {
      sv <- svd(X_list[[k]], nu = 0, nv = rk)
      VA_list[[k]] <- sv$v[, seq_len(rk), drop = FALSE]
    } else {
      VA_list[[k]] <- zero_mat(n, 0)
    }
  }

  ## Exact inversion used in MSSAT Eq. (11)-(12):
  ## z = s^2 / (n * sigma^2), d_hat^2 = n*sigma^2/2 * [z-1-y+sqrt((z-1-y)^2-4y)],
  ## and lambda_hat = d_hat / (sqrt(n) * sigma).  pmax() only guards finite-sample
  ## roundoff / sub-edge values; it is not an additional approximation.
  inv_bgn <- function(s2, sigma2, y) {
    if (!is.finite(sigma2) || sigma2 <= 0) return(rep(0, length(s2)))
    z <- s2 / (n * sigma2)
    disc <- pmax((z - 1 - y)^2 - 4 * y, 0)
    (n * sigma2 / 2) * (z - 1 - y + sqrt(disc))
  }
  ell_hat <- function(s, sigma2, y) {
    if (!is.finite(sigma2) || sigma2 <= 0) return(rep(0, length(s)))
    pmax(sqrt(pmax(inv_bgn(s^2, sigma2, y), 0)) / (sqrt(n) * sqrt(sigma2)), 0)
  }
  a_y <- function(ell, y) {
    (ell^4 - y) / (ell^2 * (ell^2 + 1))
  }
  theta_y <- function(ell, y) {
    (ell^4 + 2 * y * ell^2 + y) / (ell^3 * (ell^2 + 1)^2)
  }
  psi_y <- function(ell, y) {
    (ell^6 - 3 * y * ell^2 - 2 * y) / (ell^3 * (ell^2 + 1)^2)
  }
  V_E_y <- function(ell, y) {
    (2 / (ell^4 - y)) * (
      2 * y * (y + 1) * theta_y(ell, y)^2 -
        (y * (y - 1) * (5 * y + 1)) / (ell * (ell^2 + 1)^2) * theta_y(ell, y) +
        ((ell^4 + y) * (ell^2 + y)^2) / (ell^3 * (ell^2 + 1)^2) * psi_y(ell, y) +
        2 * y^2 * (y - 1)^2 / (ell^2 * (ell^2 + 1)^4)
    )
  }
  V_y <- function(ell, y) {
    (4 * theta_y(ell, y)^2 + V_E_y(ell, y)) / (4 * a_y(ell, y))
  }

  M <- do.call(cbind, VA_list)
  smax <- min(rhat_vec)

  pack_output <- function(rJ, VJ, PJ, acc_s, hist, rec) {
    list(
      joint_rank = as.numeric(rJ),
      V_joint = VJ,
      P_joint = PJ,
      accepted_s = as.numeric(acc_s),
      test_history = hist,
      blocks = rec,
      options = list(
        center = center,
        alpha = alpha,
        tol_proj = tol_proj,
        rank_method = rank_method,
        sigma_method = sigma_method,
        r_init = r_init,
        sigma_init = sigma_init,
        test_mode = test_mode,
        alt_angle_deg = alt_angle_deg,
        mode = mode
      )
    )
  }

  if (smax == 0) {
    if (verbose) cat("No singular vectors to test (smax=0).\n")

    V_joint <- zero_mat(n, 0)
    P_joint <- zero_mat(n, n)
    rec <- vector("list", K)
    V_block_joint <- vector("list", K)
    V_block_indiv <- vector("list", K)
    U_block_joint <- vector("list", K)
    U_block_indiv <- vector("list", K)

    for (k in seq_len(K)) {
      pk <- nrow(X_list[[k]])
      Vk0 <- VA_list[[k]]
      PVAk <- if (ncol(Vk0) > 0) Vk0 %*% t(Vk0) else zero_mat(n, n)
      V_block_joint[[k]] <- zero_mat(n, 0)
      V_block_indiv[[k]] <- Vk0
      U_block_joint[[k]] <- zero_mat(pk, 0)
      U_block_indiv[[k]] <- zero_mat(pk, ncol(Vk0))
      rec[[k]] <- list(
        P_VAk = PVAk,
        P_VAk_given_VJ = zero_mat(n, n),
        P_indiv = PVAk,
        V_joint_block = V_block_joint[[k]],
        V_indiv = V_block_indiv[[k]],
        U_joint_block = U_block_joint[[k]],
        U_indiv = U_block_indiv[[k]]
      )
    }

    out <- pack_output(0, V_joint, P_joint, numeric(0), list(), rec)
    out$V_block_joint <- V_block_joint
    out$V_block_indiv <- V_block_indiv
    out$U_block_joint <- U_block_joint
    out$U_block_indiv <- U_block_indiv
    out$test_path <- list(T2 = numeric(0), crit = numeric(0), mu0 = numeric(0),
                          mu_alt = numeric(0), V0 = numeric(0), Z = numeric(0),
                          accept = logical(0), rank_accept = logical(0))
    out$T2_path <- numeric(0)
    out$crit_path <- numeric(0)
    out$accept_path <- logical(0)
    out$first_reject_s <- NA_real_
    out$rank_accept_path <- logical(0)
    return(out)
  }

  svM <- svd(M, nu = smax, nv = 0)
  U_M <- svM$u[, seq_len(smax), drop = FALSE]

  T2_path <- rep(NA_real_, smax)
  crit_path <- rep(NA_real_, smax)
  mu0_path <- rep(NA_real_, smax)
  muAlt_path <- rep(NA_real_, smax)
  V0_path <- rep(NA_real_, smax)
  Z_path <- rep(NA_real_, smax)
  accept_path <- rep(FALSE, smax)

  test_history <- vector("list", 0)
  VJ_list <- vector("list", 0)
  accepted_s <- numeric(0)

  prefix_alive <- TRUE
  first_reject_s <- NA_real_
  rank_accept_path <- rep(FALSE, smax)
  tiny <- 1e-12

  for (s in seq_len(smax)) {
    vJ_hat <- U_M[, s, drop = FALSE]

    vAkJ <- vector("list", K)
    m_hat <- rep(0, K)
    V_hat <- rep(0, K)
    a_tilde <- rep(0, K)
    V_tilde <- rep(0, K)

    for (k in seq_len(K)) {
      Vk <- VA_list[[k]]
      rk <- ncol(Vk)

      if (rk == 0) {
        vAkJ[[k]] <- zero_mat(n, 1)
        next
      }

      Xk <- X_list[[k]]
      sig2 <- blks[[k]]$sigma2
      yk <- blks[[k]]$y

      w_raw <- crossprod(Vk, vJ_hat)
      nw <- sqrt(sum(w_raw^2))
      if (!is.finite(nw) || nw <= tiny) {
        vAkJ[[k]] <- zero_mat(n, 1)
        next
      }
      w_hat <- as.numeric(w_raw / nw)

      XkVk <- Xk %*% Vk
      sdir_i <- sqrt(colSums(XkVk^2))
      ell_i <- ell_hat(sdir_i, sig2, yk)

      a_i <- pmax(a_y(ell_i, yk), 0)
      V_i <- pmax(V_y(ell_i, yk), 0)
      m_i <- sqrt(pmax(a_i, 0))

      valid <- is.finite(m_i) & is.finite(V_i) & (m_i > tiny)
      if (!any(valid)) {
        vAkJ[[k]] <- zero_mat(n, 1)
        next
      }

      m_i_safe <- m_i
      m_i_safe[!valid] <- Inf
      Dinv_w <- w_hat / m_i_safe
      Dinv_w[!is.finite(Dinv_w)] <- 0
      norm_Dinv_w <- sqrt(sum(Dinv_w^2))

      if (!is.finite(norm_Dinv_w) || norm_Dinv_w <= tiny) {
        vAkJ[[k]] <- zero_mat(n, 1)
        next
      }

      w_circ <- Dinv_w / norm_Dinv_w
      vAkJ[[k]] <- Vk %*% matrix(w_circ, ncol = 1)

      m_hat[k] <- 1 / norm_Dinv_w
      V_hat[k] <- sum((w_hat^2) * (w_circ^2) * V_i)
      a_tilde[k] <- m_hat[k]^2
      V_tilde[k] <- V_hat[k]
    }

    Vsum <- Reduce(`+`, vAkJ)
    T2 <- sum(Vsum^2)

    mK <- pmin(pmax(m_hat, 0), 1)
    pair_sum <- 0
    if (K >= 2) {
      for (i in seq_len(K - 1)) {
        for (j in (i + 1):K) {
          pair_sum <- pair_sum + mK[i] * mK[j]
        }
      }
      mu0 <- K + 2 * pair_sum
    } else {
      mu0 <- K
    }

    cosphi <- cos(alt_angle_deg * pi / 180)
    mu_alt <- K + 2 * cosphi * pair_sum

    V0 <- 0
    for (i in seq_len(K)) {
      others <- setdiff(seq_len(K), i)
      s_others <- sum(mK[others])
      V0 <- V0 + 4 * V_hat[i] * (s_others^2)
    }
    if (K >= 2) {
      for (i in seq_len(K - 1)) {
        for (j in (i + 1):K) {
          V0 <- V0 + 4 * (1 - mK[i]^2) * (1 - mK[j]^2)
        }
      }
    }
    V0 <- max(V0, 0)

    Zn <- NA_real_
    crit <- NA_real_
    pval <- NA_real_
    alpha_s_used <- NA_real_

    if (test_mode == "alpha") {
      if (V0 > tiny) {
        Zn <- sqrt(n / V0) * (T2 - mu0)
        alpha_s_used <- alpha
        crit <- qnorm(alpha_s_used)
        pval <- pnorm(Zn)
      }
      accept_joint <- (!is.na(Zn)) && (Zn >= crit)
    } else if (test_mode == "alpha_seq") {
      if (V0 > tiny) {
        Zn <- sqrt(n / V0) * (T2 - mu0)
        alpha_s_used <- alpha / (2^s)
        crit <- qnorm(alpha_s_used)
        pval <- pnorm(Zn)
      }
      accept_joint <- (!is.na(Zn)) && (Zn >= crit)
    } else {
      crit <- mu_alt
      if (V0 > tiny) Zn <- sqrt(n / V0) * (T2 - mu0)
      pval <- NA_real_
      accept_joint <- T2 > crit
    }

    T2_path[s] <- T2
    crit_path[s] <- crit
    mu0_path[s] <- mu0
    muAlt_path[s] <- mu_alt
    V0_path[s] <- V0
    Z_path[s] <- Zn
    accept_path[s] <- isTRUE(accept_joint)

    entry <- list(
      s = s,
      T2 = T2,
      mu = mu0,
      mu_alt = mu_alt,
      V = V0,
      Z = Zn,
      p_value = pval,
      crit = crit,
      alpha_s = alpha_s_used,
      m_hat = m_hat,
      a_tilde = a_tilde,
      V_tilde = V_tilde,
      accept = isTRUE(accept_joint)
    )
    test_history[[length(test_history) + 1L]] <- entry

    if (verbose) {
      cat(sprintf("Scan s=%d: T2=%.3f, Z=%.3f, crit=%.3f, mode=%s/%s, accept=%s\n",
                  s, T2, Zn, crit, test_mode, mode, as.character(isTRUE(accept_joint))))
    }

    if (prefix_alive && isTRUE(accept_joint)) {
      VJ_list[[length(VJ_list) + 1L]] <- vJ_hat
      accepted_s <- c(accepted_s, s)
      rank_accept_path[s] <- TRUE
    } else if (prefix_alive && !isTRUE(accept_joint)) {
      first_reject_s <- s
      prefix_alive <- FALSE
    }
  }

  rJ_hat <- length(VJ_list)

  if (rJ_hat > 0) {
    VJ_raw <- do.call(cbind, VJ_list)
    svJ <- svd(VJ_raw, nu = rJ_hat, nv = 0)
    V_joint <- svJ$u[, seq_len(rJ_hat), drop = FALSE]
    P_joint <- V_joint %*% t(V_joint)
  } else {
    V_joint <- zero_mat(n, 0)
    P_joint <- zero_mat(n, n)
  }

  rec <- vector("list", K)
  V_block_joint <- vector("list", K)
  V_block_indiv <- vector("list", K)
  U_block_joint <- vector("list", K)
  U_block_indiv <- vector("list", K)

  for (k in seq_len(K)) {
    Vk0 <- VA_list[[k]]
    pk <- nrow(X_list[[k]])

    PVAk <- if (ncol(Vk0) > 0) Vk0 %*% t(Vk0) else zero_mat(n, n)

    if (rJ_hat > 0 && ncol(Vk0) > 0) {
      Ak <- PVAk %*% V_joint
      svA <- svd(Ak, nu = min(dim(Ak)), nv = 0)
      r_eff <- sum(svA$d > tol_proj)
      if (r_eff > 0) {
        Vj_block <- svA$u[, seq_len(r_eff), drop = FALSE]
        Pcond <- Vj_block %*% t(Vj_block)
      } else {
        Vj_block <- zero_mat(n, 0)
        Pcond <- zero_mat(n, n)
      }
    } else {
      Vj_block <- zero_mat(n, 0)
      Pcond <- zero_mat(n, n)
    }

    Pind <- PVAk - Pcond

    if (any(Pind != 0)) {
      svI <- svd(Pind, nu = n, nv = 0)
      r_ind <- sum(svI$d > tol_proj)
      if (r_ind > 0) {
        Vind <- svI$u[, seq_len(r_ind), drop = FALSE]
      } else {
        Vind <- zero_mat(n, 0)
      }
    } else {
      Vind <- zero_mat(n, 0)
    }

    V_block_joint[[k]] <- Vj_block
    V_block_indiv[[k]] <- Vind

    if (ncol(Vj_block) > 0) {
      Xj_k <- X_list[[k]] %*% Vj_block
      svUj <- svd(Xj_k, nu = ncol(Vj_block), nv = 0)
      U_block_joint[[k]] <- svUj$u[, seq_len(ncol(Vj_block)), drop = FALSE]
    } else {
      U_block_joint[[k]] <- zero_mat(pk, 0)
    }

    if (ncol(Vind) > 0) {
      Xi_k <- X_list[[k]] %*% Vind
      svUi <- svd(Xi_k, nu = ncol(Vind), nv = 0)
      U_block_indiv[[k]] <- svUi$u[, seq_len(ncol(Vind)), drop = FALSE]
    } else {
      U_block_indiv[[k]] <- zero_mat(pk, 0)
    }

    rec[[k]] <- list(
      P_VAk = PVAk,
      P_VAk_given_VJ = Pcond,
      P_indiv = Pind,
      V_joint_block = Vj_block,
      V_indiv = Vind,
      U_joint_block = U_block_joint[[k]],
      U_indiv = U_block_indiv[[k]]
    )
  }

  out <- pack_output(rJ_hat, V_joint, P_joint, accepted_s, test_history, rec)
  out$V_block_joint <- V_block_joint
  out$V_block_indiv <- V_block_indiv
  out$U_block_joint <- U_block_joint
  out$U_block_indiv <- U_block_indiv

  X_joint_hat <- vector("list", K)
  X_indiv_hat <- vector("list", K)
  X_signal_hat <- vector("list", K)

  for (k in seq_len(K)) {
    Xk <- X_list[[k]]
    X_joint_hat[[k]] <- Xk %*% P_joint
    Pind <- out$blocks[[k]]$P_indiv
    X_indiv_hat[[k]] <- Xk %*% Pind
    X_signal_hat[[k]] <- X_joint_hat[[k]] + X_indiv_hat[[k]]
  }

  out$X_joint_hat <- X_joint_hat
  out$X_indiv_hat <- X_indiv_hat
  out$X_signal_hat <- X_signal_hat

  out$test_path <- list(
    T2 = as.numeric(T2_path),
    crit = as.numeric(crit_path),
    mu0 = as.numeric(mu0_path),
    mu_alt = as.numeric(muAlt_path),
    V0 = as.numeric(V0_path),
    Z = as.numeric(Z_path),
    accept = as.logical(accept_path),
    rank_accept = as.logical(rank_accept_path)
  )
  out$T2_path <- as.numeric(T2_path)
  out$crit_path <- as.numeric(crit_path)
  out$accept_path <- as.logical(accept_path)
  out$first_reject_s <- first_reject_s
  out$rank_accept_path <- as.logical(rank_accept_path)

  out
}

BEMA <- function(data, alpha = 0.2){
  n <- dim(data)[1]; p <- dim(data)[2]
  S <- t(data) %*% data  /n
  l <- eigen(S)$values

  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  predictor=qmp(k/min(p,n),max(n,p),min(n,p))*max(p,n)/n
  
  sigma2hat=lm(rev(l[k])~predictor-1)$coef[[1]]
  
  return(sigma2hat)
}
