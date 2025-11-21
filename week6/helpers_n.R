# helpers_n.R
# Utility helpers (no plotting). Numerical-stability optimized version.
# ---------------------------------------------------------------

# ===============================================================
# 0) Small utilities
# ===============================================================

# Clamp any numeric vector to (eps, 1 - eps) to avoid qnorm(0/1) = ±Inf
.clamp01 <- function(u, eps = 1e-12) pmin(1 - eps, pmax(eps, u))

# Positive value cap / winsorization:
#  - removes non-finite or negative values
#  - winsorizes at the (qcap)-quantile
#  - hard-caps at hard_cap
.cap_pos <- function(x, hard_cap = 1e3, qcap = 0.999) {
  x[!is.finite(x)] <- NA_real_
  x[x < 0] <- NA_real_
  qc <- suppressWarnings(stats::quantile(x, probs = qcap, na.rm = TRUE, names = FALSE))
  if (is.finite(qc)) x <- pmin(x, qc)
  pmin(x, hard_cap)
}

# ===============================================================
# 1) Density evaluators
# ===============================================================

# Evaluate fitted Bernstein copula density on an [0,1]^2 grid
eval_cop_density_grid <- function(fit, N = 150) {
  u <- v <- seq(0, 1, length.out = N)
  grid <- expand.grid(U = u, V = v)
  dens <- kdecopula::dkdecop(as.matrix(grid), fit)
  list(u = u, v = v, df = transform(grid, density = dens))
}

# Bivariate normal PDF with correlation rho
dnorm2 <- function(x, y, rho) {
  z <- (x^2 - 2*rho*x*y + y^2) / (1 - rho^2)
  (1 / (2*pi*sqrt(1 - rho^2))) * exp(-0.5 * z)
}

# Numerically stable Gaussian copula density c(u,v; rho)
# Uses log-form derivation to avoid overflow when |rho| ~ 1 or u,v ~ 0,1
c_gauss_copula <- function(u, v, rho) {
  u <- .clamp01(u); v <- .clamp01(v)
  x <- qnorm(u);    y <- qnorm(v)
  one_minus_r2 <- pmax(1 - rho^2, 1e-14)
  # log c(u,v; rho)
  logc <- 0.5*(x^2 + y^2) - 0.5*((x^2 - 2*rho*x*y + y^2) / one_minus_r2) - 0.5*log(one_minus_r2)
  out <- exp(pmin(logc, 700))      # exp(>709) would overflow
  out[!is.finite(out)] <- NA_real_
  out
}

# ===============================================================
# 2) Grid-based evaluation and RMSE (integrated over [0,1]^2)
# ===============================================================

# Evaluate ĉ and c_true on a K×K midpoint grid (robust version)
grid_eval_copula <- function(fit, rho, K = 12, cap_hard = 1e3, cap_q = 0.999) {
  mids <- (seq(0, 1, length.out = K + 1)[-1] + seq(0, 1, length.out = K + 1)[-(K+1)]) / 2
  g <- expand.grid(U = mids, V = mids)
  
  UVsafe <- cbind(.clamp01(g$U), .clamp01(g$V))
  c_hat  <- kdecopula::dkdecop(UVsafe, fit)
  c_true <- c_gauss_copula(g$U, g$V, rho)
  
  # Winsorize / cap extreme densities to prevent blow-ups at small m
  c_hat  <- .cap_pos(c_hat, hard_cap = cap_hard, qcap = cap_q)
  c_true <- .cap_pos(c_true, hard_cap = cap_hard, qcap = 1.0)
  
  bad <- !is.finite(c_hat) | !is.finite(c_true)
  if (any(bad)) {
    c_hat[bad]  <- NA_real_
    c_true[bad] <- NA_real_
  }
  
  g$c_hat  <- c_hat
  g$c_true <- c_true
  g$diff2  <- (c_hat - c_true)^2
  attr(g, "K") <- K
  g
}

# Compute grid-integrated ISE and RMSE (ignore NAs safely)
grid_rmse <- function(g) {
  K  <- attr(g, "K"); du <- 1 / K; dv <- 1 / K
  n_ok <- sum(is.finite(g$diff2))
  ise  <- sum(g$diff2, na.rm = TRUE) * du * dv
  rmse <- sqrt(ise)
  out <- list(ise = ise, rmse = rmse)
  attr(out, "n_used") <- n_ok
  out
}

# ===============================================================
# 3) Cell-probability comparison (Chi-square style)
# ===============================================================

cell_probs_from_samples <- function(UV, K = 12) {
  br <- seq(0, 1, length.out = K + 1)
  ci <- cut(UV[,1], br, include.lowest = TRUE, right = TRUE)
  cj <- cut(UV[,2], br, include.lowest = TRUE, right = TRUE)
  tab <- table(ci, cj)
  as.matrix(tab) / nrow(UV)
}

cell_probs_from_fit <- function(fit, K = 12) {
  mids <- (seq(0, 1, length.out = K + 1)[-1] + seq(0, 1, length.out = K + 1)[-(K+1)]) / 2
  du <- 1 / K; dv <- 1 / K
  M <- matrix(0, K, K)
  for (i in 1:K) for (j in 1:K) {
    M[i, j] <- kdecopula::dkdecop(cbind(mids[i], mids[j]), fit) * du * dv
  }
  M / sum(M)
}

chisq_cells <- function(p_hat, p_emp, eps = 1e-9) {
  sum((p_hat - p_emp)^2 / (p_emp + eps))
}

# ===============================================================
# 4) RMSE sweep over m (grid integration)
# ===============================================================

rmse_for_m <- function(UV, rho, m, K = 12, cap_hard = 1e3, cap_q = 0.999) {
  m <- max(as.integer(m), 4L)
  
  fit <- kdecopula::kdecop(UV, method = "bern", knots = m)
  
  
  g   <- grid_eval_copula(fit, rho = rho, K = K, cap_hard = cap_hard, cap_q = cap_q)
  out <- grid_rmse(g)
  c(m = m, rmse = out$rmse, ise = out$ise)
}

rmse_sweep_m <- function(UV, rho, m_vec, K = 12, cap_hard = 1e3, cap_q = 0.999) {
  ans <- lapply(m_vec, function(m) rmse_for_m(UV, rho, m, K, cap_hard, cap_q))
  do.call(rbind, ans) |> as.data.frame()
}

# ===============================================================
# 5) Point-sampled RMSE (professor’s definition) and replicates
# ===============================================================

# Generate evaluation points in [0,1]^2
make_eval_points <- function(N_eval = NULL, K_eval = NULL, scheme = c("midpoint", "random")) {
  scheme <- match.arg(scheme)
  if (!is.null(K_eval)) {
    mids <- (seq(0, 1, length.out = K_eval + 1)[-1] +
               seq(0, 1, length.out = K_eval + 1)[-(K_eval+1)]) / 2
    return(expand.grid(U = mids, V = mids))
  }
  if (is.null(N_eval)) N_eval <- 2000L
  if (scheme == "random") {
    return(data.frame(U = runif(N_eval), V = runif(N_eval)))
  } else {
    K <- ceiling(sqrt(N_eval))
    mids <- (seq(0, 1, length.out = K + 1)[-1] +
               seq(0, 1, length.out = K + 1)[-(K+1)]) / 2
    expand.grid(U = mids, V = mids)
  }
}

# RMSE at given m using fixed evaluation points (professor’s formula)
rmse_points_for_m <- function(UV_train, rho, m, UV_eval, cap_hard = 1e3, cap_q = 0.999) {
  m <- max(as.integer(m), 4L)
  fit <- kdecopula::kdecop(UV_train, method = "bern", knots = m)
  UVs <- cbind(.clamp01(UV_eval$U), .clamp01(UV_eval$V))
  c_hat  <- kdecopula::dkdecop(UVs, fit)
  c_true <- c_gauss_copula(UV_eval$U, UV_eval$V, rho)
  c_hat  <- .cap_pos(c_hat, hard_cap = cap_hard, qcap = cap_q)
  c_true <- .cap_pos(c_true, hard_cap = cap_hard, qcap = 1.0)
  dif2   <- (c_hat - c_true)^2
  dif2   <- dif2[is.finite(dif2)]
  rmse   <- sqrt(mean(dif2))
  c(m = m, rmse = rmse)
}

rmse_points_sweep <- function(UV_train, rho, m_vec, UV_eval, cap_hard = 1e3, cap_q = 0.999) {
  ans <- lapply(m_vec, function(m) rmse_points_for_m(UV_train, rho, m, UV_eval, cap_hard, cap_q))
  do.call(rbind, ans) |> as.data.frame()
}

# One replicate: simulate data and compute point-RMSE curve
one_replicate_rmse_curve <- function(n, rho, m_vec, UV_eval, dimd = 2L, cap_hard = 1e3, cap_q = 0.999) {
  nc  <- copula::normalCopula(param = rho, dim = dimd, dispstr = "ex")
  U   <- copula::rCopula(n, nc)
  UV  <- U[, 1:2, drop = FALSE]
  colnames(UV) <- c("U","V")
  rmse_df <- rmse_points_sweep(UV, rho, m_vec, UV_eval, cap_hard, cap_q)
  list(rmse_df = rmse_df, UV = UV)
}

# Repeat R times: return mean ± SE RMSE by m, m* distribution, raw matrix
replicate_rmse_experiment <- function(R, n, rho, m_vec, UV_eval,
                                      dimd = 2L, seed = 1L,
                                      cap_hard = 1e3, cap_q = 0.999) {
  set.seed(seed)
  mats <- matrix(NA_real_, nrow = R, ncol = length(m_vec))
  m_star <- integer(R)
  for (r in 1:R) {
    out  <- one_replicate_rmse_curve(n, rho, m_vec, UV_eval, dimd = dimd, cap_hard = cap_hard, cap_q = cap_q)
    vals <- out$rmse_df$rmse
    mats[r, ] <- vals
    m_star[r] <- m_vec[ which.min(vals) ]
  }
  mean_rmse <- colMeans(mats)
  se_rmse   <- apply(mats, 2, sd) / sqrt(R)
  data_curve <- data.frame(m = m_vec, mean_rmse = mean_rmse, se = se_rmse)
  list(curve = data_curve, m_star = m_star, rmse_mat = mats)
}

# Compare two m values (e.g., m0 = m*, m1 = neighbor) across replicates
# Outputs t-statistics and optional bootstrap CI for mean difference
compare_m_within_reps <- function(rmse_mat, m_vec, m0, m1,
                                  B = 0L, conf = 0.95, seed = 1L) {
  idx0 <- which(m_vec == m0); idx1 <- which(m_vec == m1)
  stopifnot(length(idx0) == 1L, length(idx1) == 1L)
  d <- rmse_mat[, idx1] - rmse_mat[, idx0]  # positive → m1 worse than m0
  R <- length(d)
  md <- mean(d); sd_d <- sd(d); se_d <- sd_d / sqrt(R)
  tval <- md / se_d
  pval <- 2 * pt(-abs(tval), df = R - 1)
  out <- list(mean_diff = md, se = se_d, t = tval, df = R - 1, p = pval)
  
  if (B > 0) {
    set.seed(seed)
    boot_means <- replicate(B, mean(sample(d, replace = TRUE)))
    a <- (1 - conf)/2; b <- 1 - a
    ci <- quantile(boot_means, probs = c(a, b), names = FALSE)
    out$boot_mean_ci <- ci
    out$boot_conf <- conf
  }
  out
}

