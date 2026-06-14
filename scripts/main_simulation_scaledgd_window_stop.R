# Conditional large-N selected-entry simulation: multi-method comparison for UNIF, ORACLE-PRO, pilot-PRO, UNIF-DB-SP one-step, and r-sample nuclear estimator.
#
# This optimized version caches uniform-sample estimators across targets: UNIF-MAT,
# UNIF-DB-SP, and NUCLEAR-R are computed once for each (draw, r) uniform sample and
# then reused for every target W_S. ORACLE-PRO and PILOT-PRO remain target-specific.
#
# Purpose:
#   Generate one full data set, compute UNIF-MAT and ORACLE-PRO probabilities once,
#   save probabilities, seeds, and sampled ids, and then repeat subsampling/estimation
#   B times on the same full data set.
#
# Compared methods:
#   UNIF-MAT       uniform subsampling probability, pi_i = 1/N, two-step WLS.
#   ORACLE-PRO     infeasible defensive probability using true epsilon and Phi.
#   PILOT-PRO      feasible defensive probability using pilot residuals and pilot Phi.
#   UNIF-DB-SP    uniform debiased spectral projection on one uniform r-sample,
#                  using the pilot estimator as an independent initial estimator, without
#                  splitting the r-sample.
#   NUCLEAR-R      nuclear-norm penalized least squares using only one uniform r-sample.
#
# MSE outputs:
#   1. true-centered empirical MSE:     || tau_hat - tau0 ||_2^2
#   2. full-data-centered empirical MSE:|| tau_hat - tau_full ||_2^2
#   3. linear oracle diagnostic MSE using tau_lin = tau_full + r^{-1} sum{zeta/(N*pi)-bar_zeta}.
#
# Defaults requested by the user for this sensitivity test:
#   N = 1e5, p = 10, q = 10, rank = 2, init_size = 200,
#   r_grid = 100,200,300,400,500,600,700,800,
#   n_rep = 500 conditional subsampling repetitions.
#   This variant scales both design covariance factors by a user-defined coefficient a:
#     Sigma_p_new = a * Sigma_p_raw, Sigma_q_new = a * Sigma_q_raw.
#   The coefficient is controlled by --cov-scale=a.
#   Default output is under 0608result/.
#   This version supports three initialization choices for two-step fitting:
#     --init-method=nuclear    nuclear-norm penalized least-squares initialization;
#     --init-method=scaled_gd  factorized ScaledGD initialization on the pilot sample;
#     --init-method=minimax    artificial oracle-like minimax initialization.
#   All three initializers are computed for diagnostics and comparison.
#
# This simplified version fixes rank=2 and singular values signal=(6,4.4), and uses
# one prespecified finite-dimensional target only. By default the target selects
# three entries of M0: M[2,2], M[2,3], and M[3,2]. The selected entries can be
# changed by --selected-entries=2:2,2:3,3:2. Uniform-sample-based estimators
# UNIF-MAT, UNIF-DB-SP, and NUCLEAR-R are cached once per (draw,r) and reused
# for this target. ORACLE-PRO and PILOT-PRO still use target-specific probabilities.
#
# Example:
#   Rscript simulation2_rank2_signal64_single_target3_multi_methods_unifdbsp_uniformcache_conditional_parallel_linear.R \
#     --p=15 --q=15 --cov-scale=0.3 --init-size=1000 \
#     --r-grid=1000,1200,1500,1800,2000 --n-rep=500 --n-cores=3

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[length(hit)])
}

get_bool_arg <- function(name, default = TRUE) {
  value <- tolower(get_arg(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
    stop("--", name, " must be true or false.")
  }
  value %in% c("true", "1", "yes")
}

parse_numeric_grid <- function(value) {
  out <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  if (any(!is.finite(out)) || length(out) == 0) stop("Invalid numeric grid: ", value)
  out
}


signal_label <- function(signal) {
  paste0("signal_", paste(gsub("\\.", "p", format(signal, trim = TRUE, scientific = FALSE)), collapse = "_"))
}

# Row/column target grid utilities.
# This experiment constructs selected-entry targets along a single row and/or
# a single column. For example, with --row-index=2 and --col-index=2,
# row targets are subsets of {(2,1),...,(2,q)} and column targets are subsets
# of {(1,2),...,(p,2)}.

parse_integer_grid <- function(value) {
  out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (any(!is.finite(out)) || length(out) == 0) stop("Invalid integer grid: ", value)
  out
}

parse_selected_entries <- function(value) {
  pieces <- strsplit(value, ",", fixed = TRUE)[[1]]
  pieces <- trimws(pieces)
  if (length(pieces) == 0 || any(nchar(pieces) == 0)) {
    stop("--selected-entries must be like 2:2,2:3,3:2")
  }
  mat <- do.call(rbind, lapply(pieces, function(piece) {
    xy <- strsplit(piece, ":", fixed = TRUE)[[1]]
    if (length(xy) != 2L) stop("Invalid entry in --selected-entries: ", piece)
    as.integer(xy)
  }))
  if (ncol(mat) != 2L || any(!is.finite(mat))) {
    stop("--selected-entries must be like 2:2,2:3,3:2")
  }
  colnames(mat) <- c("row", "col")
  attr(mat, "target_kind") <- "custom"
  attr(mat, "target_anchor") <- NA_integer_
  attr(mat, "target_size") <- nrow(mat)
  mat
}

make_axis_order <- function(m, order = "sequential", seed = 20260604L) {
  idx <- seq_len(m)
  if (order == "sequential") {
    idx
  } else if (order == "spread") {
    center <- (m + 1) / 2
    idx[order(abs(idx - center), idx)]
  } else if (order == "random") {
    set.seed(seed)
    sample(idx, length(idx))
  } else {
    stop("--target-order must be one of: sequential, spread, random.")
  }
}

make_rowcol_target_grid <- function(p, q, sizes, row_index = 2L, col_index = 2L,
                                    mode = "both", order = "sequential", seed = 20260604L) {
  if (!mode %in% c("both", "row", "column")) {
    stop("--target-mode must be one of: both, row, column.")
  }
  if (row_index < 1 || row_index > p) stop("--row-index must be between 1 and p.")
  if (col_index < 1 || col_index > q) stop("--col-index must be between 1 and q.")
  if (any(sizes < 1) || any(sizes > max(p, q))) {
    stop("Every target size must be between 1 and max(p,q).")
  }

  targets <- list()
  col_order <- make_axis_order(q, order = order, seed = seed)
  row_order <- make_axis_order(p, order = order, seed = seed + 1000L)

  if (mode %in% c("both", "row")) {
    for (d in sizes) {
      if (d > q) stop("Row target size ", d, " exceeds q=", q, ".")
      mat <- cbind(row = rep(row_index, d), col = col_order[seq_len(d)])
      attr(mat, "target_kind") <- "row"
      attr(mat, "target_anchor") <- row_index
      attr(mat, "target_size") <- d
      targets[[paste0("row_", row_index, "_d", d)]] <- mat
    }
  }
  if (mode %in% c("both", "column")) {
    for (d in sizes) {
      if (d > p) stop("Column target size ", d, " exceeds p=", p, ".")
      mat <- cbind(row = row_order[seq_len(d)], col = rep(col_index, d))
      attr(mat, "target_kind") <- "column"
      attr(mat, "target_anchor") <- col_index
      attr(mat, "target_size") <- d
      targets[[paste0("col_", col_index, "_d", d)]] <- mat
    }
  }
  targets
}

target_label <- function(selected_pos) {
  kind <- attr(selected_pos, "target_kind")
  anchor <- attr(selected_pos, "target_anchor")
  if (is.null(kind) || is.null(anchor)) {
    if (length(unique(selected_pos[, 1])) == 1L && nrow(selected_pos) > 1L) {
      kind <- "row"; anchor <- unique(selected_pos[, 1])
    } else if (length(unique(selected_pos[, 2])) == 1L && nrow(selected_pos) > 1L) {
      kind <- "column"; anchor <- unique(selected_pos[, 2])
    } else {
      kind <- "custom"; anchor <- NA_integer_
    }
  }
  checksum <- sum(selected_pos[, 1] * 1000L + selected_pos[, 2])
  anchor_text <- if (is.na(anchor)) "xx" else sprintf("%02d", as.integer(anchor))
  paste0(kind, anchor_text, "_", sprintf("%03d", nrow(selected_pos)), "entries_ck", checksum)
}

N_arg <- get_arg("N", get_arg("n", "100000"))
N_default <- as.integer(N_arg)

if (!requireNamespace("parallel", quietly = TRUE)) {
  stop("The base R package 'parallel' is required.")
}
available_cores <- max(1L, parallel::detectCores(logical = TRUE))
default_cores <- max(1L, available_cores - 1L)

cfg <- list(
  seed = as.integer(get_arg("seed", "20260527")),
  n = N_default,
  p = as.integer(get_arg("p", "10")),
  q = as.integer(get_arg("q", "10")),
  rank = as.integer(get_arg("rank", "2")),
  sigma = as.numeric(get_arg("sigma", "0.5")),
  signal = parse_numeric_grid(get_arg("signal", "6.0,4.4")),
  init_size = as.integer(get_arg("init-size", "200")),
  init_method = get_arg("init-method", "scaled_gd"),
  nuclear_lambda = get_arg("nuclear-lambda", "auto"),
  nuclear_lambda_factor = as.numeric(get_arg("nuclear-lambda-factor", "1.0")),
  nuclear_maxit = as.integer(get_arg("nuclear-maxit", "300")),
  nuclear_tol = as.numeric(get_arg("nuclear-tol", "1e-6")),
  scaledgd_eta = as.numeric(get_arg("scaledgd-eta", "0.5")),
  scaledgd_maxit = as.integer(get_arg("scaledgd-maxit", "500")),
  scaledgd_tol = as.numeric(get_arg("scaledgd-tol", "1e-6")),
  scaledgd_frob_tol = as.numeric(get_arg("scaledgd-frob-tol", "1e-3")),
  scaledgd_stop_rule = get_arg("scaledgd-stop-rule", "objective"),
  scaledgd_stop_window = as.integer(get_arg("scaledgd-stop-window", "1")),
  scaledgd_ridge = as.numeric(get_arg("scaledgd-ridge", "1e-6")),
  scaledgd_backtrack = get_bool_arg("scaledgd-backtrack", TRUE),
  r_nuclear_lambda = get_arg("r-nuclear-lambda", get_arg("nuclear-r-lambda", "auto")),
  r_nuclear_lambda_factor = as.numeric(get_arg("r-nuclear-lambda-factor", get_arg("nuclear-r-lambda-factor", "1.0"))),
  run_unif_db_sp = local({
    v_new <- get_arg("run-unif-db-sp", NULL)
    if (is.null(v_new)) {
      get_bool_arg("run-xia-yuan", TRUE)
    } else {
      v <- tolower(v_new)
      if (!v %in% c("true", "false", "1", "0", "yes", "no")) {
        stop("--run-unif-db-sp must be true or false.")
      }
      v %in% c("true", "1", "yes")
    }
  }),
  run_r_nuclear = get_bool_arg("run-r-nuclear", TRUE),
  r_grid = as.integer(parse_numeric_grid(get_arg("r-grid", "100,200,300,400,500,600,700,800"))),
  n_rep = as.integer(get_arg("n-rep", "500")),
  n_cores = as.integer(get_arg("n-cores", as.character(default_cores))),
  defensive_gamma = as.numeric(get_arg("defensive-gamma", "0.10")),
  rho_p = as.numeric(get_arg("rho-p", "0.8")),
  rho_q = as.numeric(get_arg("rho-q", "0.7")),
  design = get_arg("design", "kronecker"),
  add_full_component = get_bool_arg("add-full-component", TRUE),
  save_sampled_ids = get_bool_arg("save-sampled-ids", TRUE),
  save_probabilities = get_bool_arg("save-probabilities", TRUE),
  out_root = normalizePath(get_arg("out-root", file.path("0608result", "results_rank2_signal_6_4p4_single_target3_covscale_custom_scaledgd_multi_methods_unifdbsp_uniformcache_conditional")), winslash = "/", mustWork = FALSE),
  cov_scale = as.numeric(get_arg("cov-scale", "1.0")),
  target_grid = NULL,
  selected_entries = get_arg("selected-entries", "2:2,2:3,3:2"),
  target_size_grid = parse_integer_grid(get_arg("target-size-grid", "1,2,3,5,8,10")),  # unused in this single-target version; kept for compatibility
  target_order = get_arg("target-order", "sequential"),
  target_seed = as.integer(get_arg("target-seed", "20260604")),
  target_mode = get_arg("target-mode", "both"),
  row_index = as.integer(get_arg("row-index", "2")),
  col_index = as.integer(get_arg("col-index", "2")),
  selected_pos = NULL,
  target_name = NULL
)

# Normalize aliases for initialization method.
cfg$init_method <- tolower(gsub("-", "_", cfg$init_method))
if (cfg$init_method == "scaledgd") cfg$init_method <- "scaled_gd"

# Build one fixed selected-entry target.
# Default: W selects vec(M0) entries M[2,2], M[2,3], and M[3,2].
cfg$selected_pos <- parse_selected_entries(cfg$selected_entries)
if (any(cfg$selected_pos[, 1] < 1 | cfg$selected_pos[, 1] > cfg$p |
        cfg$selected_pos[, 2] < 1 | cfg$selected_pos[, 2] > cfg$q)) {
  stop("--selected-entries contains entries outside the p x q matrix.")
}
cfg$target_grid <- list(cfg$selected_pos)
names(cfg$target_grid) <- target_label(cfg$selected_pos)

if (!cfg$design %in% c("isotropic", "kronecker")) {
  stop("--design must be isotropic or kronecker.")
}
if (cfg$p < 1 || cfg$q < 1) stop("p and q must be positive integers.")
if (cfg$rank < 1 || cfg$rank > min(cfg$p, cfg$q)) stop("Invalid rank.")
if (length(cfg$signal) != cfg$rank) stop("--signal must have length equal to rank. For this test, use --rank=2 --signal=6.0,4.4.")
if (cfg$n < 10) stop("N must be at least 10.")
if (cfg$init_size <= 0) stop("--init-size must be positive.")
if (cfg$init_size >= cfg$n) stop("--init-size must be smaller than N so that a nonempty main pool remains.")
if (any(cfg$r_grid < 2)) stop("All r values must be at least 2.")
if (any(cfg$r_grid < cfg$init_size)) {
  stop("Require r >= init_size. Please set --r-grid with all values >= --init-size.")
}
if (any(cfg$r_grid > cfg$n)) stop("All r values must be <= N.")
if (cfg$n_rep < 1) stop("--n-rep must be positive.")
if (!is.finite(cfg$n_cores) || cfg$n_cores < 1) stop("--n-cores must be a positive integer.")
cfg$n_cores <- max(1L, min(as.integer(cfg$n_cores), available_cores, cfg$n_rep))
if (cfg$defensive_gamma < 0 || cfg$defensive_gamma >= 1) {
  stop("--defensive-gamma must be in [0,1).")
}
if (!is.finite(cfg$cov_scale) || cfg$cov_scale <= 0) {
  stop("--cov-scale must be a positive finite number.")
}
if (!cfg$init_method %in% c("nuclear", "minimax", "scaled_gd")) {
  stop("--init-method must be one of: nuclear, scaled_gd, minimax.")
}
if (!is.finite(cfg$scaledgd_eta) || cfg$scaledgd_eta <= 0) {
  stop("--scaledgd-eta must be a positive finite number.")
}
if (!is.finite(cfg$scaledgd_maxit) || cfg$scaledgd_maxit < 1) {
  stop("--scaledgd-maxit must be a positive integer.")
}
if (!is.finite(cfg$scaledgd_tol) || cfg$scaledgd_tol <= 0) {
  stop("--scaledgd-tol must be a positive finite number.")
}
if (!is.finite(cfg$scaledgd_frob_tol) || cfg$scaledgd_frob_tol <= 0) {
  stop("--scaledgd-frob-tol must be a positive finite number.")
}
if (!cfg$scaledgd_stop_rule %in% c("objective", "frob", "both", "either", "window")) {
  stop("--scaledgd-stop-rule must be one of: objective, frob, both, either, window.")
}
if (!is.finite(cfg$scaledgd_stop_window) || cfg$scaledgd_stop_window < 1L) {
  stop("--scaledgd-stop-window must be a positive integer.")
}
if (!is.finite(cfg$scaledgd_ridge) || cfg$scaledgd_ridge < 0) {
  stop("--scaledgd-ridge must be a nonnegative finite number.")
}
if (!is.finite(cfg$nuclear_lambda_factor) || cfg$nuclear_lambda_factor <= 0) {
  stop("--nuclear-lambda-factor must be a positive finite number.")
}
if (!is.finite(cfg$nuclear_maxit) || cfg$nuclear_maxit < 1) {
  stop("--nuclear-maxit must be a positive integer.")
}
if (!is.finite(cfg$nuclear_tol) || cfg$nuclear_tol <= 0) {
  stop("--nuclear-tol must be a positive finite number.")
}
if (!is.finite(cfg$r_nuclear_lambda_factor) || cfg$r_nuclear_lambda_factor <= 0) {
  stop("--r-nuclear-lambda-factor must be a positive finite number.")
}

mat_fnorm <- function(A) sqrt(sum(A * A))
row_l2 <- function(A) sqrt(rowSums(A * A))

ar1_cov <- function(d, rho) {
  idx <- seq_len(d)
  outer(idx, idx, function(i, j) rho^abs(i - j))
}

exchangeable_cov <- function(d, rho) {
  (1 - rho) * diag(d) + rho * matrix(1, d, d)
}

rand_orth <- function(n, k) {
  qr.Q(qr(matrix(rnorm(n * k), n, k)))[, seq_len(k), drop = FALSE]
}

make_truth <- function(p, q, rank0, signal) {
  U <- rand_orth(p, rank0)
  V <- rand_orth(q, rank0)
  if (length(signal) < rank0) {
    signal <- seq(signal[1], signal[1] * 0.7, length.out = rank0)
  }
  D <- diag(signal[seq_len(rank0)], nrow = rank0)
  list(M0 = U %*% D %*% t(V), U0 = U, V0 = V, D0 = D)
}

make_xmat <- function(n, p, q, Sigmap, Sigmaq) {
  z <- matrix(rnorm(n * p * q), n, p * q)
  K <- kronecker(chol(Sigmaq), chol(Sigmap))
  z %*% K
}

transpose_index <- function(p, q) {
  idx <- integer(p * q)
  for (j in seq_len(q)) {
    for (i in seq_len(p)) {
      old <- i + (j - 1L) * p
      new <- j + (i - 1L) * q
      idx[new] <- old
    }
  }
  idx
}

entry_index <- function(row, col, p) row + (col - 1L) * p

safe_solve <- function(A, ridge = 1e-9) {
  A <- as.matrix(A)
  scale <- mean(abs(diag(A)))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  solve(A + ridge * scale * diag(nrow(A)))
}

symmetrize <- function(A) (A + t(A)) / 2

phi_from_factors <- function(L, R, Sigmap, Sigmaq) {
  PL <- L %*% safe_solve(t(L) %*% Sigmap %*% L) %*% t(L)
  PR <- R %*% safe_solve(t(R) %*% Sigmaq %*% R) %*% t(R)
  Phi <- kronecker(PR, safe_solve(Sigmap)) +
    kronecker(safe_solve(Sigmaq), PL) -
    kronecker(PR, PL)
  list(Phi = symmetrize(Phi), L = L, R = R)
}

phi_from_matrix <- function(M, rank0, Sigmap, Sigmaq) {
  S <- svd(M, nu = rank0, nv = rank0)
  U <- S$u[, seq_len(rank0), drop = FALSE]
  V <- S$v[, seq_len(rank0), drop = FALSE]
  L <- U %*% diag(S$d[seq_len(rank0)], nrow = rank0)
  R <- V
  phi_from_factors(L, R, Sigmap, Sigmaq)
}

artificial_minimax_initial <- function(M0, sigma, rank0, init_size) {
  p <- nrow(M0)
  q <- ncol(M0)
  target_error <- sigma * sqrt(rank0 * (p + q) / init_size)
  perturbation <- matrix(rnorm(p * q), p, q)
  perturbation <- perturbation / mat_fnorm(perturbation)
  M0 + target_error * perturbation
}

svd_soft_threshold <- function(M, tau) {
  S <- svd(M)
  d <- pmax(S$d - tau, 0)
  keep <- which(d > 0)
  if (length(keep) == 0) return(matrix(0, nrow(M), ncol(M)))
  S$u[, keep, drop = FALSE] %*% diag(d[keep], nrow = length(keep)) %*% t(S$v[, keep, drop = FALSE])
}

nuclear_norm_value <- function(M) {
  sum(svd(M, nu = 0, nv = 0)$d)
}

nuclear_ls_objective <- function(M, xmat, y, lambda) {
  resid <- as.vector(xmat %*% as.vector(M) - y)
  0.5 * mean(resid^2) + lambda * nuclear_norm_value(M)
}

estimate_lipschitz_ls <- function(xmat) {
  G <- crossprod(xmat) / nrow(xmat)
  L <- max(eigen((G + t(G)) / 2, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(L) || L <= 0) L <- 1
  L
}

nuclear_norm_initial <- function(xmat, y, p, q, sigma, lambda = "auto",
                                 lambda_factor = 1.0, maxit = 300, tol = 1e-6) {
  m <- nrow(xmat)
  if (lambda == "auto") {
    lambda_value <- lambda_factor * sigma * sqrt((p + q) / m)
  } else {
    lambda_value <- as.numeric(lambda)
    if (!is.finite(lambda_value) || lambda_value < 0) {
      stop("--nuclear-lambda must be 'auto' or a nonnegative number.")
    }
  }

  L <- estimate_lipschitz_ls(xmat)
  M_old <- matrix(0, p, q)
  Y <- M_old
  t_old <- 1
  obj_old <- nuclear_ls_objective(M_old, xmat, y, lambda_value)
  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(maxit)) {
    resid <- as.vector(xmat %*% as.vector(Y) - y)
    grad <- matrix(crossprod(xmat, resid) / m, p, q)
    M_new <- svd_soft_threshold(Y - grad / L, lambda_value / L)
    obj_new <- nuclear_ls_objective(M_new, xmat, y, lambda_value)
    rel_change <- abs(obj_old - obj_new) / (1 + abs(obj_old))

    if (is.finite(rel_change) && rel_change < tol) {
      converged <- TRUE
      M_old <- M_new
      obj_old <- obj_new
      break
    }

    t_new <- (1 + sqrt(1 + 4 * t_old^2)) / 2
    Y <- M_new + ((t_old - 1) / t_new) * (M_new - M_old)
    M_old <- M_new
    t_old <- t_new
    obj_old <- obj_new
  }

  svals <- svd(M_old, nu = 0, nv = 0)$d
  list(
    M = M_old,
    lambda = lambda_value,
    lipschitz = L,
    iter = iter,
    converged = converged,
    objective = obj_old,
    rank_1e8 = sum(svals > 1e-8),
    rank_1e4 = sum(svals > 1e-4),
    singular_values = svals
  )
}


factor_objective <- function(xmat, y, Lfac, Rfac) {
  M <- Lfac %*% t(Rfac)
  resid <- as.vector(xmat %*% as.vector(M) - y)
  0.5 * mean(resid^2)
}

scaledgd_precond <- function(G, ridge = 1e-6) {
  G <- symmetrize(as.matrix(G))
  scale <- mean(abs(diag(G)))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  safe_solve(G + ridge * scale * diag(nrow(G)), ridge = 0)
}

moment_initial_matrix <- function(xmat, y, p, q, S1_raw = NULL, S2_raw = NULL, avefro = NULL) {
  m <- nrow(xmat)
  if (!is.null(S1_raw) && !is.null(S2_raw) && !is.null(avefro)) {
    invS1 <- safe_solve(S1_raw)
    invS2scaled <- safe_solve(S2_raw / max(avefro, .Machine$double.eps))
    M0 <- matrix(0, p, q)
    for (ii in seq_len(m)) {
      Xi <- matrix(xmat[ii, ], p, q)
      M0 <- M0 + invS1 %*% Xi %*% invS2scaled * y[ii]
    }
    M0 / m
  } else {
    matrix(crossprod(xmat, y) / m, p, q)
  }
}

scaledgd_should_stop <- function(stop_rule, rel_obj, rel_frob, obj_tol, frob_tol) {
  obj_ok <- is.finite(rel_obj) && rel_obj < obj_tol
  frob_ok <- is.finite(rel_frob) && rel_frob < frob_tol
  if (stop_rule == "objective") {
    obj_ok
  } else if (stop_rule == "frob") {
    frob_ok
  } else if (stop_rule == "both") {
    obj_ok && frob_ok
  } else if (stop_rule == "either") {
    obj_ok || frob_ok
  } else {
    stop("--scaledgd-stop-rule must be one of: objective, frob, both, either, window.")
  }
}

scaledgd_window_stop <- function(obj_history, rel_frob_history, iter, window, obj_tol, frob_tol) {
  if (iter < window) {
    return(list(stop = FALSE, rel_objective = NA_real_, rel_frob = NA_real_))
  }
  obj_start <- obj_history[iter - window + 1L]
  obj_end <- obj_history[iter + 1L]
  rel_obj_window <- (obj_start - obj_end) / ((1 + abs(obj_start)) * window)
  idx <- (iter - window + 1L):iter
  rel_frob_window <- mean(rel_frob_history[idx], na.rm = TRUE)
  stop <- is.finite(rel_obj_window) && is.finite(rel_frob_window) &&
    rel_obj_window < obj_tol && rel_frob_window < frob_tol
  list(stop = stop, rel_objective = rel_obj_window, rel_frob = rel_frob_window)
}

scaledgd_initial <- function(xmat, y, p, q, rank0, sigma = NULL,
                             S1_raw = NULL, S2_raw = NULL, avefro = NULL,
                             eta = 0.5, maxit = 500, tol = 1e-6,
                             frob_tol = 1e-3, stop_rule = "objective",
                             stop_window = 1L,
                             ridge = 1e-6, backtrack = TRUE) {
  m <- nrow(xmat)
  M_start <- moment_initial_matrix(xmat, y, p, q, S1_raw = S1_raw, S2_raw = S2_raw, avefro = avefro)
  init <- rank_truncate(M_start, rank0)
  d0 <- pmax(init$d[seq_len(rank0)], 0)
  sqrtD <- diag(sqrt(d0), nrow = rank0)
  Lfac <- init$U %*% sqrtD
  Rfac <- init$V %*% sqrtD

  obj_old <- factor_objective(xmat, y, Lfac, Rfac)
  obj_start <- obj_old
  converged <- FALSE
  iter <- 0L
  eta_last <- eta
  backtrack_count_total <- 0L
  rel_change <- NA_real_
  rel_frob <- NA_real_
  window_rel_change <- NA_real_
  window_rel_frob <- NA_real_
  stop_streak <- 0L
  obj_history <- rep(NA_real_, maxit + 1L)
  rel_frob_history <- rep(NA_real_, maxit)
  obj_history[1L] <- obj_old

  for (iter in seq_len(maxit)) {
    M_current <- Lfac %*% t(Rfac)
    resid <- as.vector(xmat %*% as.vector(M_current) - y)
    gradM <- matrix(crossprod(xmat, resid) / m, p, q)
    gradL <- gradM %*% Rfac
    gradR <- t(gradM) %*% Lfac
    preR <- scaledgd_precond(t(Rfac) %*% Rfac, ridge = ridge)
    preL <- scaledgd_precond(t(Lfac) %*% Lfac, ridge = ridge)

    step <- eta
    accepted <- FALSE
    obj_new <- NA_real_
    Lnew <- Lfac
    Rnew <- Rfac
    n_backtrack <- 0L

    while (!accepted) {
      Lcand <- Lfac - step * gradL %*% preR
      Rcand <- Rfac - step * gradR %*% preL
      obj_cand <- factor_objective(xmat, y, Lcand, Rcand)
      if (!backtrack || (is.finite(obj_cand) && obj_cand <= obj_old + 1e-12)) {
        accepted <- TRUE
        Lnew <- Lcand
        Rnew <- Rcand
        obj_new <- obj_cand
        eta_last <- step
      } else if (n_backtrack >= 20L) {
        accepted <- TRUE
        Lnew <- Lfac
        Rnew <- Rfac
        obj_new <- obj_old
        eta_last <- 0
      } else {
        step <- step / 2
        n_backtrack <- n_backtrack + 1L
      }
    }
    backtrack_count_total <- backtrack_count_total + n_backtrack

    rel_change <- abs(obj_old - obj_new) / (1 + abs(obj_old))
    M_new <- Lnew %*% t(Rnew)
    rel_frob <- mat_fnorm(M_new - M_current) / (1 + mat_fnorm(M_current))
    obj_history[iter + 1L] <- obj_new
    rel_frob_history[iter] <- rel_frob
    Lfac <- Lnew
    Rfac <- Rnew

    if (stop_rule == "window") {
      window_check <- scaledgd_window_stop(
        obj_history, rel_frob_history, iter, stop_window, tol, frob_tol
      )
      stop_condition <- window_check$stop
      window_rel_change <- window_check$rel_objective
      window_rel_frob <- window_check$rel_frob
    } else {
      stop_condition <- scaledgd_should_stop(stop_rule, rel_change, rel_frob, tol, frob_tol)
      window_rel_change <- NA_real_
      window_rel_frob <- NA_real_
    }
    if (stop_condition) {
      stop_streak <- stop_streak + 1L
    } else {
      stop_streak <- 0L
    }
    if ((stop_rule == "window" && stop_condition) ||
        (stop_rule != "window" && stop_streak >= stop_window)) {
      converged <- TRUE
      obj_old <- obj_new
      break
    }
    obj_old <- obj_new
  }

  M_final <- Lfac %*% t(Rfac)
  M_rank <- rank_truncate(M_final, rank0)
  svals <- svd(M_rank$M, nu = 0, nv = 0)$d
  list(
    M = M_rank$M,
    M_untruncated = M_final,
    L = Lfac,
    R = Rfac,
    U = M_rank$U,
    V = M_rank$V,
    init_objective = obj_start,
    objective = obj_old,
    rel_objective = rel_change,
    rel_frob = rel_frob,
    window_rel_objective = window_rel_change,
    window_rel_frob = window_rel_frob,
    iter = iter,
    converged = converged,
    eta = eta,
    eta_last = eta_last,
    backtrack_total = backtrack_count_total,
    ridge = ridge,
    frob_tol = frob_tol,
    stop_rule = stop_rule,
    stop_window = stop_window,
    stop_streak = stop_streak,
    rank_1e8 = sum(svals > 1e-8),
    rank_1e4 = sum(svals > 1e-4),
    singular_values = svals
  )
}

orth_basis <- function(A) {
  qr.Q(qr(A))[, seq_len(ncol(A)), drop = FALSE]
}

subspace_distance <- function(Rhat, R0) {
  Qhat <- orth_basis(Rhat)
  Q0 <- orth_basis(R0)
  Pdiff <- Qhat %*% t(Qhat) - Q0 %*% t(Q0)
  eig <- eigen((Pdiff + t(Pdiff)) / 2, symmetric = TRUE, only.values = TRUE)$values
  list(frob = mat_fnorm(Pdiff), op = max(abs(eig)))
}

evaluate_initial_estimator <- function(M_init, M0, R_init, R0, target, name,
                                       init_used = FALSE, extra = list()) {
  err <- M_init - M0
  tau_init <- as.vector(target$W_S %*% as.vector(M_init))
  tau0 <- as.vector(target$W_S %*% as.vector(M0))
  tau_err <- tau_init - tau0
  svals <- svd(M_init, nu = 0, nv = 0)$d
  subdist <- subspace_distance(R_init, R0)

  out <- data.frame(
    init_estimator = name,
    init_used_for_twostep = init_used,
    target_label = target$target_label,
    target_size = target$target_size,
    matrix_bias_mean = mean(err),
    matrix_abs_bias_mean = mean(abs(err)),
    matrix_mse_per_entry = mean(err^2),
    matrix_frob_sq = sum(err^2),
    matrix_frob_error = mat_fnorm(err),
    matrix_relative_frob_error = mat_fnorm(err) / max(mat_fnorm(M0), .Machine$double.eps),
    target_bias_mean = mean(tau_err),
    target_abs_bias_mean = mean(abs(tau_err)),
    target_mse_total = sum(tau_err^2),
    target_mse_per_entry = mean(tau_err^2),
    target_rmse_total = sqrt(sum(tau_err^2)),
    rank_1e8 = sum(svals > 1e-8),
    rank_1e4 = sum(svals > 1e-4),
    top_singular_values = paste(round(svals[seq_len(min(5, length(svals)))], 8), collapse = ","),
    R_subspace_frob = subdist$frob,
    R_subspace_op = subdist$op
  )

  if (length(extra) > 0) {
    for (nm in names(extra)) out[[nm]] <- extra[[nm]]
  }
  out
}


regularize_cov <- function(A, ridge = 1e-6) {
  A <- symmetrize(as.matrix(A))
  scale <- mean(abs(diag(A)))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  A + ridge * scale * diag(nrow(A))
}

pilot_covariance_estimates <- function(xmat, p, q, ridge = 1e-6) {
  m <- nrow(xmat)
  S1_raw <- matrix(0, p, p)
  S2_raw <- matrix(0, q, q)
  avefro <- 0
  for (ii in seq_len(m)) {
    Xi <- matrix(xmat[ii, ], p, q)
    S1_raw <- S1_raw + Xi %*% t(Xi)
    S2_raw <- S2_raw + t(Xi) %*% Xi
    avefro <- avefro + sum(Xi * Xi)
  }
  S1_raw <- S1_raw / m
  S2_raw <- S2_raw / m
  avefro <- avefro / m
  scale <- sqrt(max(avefro, .Machine$double.eps))
  Sigmap_hat <- regularize_cov(S1_raw / scale, ridge = ridge)
  Sigmaq_hat <- regularize_cov(S2_raw / scale, ridge = ridge)
  list(
    S1_raw = regularize_cov(S1_raw, ridge = ridge),
    S2_raw = regularize_cov(S2_raw, ridge = ridge),
    avefro = avefro,
    Sigmap_hat = Sigmap_hat,
    Sigmaq_hat = Sigmaq_hat
  )
}

rank_truncate <- function(M, rank0) {
  S <- svd(M, nu = rank0, nv = rank0)
  U <- S$u[, seq_len(rank0), drop = FALSE]
  V <- S$v[, seq_len(rank0), drop = FALSE]
  D <- diag(S$d[seq_len(rank0)], nrow = rank0)
  list(M = U %*% D %*% t(V), U = U, V = V, D = D, d = S$d)
}

unif_db_sp_one_step <- function(xmat, y, M_init, U_init, V_init,
                              S1_raw, S2_raw, avefro, p, q, rank0) {
  resid <- as.vector(y - xmat %*% as.vector(M_init))
  invS1 <- safe_solve(S1_raw)
  invS2scaled <- safe_solve(S2_raw / max(avefro, .Machine$double.eps))
  Delta <- matrix(0, p, q)
  for (ii in seq_len(nrow(xmat))) {
    Xi <- matrix(xmat[ii, ], p, q)
    Delta <- Delta + invS1 %*% Xi %*% invS2scaled * resid[ii]
  }
  M_unbs <- M_init + Delta / nrow(xmat)
  U_new <- svd(M_unbs %*% V_init, nu = rank0, nv = 0)$u[, seq_len(rank0), drop = FALSE]
  V_new <- svd(t(M_unbs) %*% U_init, nu = rank0, nv = 0)$u[, seq_len(rank0), drop = FALSE]
  M_proj <- U_new %*% t(U_new) %*% M_unbs %*% V_new %*% t(V_new)
  list(M = M_proj, M_unbs = M_unbs, U = U_new, V = V_new)
}

weighted_lm <- function(A, y, ipw, ridge = 1e-9) {
  w <- as.numeric(ipw)
  w <- w / mean(w)
  sw <- sqrt(w)
  Aw <- A * sw
  yw <- y * sw
  solve_mat <- crossprod(Aw)
  scale <- mean(diag(solve_mat))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  as.vector(solve(solve_mat + ridge * scale * diag(ncol(A)), crossprod(Aw, yw)))
}

matrix_two_step <- function(xmat, xmatT, y, ipw, R_anchor, p, q, rank0) {
  A_L <- kronecker(t(R_anchor), diag(p))
  Z_L <- xmat %*% t(A_L)
  beta_L <- weighted_lm(Z_L, y, ipw)
  Lhat <- matrix(beta_L, p, rank0)
  A_R <- kronecker(t(Lhat), diag(q))
  Z_R <- xmatT %*% t(A_R)
  beta_R <- weighted_lm(Z_R, y, ipw)
  Rhat <- matrix(beta_R, q, rank0)
  list(M = Lhat %*% t(Rhat), L = Lhat, R = Rhat)
}

make_pi <- function(score, gamma) {
  n <- length(score)
  score[!is.finite(score)] <- 0
  score <- pmax(score, 0)
  raw <- if (sum(score) > 0) score / sum(score) else rep(1 / n, n)
  pi <- (1 - gamma) * raw + gamma / n
  pi / sum(pi)
}

sampling_objective <- function(a, pi) {
  n <- length(a)
  if (length(pi) != n) stop("sampling_objective requires a and pi to have the same length.")
  den <- n * pi
  terms <- rep(NA_real_, n)
  zero_num <- a == 0
  positive_den <- den > 0 & is.finite(den)
  terms[zero_num & !positive_den] <- 0
  terms[positive_den] <- a[positive_den]^2 / den[positive_den]
  terms[!zero_num & !positive_den] <- Inf
  mean(terms)
}

conditional_gamma_trace <- function(zeta, pi) {
  # zeta is N x d with rows zeta_i = W_S Phi_true x_i epsilon_i.
  # Conditional covariance trace of one draw variable:
  #   Var_*( zeta_I/(N*pi_I) ) has trace
  #   N^{-2} sum_i ||zeta_i||^2/pi_i - ||bar_zeta||^2.
  n <- nrow(zeta)
  if (length(pi) != n) stop("conditional_gamma_trace requires nrow(zeta) == length(pi).")
  zeta_norm2 <- rowSums(zeta * zeta)
  uncentered <- mean(zeta_norm2 / (n * pi))
  zeta_bar <- colMeans(zeta)
  centered <- uncentered - sum(zeta_bar * zeta_bar)
  max(centered, 0)
}

effective_sample_size <- function(pi) {
  if (length(pi) == 0) return(NA_real_)
  1 / sum(pi * pi)
}

probability_summary <- function(pi, method) {
  data.frame(
    method = method,
    pi_min = min(pi),
    pi_q01 = as.numeric(quantile(pi, 0.01, names = FALSE)),
    pi_median = median(pi),
    pi_mean = mean(pi),
    pi_q99 = as.numeric(quantile(pi, 0.99, names = FALSE)),
    pi_max = max(pi),
    n_pi_min = length(pi) * min(pi),
    n_pi_max = length(pi) * max(pi),
    effective_sample_size = effective_sample_size(pi)
  )
}

cv <- function(x) {
  m <- mean(x)
  if (!is.finite(m) || abs(m) < .Machine$double.eps) return(NA_real_)
  sd(x) / m
}

make_selected_target <- function(p, q, selected_pos = NULL) {
  if (is.null(selected_pos)) {
    selected_pos <- matrix(c(2, 2, 2, 3, 3, 2, 3, 3, 4, 4), ncol = 2, byrow = TRUE)
  }
  selected_pos <- as.matrix(selected_pos)
  if (ncol(selected_pos) != 2) stop("selected_pos must have two columns: row and col.")
  if (any(selected_pos[, 1] < 1 | selected_pos[, 1] > p |
          selected_pos[, 2] < 1 | selected_pos[, 2] > q)) {
    stop("selected_pos contains entries outside the p x q matrix.")
  }
  theta_idx <- entry_index(selected_pos[, 1], selected_pos[, 2], p)
  W_S <- matrix(0, nrow = length(theta_idx), ncol = p * q)
  W_S[cbind(seq_along(theta_idx), theta_idx)] <- 1
  list(
    selected_pos = selected_pos,
    theta_idx = theta_idx,
    W_S = W_S,
    target_size = nrow(selected_pos),
    target_label = target_label(selected_pos),
    entry_names = paste0("M[", selected_pos[, 1], ",", selected_pos[, 2], "]")
  )
}

make_design_objects <- function(cfg) {
  Sigmap_raw <- if (cfg$design == "isotropic") diag(cfg$p) else ar1_cov(cfg$p, cfg$rho_p)
  Sigmaq_raw <- if (cfg$design == "isotropic") diag(cfg$q) else exchangeable_cov(cfg$q, cfg$rho_q)
  Sigmap <- cfg$cov_scale * Sigmap_raw
  Sigmaq <- cfg$cov_scale * Sigmaq_raw
  list(Sigmap = Sigmap, Sigmaq = Sigmaq, Sigmap_raw = Sigmap_raw, Sigmaq_raw = Sigmaq_raw)
}

seed_for_data <- function(cfg) {
  offset <- if (cfg$design == "kronecker") 100000L else 0L
  cfg$seed + offset + 1000000L
}

seed_for_init <- function(cfg) {
  offset <- if (cfg$design == "kronecker") 100000L else 0L
  cfg$seed + offset + 2000000L
}

method_names <- function(cfg = NULL) {
  out <- c("UNIF-MAT", "ORACLE-PRO", "PILOT-PRO")
  if (is.null(cfg) || isTRUE(cfg$run_unif_db_sp)) out <- c(out, "UNIF-DB-SP")
  if (is.null(cfg) || isTRUE(cfg$run_r_nuclear)) out <- c(out, "NUCLEAR-R")
  out
}

seed_for_sampling <- function(draw_id, r_sub, method, cfg) {
  method_offset <- switch(method,
    "UNIF-MAT" = 10L,
    "ORACLE-PRO" = 20L,
    "PILOT-PRO" = 30L,
    "UNIF-DB-SP" = 40L,
    "NUCLEAR-R" = 50L,
    90L
  )
  cfg$seed + 3000000L + 10000L * draw_id + 10L * r_sub + method_offset
}

generate_full_data <- function(seed_data, truth, cfg, Sigmap, Sigmaq) {
  set.seed(seed_data)
  beta0 <- as.vector(truth$M0)
  xmat <- make_xmat(cfg$n, cfg$p, cfg$q, Sigmap, Sigmaq)
  eps <- cfg$sigma * rnorm(cfg$n)
  y <- as.vector(xmat %*% beta0 + eps)
  list(xmat = xmat, eps = eps, y = y)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_se <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  n_eff <- sum(!is.na(x))
  if (n_eff <= 1) return(NA_real_)
  sd(x, na.rm = TRUE) / sqrt(n_eff)
}

safe_rmse <- function(x) {
  m <- safe_mean(x)
  if (is.na(m)) return(NA_real_)
  sqrt(m)
}

summarise_mse <- function(raw_results, raw_theory) {
  groups <- split(raw_results, list(raw_results$r, raw_results$method), drop = TRUE)
  pieces <- lapply(groups, function(z) {
    th <- raw_theory[raw_theory$r == z$r[1] & raw_theory$method == z$method[1], ]
    data.frame(
      r = z$r[1],
      method = z$method[1],
      empirical_MSE_tau0_mean = safe_mean(z$squared_error_tau0),
      empirical_MSE_tau0_sd = safe_sd(z$squared_error_tau0),
      empirical_MSE_tau0_se = safe_se(z$squared_error_tau0),
      empirical_RMSE_tau0 = safe_rmse(z$squared_error_tau0),
      empirical_MSE_taufull_mean = safe_mean(z$squared_error_taufull),
      empirical_MSE_taufull_sd = safe_sd(z$squared_error_taufull),
      empirical_MSE_taufull_se = safe_se(z$squared_error_taufull),
      empirical_RMSE_taufull = safe_rmse(z$squared_error_taufull),
      empirical_linear_MSE_tau0_mean = safe_mean(z$squared_error_linear_tau0),
      empirical_linear_MSE_tau0_sd = safe_sd(z$squared_error_linear_tau0),
      empirical_linear_MSE_tau0_se = safe_se(z$squared_error_linear_tau0),
      empirical_linear_RMSE_tau0 = safe_rmse(z$squared_error_linear_tau0),
      empirical_linear_MSE_taufull_mean = safe_mean(z$squared_error_linear_taufull),
      empirical_linear_MSE_taufull_sd = safe_sd(z$squared_error_linear_taufull),
      empirical_linear_MSE_taufull_se = safe_se(z$squared_error_linear_taufull),
      empirical_linear_RMSE_taufull = safe_rmse(z$squared_error_linear_taufull),
      theory_conditional_MSE = unique(th$theory_conditional_mse),
      theory_true_centered_MSE = unique(th$theory_true_centered_mse),
      theory_linear_conditional_MSE = unique(th$theory_linear_conditional_mse),
      theory_linear_true_centered_MSE = unique(th$theory_linear_true_centered_mse),
      average_unique_sampled = safe_mean(z$sampled_unique_n),
      average_duplicate_rate = safe_mean(z$sampled_duplicate_rate)
    )
  })
  out <- do.call(rbind, pieces)
  out$method <- as.character(out$method)
  method_order <- match(out$method, method_names(NULL))
  out[order(out$r, method_order), ]
}

summarise_ratios <- function(summary_df) {
  rs <- sort(unique(summary_df$r))
  pieces <- lapply(rs, function(r_value) {
    z <- summary_df[summary_df$r == r_value, ]
    u <- z[z$method == "UNIF-MAT", ]
    o <- z[z$method == "ORACLE-PRO", ]
    data.frame(
      r = r_value,
      empirical_tau0_ratio_oracle_over_unif = o$empirical_MSE_tau0_mean / u$empirical_MSE_tau0_mean,
      empirical_tau0_gain_oracle_vs_unif = 1 - o$empirical_MSE_tau0_mean / u$empirical_MSE_tau0_mean,
      empirical_taufull_ratio_oracle_over_unif = o$empirical_MSE_taufull_mean / u$empirical_MSE_taufull_mean,
      empirical_taufull_gain_oracle_vs_unif = 1 - o$empirical_MSE_taufull_mean / u$empirical_MSE_taufull_mean,
      empirical_linear_tau0_ratio_oracle_over_unif = o$empirical_linear_MSE_tau0_mean / u$empirical_linear_MSE_tau0_mean,
      empirical_linear_tau0_gain_oracle_vs_unif = 1 - o$empirical_linear_MSE_tau0_mean / u$empirical_linear_MSE_tau0_mean,
      empirical_linear_taufull_ratio_oracle_over_unif = o$empirical_linear_MSE_taufull_mean / u$empirical_linear_MSE_taufull_mean,
      empirical_linear_taufull_gain_oracle_vs_unif = 1 - o$empirical_linear_MSE_taufull_mean / u$empirical_linear_MSE_taufull_mean,
      theory_conditional_ratio_oracle_over_unif = o$theory_conditional_MSE / u$theory_conditional_MSE,
      theory_conditional_gain_oracle_vs_unif = 1 - o$theory_conditional_MSE / u$theory_conditional_MSE,
      theory_true_centered_ratio_oracle_over_unif = o$theory_true_centered_MSE / u$theory_true_centered_MSE,
      theory_true_centered_gain_oracle_vs_unif = 1 - o$theory_true_centered_MSE / u$theory_true_centered_MSE,
      theory_linear_conditional_ratio_oracle_over_unif = o$theory_linear_conditional_MSE / u$theory_linear_conditional_MSE,
      theory_linear_conditional_gain_oracle_vs_unif = 1 - o$theory_linear_conditional_MSE / u$theory_linear_conditional_MSE,
      theory_linear_true_centered_ratio_oracle_over_unif = o$theory_linear_true_centered_MSE / u$theory_linear_true_centered_MSE,
      theory_linear_true_centered_gain_oracle_vs_unif = 1 - o$theory_linear_true_centered_MSE / u$theory_linear_true_centered_MSE
    )
  })
  do.call(rbind, pieces)
}


summarise_method_ratios <- function(summary_df, baseline = "UNIF-MAT") {
  rs <- sort(unique(summary_df$r))
  pieces <- list()
  idx <- 1L
  for (r_value in rs) {
    z <- summary_df[summary_df$r == r_value, ]
    base <- z[z$method == baseline, ]
    if (nrow(base) == 0) next
    for (method in unique(z$method)) {
      zz <- z[z$method == method, ]
      pieces[[idx]] <- data.frame(
        r = r_value,
        baseline = baseline,
        method = method,
        empirical_tau0_ratio_over_unif = zz$empirical_MSE_tau0_mean / base$empirical_MSE_tau0_mean,
        empirical_taufull_ratio_over_unif = zz$empirical_MSE_taufull_mean / base$empirical_MSE_taufull_mean,
        empirical_tau0_gain_vs_unif = 1 - zz$empirical_MSE_tau0_mean / base$empirical_MSE_tau0_mean,
        empirical_taufull_gain_vs_unif = 1 - zz$empirical_MSE_taufull_mean / base$empirical_MSE_taufull_mean
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, pieces)
}

plot_two_method_mse <- function(summary_df, file, y_col, title, ylab = "MSE") {
  pdf(file, width = 8, height = 5.5)
  on.exit(dev.off(), add = TRUE)
  methods <- intersect(method_names(NULL), unique(summary_df$method))
  methods <- methods[sapply(methods, function(mm) {
    any(is.finite(summary_df[summary_df$method == mm, y_col]))
  })]
  ylim <- range(summary_df[[y_col]], finite = TRUE)
  plot(NA, xlim = range(summary_df$r), ylim = ylim,
       xlab = "main subsample size r", ylab = ylab, main = title)
  grid()
  cols <- seq_along(methods)
  pchs <- 15 + seq_along(methods)
  for (jj in seq_along(methods)) {
    method <- methods[jj]
    z <- summary_df[summary_df$method == method, ]
    z <- z[order(z$r), ]
    lines(z$r, z[[y_col]], type = "b", lwd = 2, col = cols[jj], pch = pchs[jj])
  }
  legend("topright", legend = methods, col = cols, pch = pchs, lty = 1, bty = "n", cex = 0.75)
}

plot_empirical_vs_theory <- function(summary_df, file, empirical_col, theory_col, title) {
  pdf(file, width = 7.5, height = 5.2)
  on.exit(dev.off(), add = TRUE)
  cols <- c("UNIF-MAT" = "#2166AC", "ORACLE-PRO" = "#B2182B", "PILOT-PRO" = "#4D9221")
  pchs <- c("UNIF-MAT" = 16, "ORACLE-PRO" = 17, "PILOT-PRO" = 18)
  y_all <- c(summary_df[[empirical_col]], summary_df[[theory_col]])
  ylim <- range(y_all, finite = TRUE)
  plot(NA, xlim = range(summary_df$r), ylim = ylim,
       xlab = "main subsample size r", ylab = "MSE", main = title)
  grid()
  for (method in intersect(c("UNIF-MAT", "ORACLE-PRO", "PILOT-PRO"), unique(summary_df$method))) {
    z <- summary_df[summary_df$method == method, ]
    z <- z[order(z$r), ]
    lines(z$r, z[[empirical_col]], type = "b", lwd = 2, col = cols[method], pch = pchs[method])
    lines(z$r, z[[theory_col]], type = "b", lwd = 2, lty = 2, col = cols[method], pch = pchs[method])
  }
  legend("topright",
         legend = c("UNIF empirical", "ORACLE empirical", "PILOT empirical",
                    "UNIF theory", "ORACLE theory", "PILOT theory"),
         col = c(cols["UNIF-MAT"], cols["ORACLE-PRO"], cols["PILOT-PRO"],
                 cols["UNIF-MAT"], cols["ORACLE-PRO"], cols["PILOT-PRO"]),
         pch = c(16, 17, 18, 16, 17, 18), lty = c(1, 1, 1, 2, 2, 2), bty = "n", cex = 0.75)
}

plot_ratio <- function(ratio_df, file, empirical_col, theory_col, title) {
  pdf(file, width = 7, height = 5)
  on.exit(dev.off(), add = TRUE)
  y_all <- c(ratio_df[[empirical_col]], ratio_df[[theory_col]])
  plot(NA, xlim = range(ratio_df$r), ylim = range(y_all, 1, finite = TRUE),
       xlab = "main subsample size r", ylab = "ORACLE / UNIF MSE ratio", main = title)
  grid()
  lines(ratio_df$r, ratio_df[[empirical_col]], type = "b", lwd = 2, pch = 16)
  lines(ratio_df$r, ratio_df[[theory_col]], type = "b", lwd = 2, pch = 17, lty = 2)
  abline(h = 1, lty = 4)
  legend("topright", legend = c("Empirical", "Theory"), pch = c(16, 17), lty = c(1, 2), bty = "n")
}

prepare_shared_context <- function(cfg) {
  design <- make_design_objects(cfg)
  Sigmap <- design$Sigmap
  Sigmaq <- design$Sigmaq
  tr_idx <- transpose_index(cfg$p, cfg$q)

  set.seed(cfg$seed)
  truth <- make_truth(cfg$p, cfg$q, cfg$rank, cfg$signal)
  beta0 <- as.vector(truth$M0)

  seed_data <- seed_for_data(cfg)
  seed_init <- seed_for_init(cfg)

  message("Available cores: ", available_cores, "; using cores: ", cfg$n_cores)
  message("Generating one shared full data set...")
  dat <- generate_full_data(seed_data, truth, cfg, Sigmap, Sigmaq)
  xmatT <- dat$xmat[, tr_idx, drop = FALSE]

  message("Constructing shared initialization estimators...")
  set.seed(seed_init)
  M_minimax <- artificial_minimax_initial(truth$M0, cfg$sigma, cfg$rank, cfg$init_size)
  sv_minimax <- svd(M_minimax, nu = cfg$rank, nv = cfg$rank)
  R_minimax <- sv_minimax$v[, seq_len(cfg$rank), drop = FALSE]
  matrix_minimax_error <- mat_fnorm(M_minimax - truth$M0)

  init_n <- cfg$init_size
  set.seed(seed_init + 777L)
  init_idx <- sample.int(cfg$n, init_n, replace = FALSE)
  main_idx <- setdiff(seq_len(cfg$n), init_idx)
  n_main <- length(main_idx)
  if (n_main <= 0) stop("The main pool is empty after removing init_idx.")

  message("Computing shared pilot covariance estimates...")
  pilot_cov <- pilot_covariance_estimates(dat$xmat[init_idx, , drop = FALSE], cfg$p, cfg$q)

  nuclear_fit <- nuclear_norm_initial(
    dat$xmat[init_idx, , drop = FALSE], dat$y[init_idx],
    p = cfg$p, q = cfg$q, sigma = cfg$sigma,
    lambda = cfg$nuclear_lambda,
    lambda_factor = cfg$nuclear_lambda_factor,
    maxit = cfg$nuclear_maxit,
    tol = cfg$nuclear_tol
  )
  M_nuclear <- nuclear_fit$M
  sv_nuclear <- svd(M_nuclear, nu = cfg$rank, nv = cfg$rank)
  R_nuclear <- sv_nuclear$v[, seq_len(cfg$rank), drop = FALSE]
  matrix_nuclear_error <- mat_fnorm(M_nuclear - truth$M0)

  scaledgd_fit <- scaledgd_initial(
    dat$xmat[init_idx, , drop = FALSE], dat$y[init_idx],
    p = cfg$p, q = cfg$q, rank0 = cfg$rank, sigma = cfg$sigma,
    S1_raw = pilot_cov$S1_raw, S2_raw = pilot_cov$S2_raw, avefro = pilot_cov$avefro,
    eta = cfg$scaledgd_eta,
    maxit = cfg$scaledgd_maxit,
    tol = cfg$scaledgd_tol,
    frob_tol = cfg$scaledgd_frob_tol,
    stop_rule = cfg$scaledgd_stop_rule,
    stop_window = cfg$scaledgd_stop_window,
    ridge = cfg$scaledgd_ridge,
    backtrack = cfg$scaledgd_backtrack
  )
  M_scaledgd <- scaledgd_fit$M
  sv_scaledgd <- svd(M_scaledgd, nu = cfg$rank, nv = cfg$rank)
  R_scaledgd <- sv_scaledgd$v[, seq_len(cfg$rank), drop = FALSE]
  matrix_scaledgd_error <- mat_fnorm(M_scaledgd - truth$M0)

  if (cfg$init_method == "nuclear") {
    M_tilde <- M_nuclear
    U_tilde <- sv_nuclear$u[, seq_len(cfg$rank), drop = FALSE]
    R_tilde <- R_nuclear
  } else if (cfg$init_method == "scaled_gd") {
    M_tilde <- M_scaledgd
    U_tilde <- sv_scaledgd$u[, seq_len(cfg$rank), drop = FALSE]
    R_tilde <- R_scaledgd
  } else {
    M_tilde <- M_minimax
    U_tilde <- sv_minimax$u[, seq_len(cfg$rank), drop = FALSE]
    R_tilde <- R_minimax
  }
  init_method_used <- cfg$init_method
  matrix_init_error <- mat_fnorm(M_tilde - truth$M0)

  message("Computing shared true Phi...")
  true_phi <- phi_from_matrix(truth$M0, cfg$rank, Sigmap, Sigmaq)$Phi

  message("Computing shared main-pool full-data two-step estimator for tau_full...")
  full_fit <- matrix_two_step(
    dat$xmat[main_idx, , drop = FALSE], xmatT[main_idx, , drop = FALSE], dat$y[main_idx],
    ipw = rep(1, n_main), R_anchor = R_tilde,
    p = cfg$p, q = cfg$q, rank0 = cfg$rank
  )

  list(
    design = design,
    Sigmap = Sigmap,
    Sigmaq = Sigmaq,
    tr_idx = tr_idx,
    truth = truth,
    beta0 = beta0,
    seed_data = seed_data,
    seed_init = seed_init,
    dat = dat,
    xmatT = xmatT,
    init_idx = init_idx,
    main_idx = main_idx,
    n_main = n_main,
    M_minimax = M_minimax,
    R_minimax = R_minimax,
    sv_minimax = sv_minimax,
    matrix_minimax_error = matrix_minimax_error,
    nuclear_fit = nuclear_fit,
    M_nuclear = M_nuclear,
    R_nuclear = R_nuclear,
    sv_nuclear = sv_nuclear,
    matrix_nuclear_error = matrix_nuclear_error,
    scaledgd_fit = scaledgd_fit,
    M_scaledgd = M_scaledgd,
    R_scaledgd = R_scaledgd,
    sv_scaledgd = sv_scaledgd,
    matrix_scaledgd_error = matrix_scaledgd_error,
    M_tilde = M_tilde,
    U_tilde = U_tilde,
    R_tilde = R_tilde,
    init_method_used = init_method_used,
    matrix_init_error = matrix_init_error,
    pilot_cov = pilot_cov,
    true_phi = true_phi,
    full_fit = full_fit
  )
}


prepare_uniform_fit_cache <- function(cfg, shared) {
  uniform_methods <- intersect(c("UNIF-MAT", "UNIF-DB-SP", "NUCLEAR-R"), method_names(cfg))
  if (length(uniform_methods) == 0L) {
    return(list(fits = list(), index = data.frame()))
  }

  dat <- shared$dat
  xmatT <- shared$xmatT
  main_idx <- shared$main_idx
  n_main <- shared$n_main
  R_tilde <- shared$R_tilde
  M_tilde <- shared$M_tilde
  U_tilde <- shared$U_tilde
  pilot_cov <- shared$pilot_cov

  pi_unif <- rep(1 / n_main, n_main)

  fit_one_uniform_draw <- function(b) {
    local_fits <- list()
    local_index <- list()
    local_id <- 1L

    for (r_sub in cfg$r_grid) {
      uniform_seed <- seed_for_sampling(b, r_sub, "UNIF-MAT", cfg)
      set.seed(uniform_seed)
      id_local <- sample.int(n_main, r_sub, replace = TRUE, prob = pi_unif)
      id <- main_idx[id_local]
      xs <- dat$xmat[id, , drop = FALSE]
      ys <- dat$y[id]
      xTs <- xmatT[id, , drop = FALSE]

      sampled_unique_n <- length(unique(id))
      sampled_duplicate_rate <- 1 - sampled_unique_n / length(id)

      if ("UNIF-MAT" %in% uniform_methods) {
        fit <- matrix_two_step(
          xs, xTs, ys,
          ipw = 1 / pi_unif[id_local], R_anchor = R_tilde,
          p = cfg$p, q = cfg$q, rank0 = cfg$rank
        )
        key <- paste(b, r_sub, "UNIF-MAT", sep = "__")
        local_fits[[key]] <- list(M = fit$M, id = id, id_local = id_local, seed = uniform_seed,
                                  sampled_unique_n = sampled_unique_n,
                                  sampled_duplicate_rate = sampled_duplicate_rate)
        local_index[[local_id]] <- data.frame(draw = b, r = r_sub, method = "UNIF-MAT",
                                              seed_sampling = uniform_seed,
                                              sampled_unique_n = sampled_unique_n,
                                              sampled_duplicate_rate = sampled_duplicate_rate,
                                              nuclear_iter = NA_integer_,
                                              nuclear_converged = NA,
                                              nuclear_lambda = NA_real_)
        local_id <- local_id + 1L
      }

      if ("UNIF-DB-SP" %in% uniform_methods) {
        fit <- unif_db_sp_one_step(
          xs, ys, M_tilde, U_tilde, R_tilde,
          pilot_cov$S1_raw, pilot_cov$S2_raw, pilot_cov$avefro,
          p = cfg$p, q = cfg$q, rank0 = cfg$rank
        )
        key <- paste(b, r_sub, "UNIF-DB-SP", sep = "__")
        local_fits[[key]] <- list(M = fit$M, id = id, id_local = id_local, seed = uniform_seed,
                                  sampled_unique_n = sampled_unique_n,
                                  sampled_duplicate_rate = sampled_duplicate_rate)
        local_index[[local_id]] <- data.frame(draw = b, r = r_sub, method = "UNIF-DB-SP",
                                              seed_sampling = uniform_seed,
                                              sampled_unique_n = sampled_unique_n,
                                              sampled_duplicate_rate = sampled_duplicate_rate,
                                              nuclear_iter = NA_integer_,
                                              nuclear_converged = NA,
                                              nuclear_lambda = NA_real_)
        local_id <- local_id + 1L
      }

      if ("NUCLEAR-R" %in% uniform_methods) {
        fit <- nuclear_norm_initial(
          xs, ys, p = cfg$p, q = cfg$q, sigma = cfg$sigma,
          lambda = cfg$r_nuclear_lambda,
          lambda_factor = cfg$r_nuclear_lambda_factor,
          maxit = cfg$nuclear_maxit,
          tol = cfg$nuclear_tol
        )
        M_hat <- rank_truncate(fit$M, cfg$rank)$M
        key <- paste(b, r_sub, "NUCLEAR-R", sep = "__")
        local_fits[[key]] <- list(M = M_hat, id = id, id_local = id_local, seed = uniform_seed,
                                  sampled_unique_n = sampled_unique_n,
                                  sampled_duplicate_rate = sampled_duplicate_rate,
                                  nuclear_iter = fit$iter,
                                  nuclear_converged = fit$converged,
                                  nuclear_lambda = fit$lambda)
        local_index[[local_id]] <- data.frame(draw = b, r = r_sub, method = "NUCLEAR-R",
                                              seed_sampling = uniform_seed,
                                              sampled_unique_n = sampled_unique_n,
                                              sampled_duplicate_rate = sampled_duplicate_rate,
                                              nuclear_iter = fit$iter,
                                              nuclear_converged = fit$converged,
                                              nuclear_lambda = fit$lambda)
        local_id <- local_id + 1L
      }
    }

    list(
      fits = local_fits,
      index = if (length(local_index) > 0L) do.call(rbind, local_index) else data.frame()
    )
  }

  message("Precomputing target-invariant uniform-sample estimators across targets: ",
          paste(uniform_methods, collapse = ", "))
  draw_ids <- seq_len(cfg$n_rep)
  if (cfg$n_cores > 1L && .Platform$OS.type != "windows") {
    chunks <- parallel::mclapply(
      draw_ids, fit_one_uniform_draw,
      mc.cores = cfg$n_cores,
      mc.preschedule = FALSE
    )
  } else if (cfg$n_cores > 1L && .Platform$OS.type == "windows") {
    warning(
      "Uniform fit cache precomputation on Windows would copy the full N x pq matrix to each worker. ",
      "Falling back to sequential execution. Use Linux/macOS for fork-based parallelism."
    )
    chunks <- lapply(draw_ids, fit_one_uniform_draw)
  } else {
    chunks <- lapply(draw_ids, fit_one_uniform_draw)
  }

  fits <- unlist(lapply(chunks, function(z) z$fits), recursive = FALSE)
  index <- do.call(rbind, lapply(chunks, function(z) z$index))
  list(fits = fits, index = index)
}

run_all <- function(cfg, shared) {
  dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(cfg$out_root, "figures"), recursive = TRUE, showWarnings = FALSE)

  Sigmap <- shared$Sigmap
  Sigmaq <- shared$Sigmaq
  target <- make_selected_target(cfg$p, cfg$q, cfg$selected_pos)
  truth <- shared$truth
  beta0 <- shared$beta0
  tau0 <- as.vector(target$W_S %*% beta0)
  seed_data <- shared$seed_data
  seed_init <- shared$seed_init
  dat <- shared$dat
  xmatT <- shared$xmatT
  init_idx <- shared$init_idx
  main_idx <- shared$main_idx
  n_main <- shared$n_main
  x_main <- dat$xmat[main_idx, , drop = FALSE]
  xT_main <- xmatT[main_idx, , drop = FALSE]
  y_main <- dat$y[main_idx]
  eps_main <- dat$eps[main_idx]

  M_minimax <- shared$M_minimax
  R_minimax <- shared$R_minimax
  M_nuclear <- shared$M_nuclear
  R_nuclear <- shared$R_nuclear
  nuclear_fit <- shared$nuclear_fit
  M_scaledgd <- shared$M_scaledgd
  R_scaledgd <- shared$R_scaledgd
  scaledgd_fit <- shared$scaledgd_fit
  matrix_minimax_error <- shared$matrix_minimax_error
  matrix_nuclear_error <- shared$matrix_nuclear_error
  matrix_scaledgd_error <- shared$matrix_scaledgd_error
  M_tilde <- shared$M_tilde
  U_tilde <- shared$U_tilde
  R_tilde <- shared$R_tilde
  init_method_used <- shared$init_method_used
  matrix_init_error <- shared$matrix_init_error
  pilot_cov <- shared$pilot_cov
  true_phi <- shared$true_phi
  full_fit <- shared$full_fit
  tau_full <- as.vector(target$W_S %*% as.vector(full_fit$M))

  initial_diagnostics <- rbind(
    evaluate_initial_estimator(
      M_nuclear, truth$M0, R_nuclear, truth$V0, target,
      name = "nuclear_norm",
      init_used = (cfg$init_method == "nuclear"),
      extra = list(
        init_sample_size = cfg$init_size,
        lambda = nuclear_fit$lambda,
        lambda_factor = cfg$nuclear_lambda_factor,
        lipschitz = nuclear_fit$lipschitz,
        objective = nuclear_fit$objective,
        init_objective = NA_real_,
        iterations = nuclear_fit$iter,
        converged = nuclear_fit$converged,
        rel_objective = NA_real_,
        rel_frob = NA_real_,
        window_rel_objective = NA_real_,
        window_rel_frob = NA_real_,
        frob_tol = NA_real_,
        stop_rule = NA_character_,
        stop_window = NA_integer_,
        stop_streak = NA_integer_,
        eta = NA_real_,
        eta_last = NA_real_,
        ridge = NA_real_,
        backtrack_total = NA_integer_
      )
    ),
    evaluate_initial_estimator(
      M_scaledgd, truth$M0, R_scaledgd, truth$V0, target,
      name = "scaled_gd",
      init_used = (cfg$init_method == "scaled_gd"),
      extra = list(
        init_sample_size = cfg$init_size,
        lambda = NA_real_,
        lambda_factor = NA_real_,
        lipschitz = NA_real_,
        objective = scaledgd_fit$objective,
        init_objective = scaledgd_fit$init_objective,
        iterations = scaledgd_fit$iter,
        converged = scaledgd_fit$converged,
        rel_objective = scaledgd_fit$rel_objective,
        rel_frob = scaledgd_fit$rel_frob,
        window_rel_objective = scaledgd_fit$window_rel_objective,
        window_rel_frob = scaledgd_fit$window_rel_frob,
        frob_tol = scaledgd_fit$frob_tol,
        stop_rule = scaledgd_fit$stop_rule,
        stop_window = scaledgd_fit$stop_window,
        stop_streak = scaledgd_fit$stop_streak,
        eta = scaledgd_fit$eta,
        eta_last = scaledgd_fit$eta_last,
        ridge = scaledgd_fit$ridge,
        backtrack_total = scaledgd_fit$backtrack_total
      )
    ),
    evaluate_initial_estimator(
      M_minimax, truth$M0, R_minimax, truth$V0, target,
      name = "artificial_minimax",
      init_used = (cfg$init_method == "minimax"),
      extra = list(
        init_sample_size = cfg$init_size,
        lambda = NA_real_,
        lambda_factor = NA_real_,
        lipschitz = NA_real_,
        objective = NA_real_,
        init_objective = NA_real_,
        iterations = NA_integer_,
        converged = NA,
        rel_objective = NA_real_,
        rel_frob = NA_real_,
        window_rel_objective = NA_real_,
        window_rel_frob = NA_real_,
        frob_tol = NA_real_,
        stop_rule = NA_character_,
        stop_window = NA_integer_,
        stop_streak = NA_integer_,
        eta = NA_real_,
        eta_last = NA_real_,
        ridge = NA_real_,
        backtrack_total = NA_integer_
      )
    )
  )

  initial_comparison <- data.frame(
    target_label = target$target_label,
    target_size = target$target_size,
    init_method_used = init_method_used,
    nuclear_over_minimax_matrix_frob = matrix_nuclear_error / max(matrix_minimax_error, .Machine$double.eps),
    scaledgd_over_minimax_matrix_frob = matrix_scaledgd_error / max(matrix_minimax_error, .Machine$double.eps),
    scaledgd_over_nuclear_matrix_frob = matrix_scaledgd_error / max(matrix_nuclear_error, .Machine$double.eps),
    nuclear_over_minimax_matrix_mse = mean((M_nuclear - truth$M0)^2) / max(mean((M_minimax - truth$M0)^2), .Machine$double.eps),
    scaledgd_over_minimax_matrix_mse = mean((M_scaledgd - truth$M0)^2) / max(mean((M_minimax - truth$M0)^2), .Machine$double.eps),
    scaledgd_over_nuclear_matrix_mse = mean((M_scaledgd - truth$M0)^2) / max(mean((M_nuclear - truth$M0)^2), .Machine$double.eps),
    nuclear_over_minimax_target_mse =
      sum((as.vector(target$W_S %*% as.vector(M_nuclear)) - tau0)^2) /
      max(sum((as.vector(target$W_S %*% as.vector(M_minimax)) - tau0)^2), .Machine$double.eps),
    scaledgd_over_minimax_target_mse =
      sum((as.vector(target$W_S %*% as.vector(M_scaledgd)) - tau0)^2) /
      max(sum((as.vector(target$W_S %*% as.vector(M_minimax)) - tau0)^2), .Machine$double.eps),
    scaledgd_over_nuclear_target_mse =
      sum((as.vector(target$W_S %*% as.vector(M_scaledgd)) - tau0)^2) /
      max(sum((as.vector(target$W_S %*% as.vector(M_nuclear)) - tau0)^2), .Machine$double.eps),
    minimax_matrix_frob_error = matrix_minimax_error,
    nuclear_matrix_frob_error = matrix_nuclear_error,
    scaledgd_matrix_frob_error = matrix_scaledgd_error,
    used_matrix_frob_error = matrix_init_error
  )

  message("Computing target-specific oracle and pilot probabilities on the main pool...")
  WSPhi_true <- target$W_S %*% true_phi
  oracle_loading <- x_main %*% t(WSPhi_true)
  zeta <- oracle_loading * eps_main
  zeta_bar <- colMeans(zeta)
  h <- row_l2(oracle_loading)
  a <- abs(eps_main) * h

  pilot_phi <- phi_from_matrix(M_tilde, cfg$rank, pilot_cov$Sigmap_hat, pilot_cov$Sigmaq_hat)$Phi
  WSPhi_pilot <- target$W_S %*% pilot_phi
  pilot_resid_proxy <- as.vector(y_main - x_main %*% as.vector(M_tilde))
  pilot_loading <- x_main %*% t(WSPhi_pilot)
  h_pilot <- row_l2(pilot_loading)
  a_pilot_score <- abs(pilot_resid_proxy) * h_pilot

  pi_unif <- rep(1 / n_main, n_main)
  pi_oracle <- make_pi(a, cfg$defensive_gamma)
  pi_pilot <- make_pi(a_pilot_score, cfg$defensive_gamma)
  pi_list <- list("UNIF-MAT" = pi_unif, "ORACLE-PRO" = pi_oracle, "PILOT-PRO" = pi_pilot)

  A_unif <- sampling_objective(a, pi_unif)
  A_opt_no_defensive <- mean(a)^2
  A_oracle_used <- sampling_objective(a, pi_oracle)
  A_pilot_used <- sampling_objective(a, pi_pilot)
  Gamma_trace_unif <- conditional_gamma_trace(zeta, pi_unif)
  Gamma_trace_oracle <- conditional_gamma_trace(zeta, pi_oracle)
  Gamma_trace_pilot <- conditional_gamma_trace(zeta, pi_pilot)
  full_component_trace <- cfg$sigma^2 * sum(diag(target$W_S %*% true_phi %*% t(target$W_S)))

  probability_diagnostics <- rbind(
    probability_summary(pi_unif, "UNIF-MAT"),
    probability_summary(pi_oracle, "ORACLE-PRO"),
    probability_summary(pi_pilot, "PILOT-PRO")
  )
  probability_diagnostics$N <- cfg$n
  probability_diagnostics$N_main <- n_main
  probability_diagnostics$init_size <- cfg$init_size
  probability_diagnostics$design <- cfg$design
  probability_diagnostics$signal <- paste(cfg$signal, collapse = ",")
  probability_diagnostics$signal_label <- signal_label(cfg$signal)
  probability_diagnostics$target_label <- target$target_label
  probability_diagnostics$target_size <- target$target_size
  probability_diagnostics$selected_entries <- paste(target$entry_names, collapse = ";")

  score_diagnostics <- data.frame(
    N = cfg$n,
    N_main = n_main,
    init_size = cfg$init_size,
    design = cfg$design,
    signal = paste(cfg$signal, collapse = ","),
    signal_label = signal_label(cfg$signal),
    target_label = target$target_label,
    target_size = target$target_size,
    selected_entries = paste(target$entry_names, collapse = ";"),
    a_mean = mean(a),
    a_sd = sd(a),
    a_cv = cv(a),
    a_q01 = as.numeric(quantile(a, 0.01, names = FALSE)),
    a_median = median(a),
    a_q99 = as.numeric(quantile(a, 0.99, names = FALSE)),
    h_mean = mean(h),
    h_sd = sd(h),
    h_cv = cv(h),
    abs_eps_mean = mean(abs(eps_main)),
    abs_eps_sd = sd(abs(eps_main)),
    abs_eps_cv = cv(abs(eps_main)),
    pilot_score_mean = mean(a_pilot_score),
    pilot_score_sd = sd(a_pilot_score),
    pilot_score_cv = cv(a_pilot_score),
    cor_a_h = suppressWarnings(cor(a, h)),
    cor_a_abs_eps = suppressWarnings(cor(a, abs(eps_main))),
    cor_a_pilot_score = suppressWarnings(cor(a, a_pilot_score)),
    cor_h_h_pilot = suppressWarnings(cor(h, h_pilot)),
    cor_abs_eps_pilot_resid = suppressWarnings(cor(abs(eps_main), abs(pilot_resid_proxy))),
    A_UNIF = A_unif,
    A_OPT_no_defensive = A_opt_no_defensive,
    A_ORACLE_used = A_oracle_used,
    A_PILOT_used = A_pilot_used,
    leading_ratio_oracle_used_over_unif = A_oracle_used / A_unif,
    leading_ratio_pilot_used_over_unif = A_pilot_used / A_unif,
    leading_ratio_opt_no_def_over_unif = A_opt_no_defensive / A_unif,
    leading_gain_oracle_used_vs_unif = 1 - A_oracle_used / A_unif,
    leading_gain_pilot_used_vs_unif = 1 - A_pilot_used / A_unif,
    leading_gain_opt_no_def_vs_unif = 1 - A_opt_no_defensive / A_unif,
    leading_rmse_gain_oracle_used_vs_unif = 1 - sqrt(A_oracle_used / A_unif),
    full_component_trace = full_component_trace,
    Gamma_trace_UNIF = Gamma_trace_unif,
    Gamma_trace_ORACLE = Gamma_trace_oracle,
    Gamma_trace_PILOT = Gamma_trace_pilot,
    Gamma_trace_ratio_oracle_over_unif = Gamma_trace_oracle / Gamma_trace_unif,
    Gamma_trace_ratio_pilot_over_unif = Gamma_trace_pilot / Gamma_trace_unif,
    init_method_used = init_method_used,
    matrix_init_error = matrix_init_error,
    matrix_minimax_error = matrix_minimax_error,
    matrix_nuclear_error = matrix_nuclear_error,
    matrix_scaledgd_error = matrix_scaledgd_error,
    nuclear_over_minimax_matrix_frob = matrix_nuclear_error / max(matrix_minimax_error, .Machine$double.eps),
    scaledgd_over_minimax_matrix_frob = matrix_scaledgd_error / max(matrix_minimax_error, .Machine$double.eps),
    scaledgd_over_nuclear_matrix_frob = matrix_scaledgd_error / max(matrix_nuclear_error, .Machine$double.eps),
    nuclear_lambda = nuclear_fit$lambda,
    nuclear_iterations = nuclear_fit$iter,
    nuclear_converged = nuclear_fit$converged,
    seed_data = seed_data,
    seed_init = seed_init
  )

  full_data_diagnostics <- data.frame(
    N = cfg$n,
    N_main = n_main,
    init_size = cfg$init_size,
    design = cfg$design,
    signal = paste(cfg$signal, collapse = ","),
    signal_label = signal_label(cfg$signal),
    target_label = target$target_label,
    target_size = target$target_size,
    selected_entries = paste(target$entry_names, collapse = ";"),
    full_squared_error_tau0 = sum((tau_full - tau0)^2),
    full_rmse_tau0 = sqrt(sum((tau_full - tau0)^2)),
    full_estimator_type = "main_pool_two_step",
    init_method_used = init_method_used,
    matrix_init_error = matrix_init_error,
    matrix_minimax_error = matrix_minimax_error,
    matrix_nuclear_error = matrix_nuclear_error,
    matrix_scaledgd_error = matrix_scaledgd_error,
    nuclear_over_minimax_matrix_frob = matrix_nuclear_error / max(matrix_minimax_error, .Machine$double.eps),
    scaledgd_over_minimax_matrix_frob = matrix_scaledgd_error / max(matrix_minimax_error, .Machine$double.eps),
    scaledgd_over_nuclear_matrix_frob = matrix_scaledgd_error / max(matrix_nuclear_error, .Machine$double.eps),
    nuclear_lambda = nuclear_fit$lambda,
    nuclear_iterations = nuclear_fit$iter,
    nuclear_converged = nuclear_fit$converged,
    scaledgd_eta = scaledgd_fit$eta,
    scaledgd_eta_last = scaledgd_fit$eta_last,
    scaledgd_iterations = scaledgd_fit$iter,
    scaledgd_converged = scaledgd_fit$converged,
    scaledgd_objective = scaledgd_fit$objective,
    scaledgd_init_objective = scaledgd_fit$init_objective,
    scaledgd_backtrack_total = scaledgd_fit$backtrack_total
  )

  raw_results <- list()
  sampled_ids <- list()

  theory_records <- list()
  theory_id <- 1L
  for (r_sub in cfg$r_grid) {
    for (method in method_names(cfg)) {
      if (method == "UNIF-MAT") {
        A_used <- A_unif
        Gamma_trace_used <- Gamma_trace_unif
      } else if (method == "ORACLE-PRO") {
        A_used <- A_oracle_used
        Gamma_trace_used <- Gamma_trace_oracle
      } else if (method == "PILOT-PRO") {
        A_used <- A_pilot_used
        Gamma_trace_used <- Gamma_trace_pilot
      } else {
        A_used <- NA_real_
        Gamma_trace_used <- NA_real_
      }
      theory_conditional_mse <- if (is.finite(A_used)) A_used / r_sub else NA_real_
      theory_true_centered_mse <- if (is.finite(A_used) && cfg$add_full_component) {
        A_used / r_sub + full_component_trace / n_main
      } else if (is.finite(A_used)) {
        A_used / r_sub
      } else {
        NA_real_
      }
      theory_linear_conditional_mse <- if (is.finite(Gamma_trace_used)) Gamma_trace_used / r_sub else NA_real_
      theory_linear_true_centered_mse <- if (is.finite(Gamma_trace_used)) Gamma_trace_used / r_sub + sum((tau_full - tau0)^2) else NA_real_
      theory_records[[theory_id]] <- data.frame(
        r = r_sub,
        method = method,
        signal = paste(cfg$signal, collapse = ","),
        signal_label = signal_label(cfg$signal),
        target_label = target$target_label,
        target_size = target$target_size,
        selected_entries = paste(target$entry_names, collapse = ";"),
        A_used = A_used,
        A_UNIF = A_unif,
        A_OPT_no_defensive = A_opt_no_defensive,
        A_ORACLE_used = A_oracle_used,
        A_PILOT_used = A_pilot_used,
        Gamma_trace_used = Gamma_trace_used,
        Gamma_trace_UNIF = Gamma_trace_unif,
        Gamma_trace_ORACLE = Gamma_trace_oracle,
        Gamma_trace_PILOT = Gamma_trace_pilot,
        full_component_trace = full_component_trace,
        theory_conditional_mse = theory_conditional_mse,
        theory_true_centered_mse = theory_true_centered_mse,
        theory_linear_conditional_mse = theory_linear_conditional_mse,
        theory_linear_true_centered_mse = theory_linear_true_centered_mse
      )
      theory_id <- theory_id + 1L
    }
  }
  raw_theory <- do.call(rbind, theory_records)

  uniform_draw_methods <- intersect(c("UNIF-MAT", "UNIF-DB-SP", "NUCLEAR-R"), method_names(cfg))
  target_specific_methods <- intersect(c("ORACLE-PRO", "PILOT-PRO"), method_names(cfg))
  two_step_methods <- c("UNIF-MAT", "ORACLE-PRO", "PILOT-PRO")
  uniform_fit_cache <- shared$uniform_fit_cache
  if (is.null(uniform_fit_cache)) {
    stop("shared$uniform_fit_cache is missing. Please run prepare_uniform_fit_cache(cfg, shared) before target loop.")
  }

  fit_one_draw <- function(b) {
    local_results <- list()
    local_sampled_ids <- list()
    local_id <- 1L

    for (r_sub in cfg$r_grid) {
      # Target-invariant uniform-based methods are retrieved from the cache.
      # They are computed once per (draw, r) and then reused for every target W_S.
      for (method in uniform_draw_methods) {
        key <- paste(b, r_sub, method, sep = "__")
        uf <- uniform_fit_cache$fits[[key]]
        if (is.null(uf)) stop("Missing uniform fit cache entry: ", key)

        M_hat <- uf$M
        id <- uf$id
        id_local <- uf$id_local
        this_seed <- uf$seed

        if (cfg$save_sampled_ids) {
          local_sampled_ids[[paste(b, method, r_sub, sep = "__")]] <- id
        }

        tau_hat <- as.vector(target$W_S %*% as.vector(M_hat))
        err_tau0 <- tau_hat - tau0
        if (method %in% two_step_methods) {
          err_taufull <- tau_hat - tau_full
          squared_error_taufull <- sum(err_taufull * err_taufull)
        } else {
          squared_error_taufull <- NA_real_
        }

        if (method %in% names(pi_list)) {
          pi <- pi_list[[method]]
          eta_lin <- zeta[id_local, , drop = FALSE] / (n_main * pi[id_local])
          lin_delta <- colMeans(eta_lin) - zeta_bar
          tau_lin <- tau_full + as.vector(lin_delta)
          err_lin_tau0 <- tau_lin - tau0
          err_lin_taufull <- tau_lin - tau_full
          squared_error_linear_tau0 <- sum(err_lin_tau0 * err_lin_tau0)
          squared_error_linear_taufull <- sum(err_lin_taufull * err_lin_taufull)
        } else {
          squared_error_linear_tau0 <- NA_real_
          squared_error_linear_taufull <- NA_real_
        }

        local_results[[local_id]] <- data.frame(
          draw = b,
          r = r_sub,
          method = method,
          signal = paste(cfg$signal, collapse = ","),
          signal_label = signal_label(cfg$signal),
          target_label = target$target_label,
          target_size = target$target_size,
          selected_entries = paste(target$entry_names, collapse = ";"),
          squared_error_tau0 = sum(err_tau0 * err_tau0),
          squared_error_taufull = squared_error_taufull,
          squared_error_linear_tau0 = squared_error_linear_tau0,
          squared_error_linear_taufull = squared_error_linear_taufull,
          init_method_used = init_method_used,
          matrix_init_error = matrix_init_error,
          matrix_minimax_error = matrix_minimax_error,
          matrix_nuclear_error = matrix_nuclear_error,
          matrix_scaledgd_error = matrix_scaledgd_error,
          sampled_unique_n = uf$sampled_unique_n,
          sampled_duplicate_rate = uf$sampled_duplicate_rate,
          seed_data = seed_data,
          seed_sampling = this_seed,
          sampled_from_main_pool = TRUE,
          shared_uniform_ids = TRUE,
          uniform_fit_cached_across_targets = TRUE
        )
        local_id <- local_id + 1L
      }

      # Target-specific probability methods still have to be sampled and fitted for each target.
      for (method in target_specific_methods) {
        pi <- pi_list[[method]]
        this_seed <- seed_for_sampling(b, r_sub, method, cfg)
        set.seed(this_seed)
        id_local <- sample.int(n_main, r_sub, replace = TRUE, prob = pi)
        id <- main_idx[id_local]
        if (cfg$save_sampled_ids) {
          local_sampled_ids[[paste(b, method, r_sub, sep = "__")]] <- id
        }

        xs <- dat$xmat[id, , drop = FALSE]
        ys <- dat$y[id]
        fit <- matrix_two_step(
          xs, xmatT[id, , drop = FALSE], ys,
          ipw = 1 / pi[id_local], R_anchor = R_tilde,
          p = cfg$p, q = cfg$q, rank0 = cfg$rank
        )
        M_hat <- fit$M

        tau_hat <- as.vector(target$W_S %*% as.vector(M_hat))
        err_tau0 <- tau_hat - tau0
        err_taufull <- tau_hat - tau_full
        squared_error_taufull <- sum(err_taufull * err_taufull)

        eta_lin <- zeta[id_local, , drop = FALSE] / (n_main * pi[id_local])
        lin_delta <- colMeans(eta_lin) - zeta_bar
        tau_lin <- tau_full + as.vector(lin_delta)
        err_lin_tau0 <- tau_lin - tau0
        err_lin_taufull <- tau_lin - tau_full
        squared_error_linear_tau0 <- sum(err_lin_tau0 * err_lin_tau0)
        squared_error_linear_taufull <- sum(err_lin_taufull * err_lin_taufull)

        local_results[[local_id]] <- data.frame(
          draw = b,
          r = r_sub,
          method = method,
          signal = paste(cfg$signal, collapse = ","),
          signal_label = signal_label(cfg$signal),
          target_label = target$target_label,
          target_size = target$target_size,
          selected_entries = paste(target$entry_names, collapse = ";"),
          squared_error_tau0 = sum(err_tau0 * err_tau0),
          squared_error_taufull = squared_error_taufull,
          squared_error_linear_tau0 = squared_error_linear_tau0,
          squared_error_linear_taufull = squared_error_linear_taufull,
          init_method_used = init_method_used,
          matrix_init_error = matrix_init_error,
          matrix_minimax_error = matrix_minimax_error,
          matrix_nuclear_error = matrix_nuclear_error,
          matrix_scaledgd_error = matrix_scaledgd_error,
          sampled_unique_n = length(unique(id)),
          sampled_duplicate_rate = 1 - length(unique(id)) / length(id),
          seed_data = seed_data,
          seed_sampling = this_seed,
          sampled_from_main_pool = TRUE,
          shared_uniform_ids = FALSE,
          uniform_fit_cached_across_targets = FALSE
        )
        local_id <- local_id + 1L
      }
    }

    list(
      raw_results = do.call(rbind, local_results),
      sampled_ids = if (cfg$save_sampled_ids) local_sampled_ids else NULL
    )
  }

  message("Running conditional subsampling repetitions in parallel: B=", cfg$n_rep,
          ", n_cores=", cfg$n_cores)
  draw_ids <- seq_len(cfg$n_rep)

  if (cfg$n_cores > 1L && .Platform$OS.type != "windows") {
    chunks <- parallel::mclapply(
      draw_ids, fit_one_draw,
      mc.cores = cfg$n_cores,
      mc.preschedule = FALSE
    )
  } else if (cfg$n_cores > 1L && .Platform$OS.type == "windows") {
    warning(
      "Parallel execution on Windows would copy the full N x pq matrix to each worker. ",
      "Falling back to sequential execution. Use Linux/macOS for fork-based parallelism, ",
      "or set --n-cores=1 to suppress this warning."
    )
    chunks <- lapply(draw_ids, fit_one_draw)
  } else {
    chunks <- lapply(draw_ids, fit_one_draw)
  }

  raw_results <- do.call(rbind, lapply(chunks, function(z) z$raw_results))
  if (cfg$save_sampled_ids) {
    sampled_ids <- unlist(lapply(chunks, function(z) z$sampled_ids), recursive = FALSE)
  } else {
    sampled_ids <- NULL
  }

  summary_mse <- summarise_mse(raw_results, raw_theory)
  ratio_summary <- summarise_ratios(summary_mse)
  method_ratio_summary <- summarise_method_ratios(summary_mse)

  config_lines <- c(
    "Conditional large-N selected-entry target-sensitivity simulation: multi-method MSE comparison",
    "One shared full data set is generated once and reused across all targets.",
    "The pilot initialization sample is removed from the main subsampling pool.",
    "UNIF-MAT, UNIF-DB-SP, and NUCLEAR-R share the same uniform sampled IDs for each draw and r.",
    "Uniform-based estimators are cached once per (draw,r) and reused across all targets W_S.",
    "ORACLE-PRO probability uses true epsilon and true Phi; PILOT-PRO uses pilot residuals and pilot Phi.",
    "Initialization used for two-step fitting is controlled by --init-method; default is scaled_gd in this file.",
    "Nuclear initialization solves squared linear loss plus lambda*||M||_* on an initialization subset of size --init-size.",
    "ScaledGD initialization factorizes M = L R^T and runs covariance-calibrated moment initialization plus scaled gradient updates on the pilot sample.",
    "Artificial minimax initial estimator is still computed for diagnostics and comparison.",
    "Full-data-centered MSE is reported only for two-step subsampling methods: UNIF-MAT, ORACLE-PRO, and PILOT-PRO.",
    "For UNIF-DB-SP and NUCLEAR-R, tau_full-centered MSE is set to NA because these are different estimator classes.",
    "Linear oracle diagnostic is added only for probability-based two-step methods.",
    paste("N =", cfg$n),
    paste("N_main =", n_main),
    paste("p =", cfg$p),
    paste("q =", cfg$q),
    paste("rank =", cfg$rank),
    paste("sigma =", cfg$sigma),
    paste("signal =", paste(cfg$signal, collapse = ",")),
    paste("signal_label =", signal_label(cfg$signal)),
    paste("init_size =", cfg$init_size),
    paste("init_method =", cfg$init_method),
    paste("nuclear_lambda =", cfg$nuclear_lambda),
    paste("nuclear_lambda_factor =", cfg$nuclear_lambda_factor),
    paste("nuclear_maxit =", cfg$nuclear_maxit),
    paste("nuclear_tol =", cfg$nuclear_tol),
    paste("scaledgd_eta =", cfg$scaledgd_eta),
    paste("scaledgd_maxit =", cfg$scaledgd_maxit),
    paste("scaledgd_tol =", cfg$scaledgd_tol),
    paste("scaledgd_frob_tol =", cfg$scaledgd_frob_tol),
    paste("scaledgd_stop_rule =", cfg$scaledgd_stop_rule),
    paste("scaledgd_stop_window =", cfg$scaledgd_stop_window),
    paste("scaledgd_ridge =", cfg$scaledgd_ridge),
    paste("scaledgd_backtrack =", cfg$scaledgd_backtrack),
    paste("r_nuclear_lambda =", cfg$r_nuclear_lambda),
    paste("r_nuclear_lambda_factor =", cfg$r_nuclear_lambda_factor),
    paste("run_unif_db_sp =", cfg$run_unif_db_sp),
    paste("run_r_nuclear =", cfg$run_r_nuclear),
    paste("nuclear_lambda_used =", nuclear_fit$lambda),
    paste("nuclear_iterations =", nuclear_fit$iter),
    paste("nuclear_converged =", nuclear_fit$converged),
    paste("scaledgd_objective =", scaledgd_fit$objective),
    paste("scaledgd_init_objective =", scaledgd_fit$init_objective),
    paste("scaledgd_iterations =", scaledgd_fit$iter),
    paste("scaledgd_converged =", scaledgd_fit$converged),
    paste("scaledgd_rel_objective =", scaledgd_fit$rel_objective),
    paste("scaledgd_rel_frob =", scaledgd_fit$rel_frob),
    paste("scaledgd_window_rel_objective =", scaledgd_fit$window_rel_objective),
    paste("scaledgd_window_rel_frob =", scaledgd_fit$window_rel_frob),
    paste("scaledgd_stop_streak =", scaledgd_fit$stop_streak),
    paste("scaledgd_eta_last =", scaledgd_fit$eta_last),
    paste("scaledgd_backtrack_total =", scaledgd_fit$backtrack_total),
    paste("r_grid =", paste(cfg$r_grid, collapse = ",")),
    paste("n_rep_conditional_draws =", cfg$n_rep),
    paste("n_cores =", cfg$n_cores),
    paste("available_cores =", available_cores),
    paste("design =", cfg$design),
    paste("rho_p =", cfg$rho_p),
    paste("rho_q =", cfg$rho_q),
    paste("cov_scale =", cfg$cov_scale),
    "Scaled design covariance: Sigmap = cov_scale * Sigmap_raw and Sigmaq = cov_scale * Sigmaq_raw; set cov_scale by --cov-scale=a.",
    paste("defensive_gamma =", cfg$defensive_gamma),
    paste("add_full_component =", cfg$add_full_component),
    paste("save_probabilities =", cfg$save_probabilities),
    paste("save_sampled_ids =", cfg$save_sampled_ids),
    paste("target_label =", target$target_label),
    paste("target_size =", target$target_size),
    paste("selected_entries =", paste(target$entry_names, collapse = ",")),
    paste("target_order =", cfg$target_order),
    paste("target_mode =", cfg$target_mode),
    paste("row_index =", cfg$row_index),
    paste("col_index =", cfg$col_index),
    paste("selected_entries_argument =", cfg$selected_entries),
    paste("target_size_grid =", paste(cfg$target_size_grid, collapse = ",")),  # unused in this single-target version
    paste("seed_data =", seed_data),
    paste("seed_init =", seed_init),
    paste("out_root =", cfg$out_root)
  )

  writeLines(config_lines, file.path(cfg$out_root, "config.txt"))
  write.csv(raw_results, file.path(cfg$out_root, "raw_empirical_mse_by_draw.csv"), row.names = FALSE)
  write.csv(raw_theory, file.path(cfg$out_root, "raw_theory_mse.csv"), row.names = FALSE)
  write.csv(probability_diagnostics, file.path(cfg$out_root, "probability_diagnostics.csv"), row.names = FALSE)
  write.csv(score_diagnostics, file.path(cfg$out_root, "score_theory_diagnostics.csv"), row.names = FALSE)
  write.csv(full_data_diagnostics, file.path(cfg$out_root, "full_data_diagnostics.csv"), row.names = FALSE)
  write.csv(initial_diagnostics, file.path(cfg$out_root, "initial_estimator_diagnostics.csv"), row.names = FALSE)
  write.csv(initial_comparison, file.path(cfg$out_root, "initial_estimator_comparison.csv"), row.names = FALSE)
  write.csv(summary_mse, file.path(cfg$out_root, "summary_mse_empirical_and_theory.csv"), row.names = FALSE)
  write.csv(ratio_summary, file.path(cfg$out_root, "summary_oracle_unif_ratios.csv"), row.names = FALSE)
  write.csv(method_ratio_summary, file.path(cfg$out_root, "summary_method_ratios_vs_unif.csv"), row.names = FALSE)

  cache <- list(
    cfg_small = cfg[setdiff(names(cfg), c("out_root"))],
    seed_data = seed_data,
    seed_init = seed_init,
    init_idx = init_idx,
    main_idx = main_idx,
    R_tilde = R_tilde,
    M_tilde_fnorm_error = matrix_init_error,
    scaledgd_fit_summary = list(
      objective = scaledgd_fit$objective,
      init_objective = scaledgd_fit$init_objective,
      iter = scaledgd_fit$iter,
      converged = scaledgd_fit$converged,
      rel_objective = scaledgd_fit$rel_objective,
      rel_frob = scaledgd_fit$rel_frob,
      window_rel_objective = scaledgd_fit$window_rel_objective,
      window_rel_frob = scaledgd_fit$window_rel_frob,
      frob_tol = scaledgd_fit$frob_tol,
      stop_rule = scaledgd_fit$stop_rule,
      stop_window = scaledgd_fit$stop_window,
      stop_streak = scaledgd_fit$stop_streak,
      eta_last = scaledgd_fit$eta_last,
      backtrack_total = scaledgd_fit$backtrack_total
    ),
    tau0 = tau0,
    tau_full = tau_full,
    full_data_diagnostics = full_data_diagnostics,
    initial_diagnostics = initial_diagnostics,
    initial_comparison = initial_comparison,
    score_diagnostics = score_diagnostics,
    probability_diagnostics = probability_diagnostics,
    theory = list(
      A_UNIF = A_unif,
      A_OPT_no_defensive = A_opt_no_defensive,
      A_ORACLE_used = A_oracle_used,
      A_PILOT_used = A_pilot_used,
      Gamma_trace_UNIF = Gamma_trace_unif,
      Gamma_trace_ORACLE = Gamma_trace_oracle,
      Gamma_trace_PILOT = Gamma_trace_pilot,
      full_component_trace = full_component_trace
    ),
    sampled_ids = if (cfg$save_sampled_ids) sampled_ids else NULL
  )
  if (cfg$save_probabilities) {
    cache$probabilities_main_pool <- pi_list
  }
  saveRDS(cache, file.path(cfg$out_root, "conditional_cache.rds"), compress = "gzip")

  plot_two_method_mse(summary_mse, file.path(cfg$out_root, "figures", "fig_empirical_mse_tau0.pdf"),
                      "empirical_MSE_tau0_mean", "Selected-entry MSE centered at tau0")
  plot_two_method_mse(summary_mse, file.path(cfg$out_root, "figures", "fig_empirical_mse_taufull.pdf"),
                      "empirical_MSE_taufull_mean", "Two-step selected-entry MSE centered at tau_full")
  plot_empirical_vs_theory(summary_mse, file.path(cfg$out_root, "figures", "fig_empirical_vs_theory_tau0.pdf"),
                           "empirical_MSE_tau0_mean", "theory_true_centered_MSE", "MSE centered at tau0: empirical vs theory")
  plot_empirical_vs_theory(summary_mse, file.path(cfg$out_root, "figures", "fig_empirical_vs_theory_taufull.pdf"),
                           "empirical_MSE_taufull_mean", "theory_conditional_MSE", "Two-step MSE centered at tau_full: empirical vs conditional theory")
  plot_empirical_vs_theory(summary_mse, file.path(cfg$out_root, "figures", "fig_linear_vs_theory_taufull.pdf"),
                           "empirical_linear_MSE_taufull_mean", "theory_linear_conditional_MSE", "Linear diagnostic centered at tau_full: empirical vs theory")
  plot_empirical_vs_theory(summary_mse, file.path(cfg$out_root, "figures", "fig_linear_vs_theory_tau0.pdf"),
                           "empirical_linear_MSE_tau0_mean", "theory_linear_true_centered_MSE", "Linear diagnostic centered at tau0: empirical vs theory")
  plot_ratio(ratio_summary, file.path(cfg$out_root, "figures", "fig_oracle_unif_ratio_tau0.pdf"),
             "empirical_tau0_ratio_oracle_over_unif", "theory_true_centered_ratio_oracle_over_unif", "MSE ratio centered at tau0")
  plot_ratio(ratio_summary, file.path(cfg$out_root, "figures", "fig_oracle_unif_ratio_taufull.pdf"),
             "empirical_taufull_ratio_oracle_over_unif", "theory_conditional_ratio_oracle_over_unif", "MSE ratio centered at tau_full")
  plot_ratio(ratio_summary, file.path(cfg$out_root, "figures", "fig_linear_ratio_taufull.pdf"),
             "empirical_linear_taufull_ratio_oracle_over_unif", "theory_linear_conditional_ratio_oracle_over_unif", "Linear MSE ratio centered at tau_full")
  plot_ratio(ratio_summary, file.path(cfg$out_root, "figures", "fig_linear_ratio_tau0.pdf"),
             "empirical_linear_tau0_ratio_oracle_over_unif", "theory_linear_true_centered_ratio_oracle_over_unif", "Linear MSE ratio centered at tau0")

  message("Finished. Main summary: ", file.path(cfg$out_root, "summary_mse_empirical_and_theory.csv"))
  invisible(list(
    raw_results = raw_results,
    raw_theory = raw_theory,
    probability_diagnostics = probability_diagnostics,
    score_diagnostics = score_diagnostics,
    full_data_diagnostics = full_data_diagnostics,
    initial_diagnostics = initial_diagnostics,
    initial_comparison = initial_comparison,
    summary_mse = summary_mse,
    ratio_summary = ratio_summary,
    method_ratio_summary = method_ratio_summary
  ))
}

run_single_target <- function(cfg) {
  top_out <- cfg$out_root
  dir.create(top_out, recursive = TRUE, showWarnings = FALSE)

  shared <- prepare_shared_context(cfg)
  shared$uniform_fit_cache <- prepare_uniform_fit_cache(cfg, shared)
  write.csv(shared$uniform_fit_cache$index, file.path(top_out, "uniform_fit_cache_index.csv"), row.names = FALSE)

  all_summary <- list()
  all_ratios <- list()
  all_method_ratios <- list()
  all_scores <- list()
  all_full <- list()
  all_initial <- list()
  all_initial_comparison <- list()

  for (selected_pos in cfg$target_grid) {  # only one target in this simplified version
    cfg_i <- cfg
    cfg_i$selected_pos <- selected_pos
    cfg_i$target_name <- target_label(selected_pos)
    cfg_i$out_root <- file.path(top_out, cfg_i$target_name)
    message("============================================================")
    message("Running rank-", cfg_i$rank, " signal=", paste(cfg_i$signal, collapse = ","),
            " target: ", cfg_i$target_name)
    message("Selected entries: ", paste(apply(selected_pos, 1, function(z) paste0("M[", z[1], ",", z[2], "]")), collapse = ","))
    message("Output directory: ", cfg_i$out_root)
    res <- run_all(cfg_i, shared = shared)

    sm <- res$summary_mse
    sm$target_label <- cfg_i$target_name
    sm$target_size <- nrow(selected_pos)
    sm$selected_entries <- paste(apply(selected_pos, 1, function(z) paste0("M[", z[1], ",", z[2], "]")), collapse = ";")
    all_summary[[cfg_i$target_name]] <- sm

    rt <- res$ratio_summary
    rt$target_label <- cfg_i$target_name
    rt$target_size <- nrow(selected_pos)
    rt$selected_entries <- paste(apply(selected_pos, 1, function(z) paste0("M[", z[1], ",", z[2], "]")), collapse = ";")
    all_ratios[[cfg_i$target_name]] <- rt

    mrt <- res$method_ratio_summary
    mrt$target_label <- cfg_i$target_name
    mrt$target_size <- nrow(selected_pos)
    mrt$selected_entries <- paste(apply(selected_pos, 1, function(z) paste0("M[", z[1], ",", z[2], "]")), collapse = ";")
    all_method_ratios[[cfg_i$target_name]] <- mrt

    all_scores[[cfg_i$target_name]] <- res$score_diagnostics
    all_full[[cfg_i$target_name]] <- res$full_data_diagnostics
    all_initial[[cfg_i$target_name]] <- res$initial_diagnostics
    all_initial_comparison[[cfg_i$target_name]] <- res$initial_comparison
  }

  combined_summary <- do.call(rbind, all_summary)
  combined_ratios <- do.call(rbind, all_ratios)
  combined_method_ratios <- do.call(rbind, all_method_ratios)
  combined_scores <- do.call(rbind, all_scores)
  combined_full <- do.call(rbind, all_full)
  combined_initial <- do.call(rbind, all_initial)
  combined_initial_comparison <- do.call(rbind, all_initial_comparison)

  write.csv(combined_summary, file.path(top_out, "combined_summary_mse_empirical_and_theory.csv"), row.names = FALSE)
  write.csv(combined_ratios, file.path(top_out, "combined_summary_oracle_unif_ratios.csv"), row.names = FALSE)
  write.csv(combined_method_ratios, file.path(top_out, "combined_summary_method_ratios_vs_unif.csv"), row.names = FALSE)
  write.csv(combined_scores, file.path(top_out, "combined_score_theory_diagnostics.csv"), row.names = FALSE)
  write.csv(combined_full, file.path(top_out, "combined_full_data_diagnostics.csv"), row.names = FALSE)
  write.csv(combined_initial, file.path(top_out, "combined_initial_estimator_diagnostics.csv"), row.names = FALSE)
  write.csv(combined_initial_comparison, file.path(top_out, "combined_initial_estimator_comparison.csv"), row.names = FALSE)

  target_check <- combined_scores[, c(
    "signal", "signal_label", "target_label", "target_size", "selected_entries",
    "N", "N_main", "init_size",
    "A_UNIF", "A_ORACLE_used", "A_PILOT_used",
    "leading_ratio_oracle_used_over_unif", "leading_ratio_pilot_used_over_unif",
    "Gamma_trace_UNIF", "Gamma_trace_ORACLE", "Gamma_trace_PILOT",
    "Gamma_trace_ratio_oracle_over_unif", "Gamma_trace_ratio_pilot_over_unif",
    "cor_a_pilot_score", "cor_h_h_pilot", "cor_abs_eps_pilot_resid",
    "full_component_trace", "matrix_init_error"
  )]
  write.csv(target_check, file.path(top_out, "combined_target_theory_check.csv"), row.names = FALSE)

  pdf(file.path(top_out, "fig_target_theory_ratio_check.pdf"), width = 8, height = 5)
  on.exit(dev.off(), add = TRUE)
  op <- par(mar = c(8, 4, 3, 1))
  labels <- target_check$target_label
  y <- target_check$leading_ratio_oracle_used_over_unif
  plot(seq_along(y), y, type = "b", xaxt = "n", xlab = "", ylab = "ORACLE / UNIF leading ratio",
       main = "First-order theory for the selected-entry target", ylim = range(y, finite = TRUE))
  axis(1, at = seq_along(labels), labels = labels, las = 2, cex.axis = 0.65)
  grid()
  par(op)

  message("Single selected-entry target finished. Combined outputs written to: ", top_out)
  invisible(list(
    combined_summary = combined_summary,
    combined_ratios = combined_ratios,
    combined_method_ratios = combined_method_ratios,
    combined_scores = combined_scores,
    combined_full = combined_full,
    combined_initial = combined_initial,
    combined_initial_comparison = combined_initial_comparison,
    target_check = target_check
  ))
}

run_single_target(cfg)
