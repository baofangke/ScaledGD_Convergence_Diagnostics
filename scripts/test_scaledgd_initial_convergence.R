# Trace scaledgd_initial convergence for the default maxit=500 setting.
#
# This standalone smoke test mirrors the ScaledGD initialization loop used by
# scaledgd_initial() in
# simulation2_rank2_signal64_single_target3_scaledgd_multi_methods_unifdbsp_uniformcache_conditional_parallel_linear.R.
# It records the moment/SVD initial point and every accepted ScaledGD iteration.
#
# Outputs:
#   scaledgd_initial_convergence_history.csv
#   scaledgd_initial_convergence_summary.csv
#   config.txt
#   figures/fig_scaledgd_initial_convergence.pdf
#   figures/fig_scaledgd_initial_convergence.png
#
# Example:
#   Rscript R/test_scaledgd_initial_convergence.R
#   Rscript R/test_scaledgd_initial_convergence.R --sample-size=800 --scaledgd-maxit=500

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
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
  if (length(out) == 0L || any(!is.finite(out))) stop("Invalid numeric grid: ", value)
  out
}

mat_fnorm <- function(A) sqrt(sum(A * A))

safe_solve <- function(A, ridge = 1e-9) {
  A <- as.matrix(A)
  scale <- mean(abs(diag(A)))
  if (!is.finite(scale) || scale <= 0) scale <- 1
  solve(A + ridge * scale * diag(nrow(A)))
}

symmetrize <- function(A) (A + t(A)) / 2

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

make_design_objects <- function(cfg) {
  Sigmap_raw <- if (cfg$design == "isotropic") diag(cfg$p) else ar1_cov(cfg$p, cfg$rho_p)
  Sigmaq_raw <- if (cfg$design == "isotropic") diag(cfg$q) else exchangeable_cov(cfg$q, cfg$rho_q)
  list(
    Sigmap = cfg$cov_scale * Sigmap_raw,
    Sigmaq = cfg$cov_scale * Sigmaq_raw,
    Sigmap_raw = Sigmap_raw,
    Sigmaq_raw = Sigmaq_raw
  )
}

make_xmat <- function(n, p, q, Sigmap, Sigmaq) {
  z <- matrix(rnorm(n * p * q), n, p * q)
  K <- kronecker(chol(Sigmaq), chol(Sigmap))
  z %*% K
}

generate_data <- function(n, truth, cfg, Sigmap, Sigmaq) {
  beta0 <- as.vector(truth$M0)
  xmat <- make_xmat(n, cfg$p, cfg$q, Sigmap, Sigmaq)
  eps <- cfg$sigma * rnorm(n)
  y <- as.vector(xmat %*% beta0 + eps)
  list(xmat = xmat, y = y, eps = eps)
}

rank_truncate <- function(M, rank0) {
  S <- svd(M, nu = rank0, nv = rank0)
  U <- S$u[, seq_len(rank0), drop = FALSE]
  V <- S$v[, seq_len(rank0), drop = FALSE]
  D <- diag(S$d[seq_len(rank0)], nrow = rank0)
  list(M = U %*% D %*% t(V), U = U, V = V, D = D, d = S$d)
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
  list(
    S1_raw = regularize_cov(S1_raw, ridge = ridge),
    S2_raw = regularize_cov(S2_raw, ridge = ridge),
    avefro = avefro,
    Sigmap_hat = regularize_cov(S1_raw / scale, ridge = ridge),
    Sigmaq_hat = regularize_cov(S2_raw / scale, ridge = ridge)
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

scaledgd_initial_trace <- function(xmat, y, p, q, rank0, M_truth = NULL,
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
  M_old <- Lfac %*% t(Rfac)
  truth_norm <- if (is.null(M_truth)) NA_real_ else max(mat_fnorm(M_truth), .Machine$double.eps)

  history <- vector("list", maxit + 1L)
  history[[1L]] <- data.frame(
    iter = 0L,
    objective = obj_old,
    rel_objective_change = NA_real_,
    step_frob = NA_real_,
    rel_step_frob = NA_real_,
    truth_frob = if (is.null(M_truth)) NA_real_ else mat_fnorm(M_old - M_truth),
    rel_truth_frob = if (is.null(M_truth)) NA_real_ else mat_fnorm(M_old - M_truth) / truth_norm,
    truth_matrix_mse = if (is.null(M_truth)) NA_real_ else mean((M_old - M_truth)^2),
    eta_step = NA_real_,
    backtrack_count = NA_integer_,
    backtrack_total = 0L,
    accepted_no_move = FALSE,
    stop_condition = FALSE,
    stop_streak = 0L,
    window_rel_objective = NA_real_,
    window_rel_frob = NA_real_,
    rank_1e8 = sum(svd(M_old, nu = 0, nv = 0)$d > 1e-8),
    rank_1e4 = sum(svd(M_old, nu = 0, nv = 0)$d > 1e-4),
    converged = FALSE
  )

  converged <- FALSE
  iter <- 0L
  eta_last <- eta
  backtrack_count_total <- 0L
  stop_streak <- 0L
  window_rel_change <- NA_real_
  window_rel_frob <- NA_real_
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
    accepted_no_move <- FALSE

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
        accepted_no_move <- TRUE
      } else {
        step <- step / 2
        n_backtrack <- n_backtrack + 1L
      }
    }
    backtrack_count_total <- backtrack_count_total + n_backtrack

    M_new <- Lnew %*% t(Rnew)
    rel_change <- abs(obj_old - obj_new) / (1 + abs(obj_old))
    step_frob <- mat_fnorm(M_new - M_current)
    rel_step_frob <- step_frob / (1 + mat_fnorm(M_current))
    obj_history[iter + 1L] <- obj_new
    rel_frob_history[iter] <- rel_step_frob
    svals <- svd(M_new, nu = 0, nv = 0)$d
    if (stop_rule == "window") {
      window_check <- scaledgd_window_stop(
        obj_history, rel_frob_history, iter, stop_window, tol, frob_tol
      )
      stop_condition <- window_check$stop
      window_rel_change <- window_check$rel_objective
      window_rel_frob <- window_check$rel_frob
    } else {
      stop_condition <- scaledgd_should_stop(stop_rule, rel_change, rel_step_frob, tol, frob_tol)
      window_rel_change <- NA_real_
      window_rel_frob <- NA_real_
    }
    if (stop_condition) {
      stop_streak <- stop_streak + 1L
    } else {
      stop_streak <- 0L
    }

    history[[iter + 1L]] <- data.frame(
      iter = iter,
      objective = obj_new,
      rel_objective_change = rel_change,
      step_frob = step_frob,
      rel_step_frob = rel_step_frob,
      truth_frob = if (is.null(M_truth)) NA_real_ else mat_fnorm(M_new - M_truth),
      rel_truth_frob = if (is.null(M_truth)) NA_real_ else mat_fnorm(M_new - M_truth) / truth_norm,
      truth_matrix_mse = if (is.null(M_truth)) NA_real_ else mean((M_new - M_truth)^2),
      eta_step = eta_last,
      backtrack_count = n_backtrack,
      backtrack_total = backtrack_count_total,
      accepted_no_move = accepted_no_move,
      stop_condition = stop_condition,
      stop_streak = stop_streak,
      window_rel_objective = window_rel_change,
      window_rel_frob = window_rel_frob,
      rank_1e8 = sum(svals > 1e-8),
      rank_1e4 = sum(svals > 1e-4),
      converged = FALSE
    )

    Lfac <- Lnew
    Rfac <- Rnew
    if ((stop_rule == "window" && stop_condition) ||
        (stop_rule != "window" && stop_streak >= stop_window)) {
      converged <- TRUE
      history[[iter + 1L]]$converged <- TRUE
      obj_old <- obj_new
      break
    }
    obj_old <- obj_new
    M_old <- M_new
  }

  history <- do.call(rbind, history[seq_len(iter + 1L)])
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
    iter = iter,
    converged = converged,
    eta = eta,
    eta_last = eta_last,
    backtrack_total = backtrack_count_total,
    ridge = ridge,
    history = history,
    stop_window = stop_window,
    stop_streak = stop_streak,
    window_rel_objective = window_rel_change,
    window_rel_frob = window_rel_frob,
    rank_1e8 = sum(svals > 1e-8),
    rank_1e4 = sum(svals > 1e-4),
    singular_values = svals
  )
}

plot_convergence <- function(history, cfg, fit, pdf_file, png_file) {
  draw <- function() {
    op <- par(mfrow = c(2, 3), mar = c(4.2, 4.4, 3.0, 1.0))
    on.exit(par(op), add = TRUE)

    z <- history[history$iter > 0, , drop = FALSE]
    stop_iter <- if (fit$converged) fit$iter else NA_integer_
    main_suffix <- paste0("m=", cfg$sample_size, ", maxit=", cfg$scaledgd_maxit)

    plot(z$iter, pmax(z$rel_step_frob, .Machine$double.eps),
         type = "l", log = "y", lwd = 2, col = "#2166AC",
         xlab = "iteration", ylab = expression("relative " * "||" * M[k] - M[k-1] * "||"[F]),
         main = paste("Relative Frobenius step,", main_suffix))
    if (any(is.finite(z$window_rel_frob))) {
      lines(z$iter, pmax(z$window_rel_frob, .Machine$double.eps),
            lwd = 2, lty = 2, col = "#D6604D")
      legend("topright", bty = "n", lwd = c(2, 2), lty = c(1, 2),
             col = c("#2166AC", "#D6604D"),
             legend = c("single step", paste0(cfg$scaledgd_stop_window, "-step mean")),
             cex = 0.8)
    }
    abline(h = cfg$scaledgd_frob_tol, lty = 3, col = "#B2182B")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(history$iter, history$truth_frob,
         type = "l", lwd = 2, col = "#1B7837",
         xlab = "iteration", ylab = expression("||" * M[k] - M[0] * "||"[F]),
         main = "Frobenius error to truth")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(history$iter, history$objective,
         type = "l", lwd = 2, col = "#762A83",
         xlab = "iteration", ylab = "objective",
         main = "Least-squares objective")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(z$iter, pmax(z$rel_objective_change, .Machine$double.eps),
         type = "l", log = "y", lwd = 2, col = "#B35806",
         xlab = "iteration", ylab = "relative objective change",
         main = "Original stopping criterion")
    if (any(is.finite(z$window_rel_objective))) {
      lines(z$iter, pmax(z$window_rel_objective, .Machine$double.eps),
            lwd = 2, lty = 2, col = "#4393C3")
      legend("topright", bty = "n", lwd = c(2, 2), lty = c(1, 2),
             col = c("#B35806", "#4393C3"),
             legend = c("single step", paste0(cfg$scaledgd_stop_window, "-step mean")),
             cex = 0.8)
    }
    abline(h = cfg$scaledgd_tol, lty = 3, col = "#B2182B")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(z$iter, z$eta_step,
         type = "s", lwd = 2, col = "#4D4D4D",
         xlab = "iteration", ylab = "accepted step size",
         main = "Backtracking step size")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")

    plot(z$iter, z$backtrack_count,
         type = "h", lwd = 2, col = "#5AAE61",
         xlab = "iteration", ylab = "backtracks",
         main = "Backtracks per iteration")
    grid()
    if (is.finite(stop_iter)) abline(v = stop_iter, lty = 2, col = "#B2182B")
  }

  pdf(pdf_file, width = 11, height = 7)
  draw()
  dev.off()

  png(png_file, width = 1650, height = 1050, res = 150)
  draw()
  dev.off()
}

cfg <- list(
  seed = as.integer(get_arg("seed", "20260527")),
  p = as.integer(get_arg("p", "10")),
  q = as.integer(get_arg("q", "10")),
  rank = as.integer(get_arg("rank", "2")),
  sigma = as.numeric(get_arg("sigma", "0.5")),
  signal = parse_numeric_grid(get_arg("signal", "6.0,4.4")),
  design = get_arg("design", "kronecker"),
  rho_p = as.numeric(get_arg("rho-p", "0.8")),
  rho_q = as.numeric(get_arg("rho-q", "0.7")),
  cov_scale = as.numeric(get_arg("cov-scale", "1.0")),
  sample_size = as.integer(get_arg("sample-size", get_arg("m", "200"))),
  scaledgd_eta = as.numeric(get_arg("scaledgd-eta", "0.5")),
  scaledgd_maxit = as.integer(get_arg("scaledgd-maxit", "500")),
  scaledgd_tol = as.numeric(get_arg("scaledgd-tol", "1e-6")),
  scaledgd_frob_tol = as.numeric(get_arg("scaledgd-frob-tol", "1e-3")),
  scaledgd_stop_rule = get_arg("scaledgd-stop-rule", "objective"),
  scaledgd_stop_window = as.integer(get_arg("scaledgd-stop-window", "1")),
  scaledgd_ridge = as.numeric(get_arg("scaledgd-ridge", "1e-6")),
  scaledgd_backtrack = get_bool_arg("scaledgd-backtrack", TRUE),
  reference_maxit = as.integer(get_arg("reference-maxit", "2000")),
  reference_tol = as.numeric(get_arg("reference-tol", "0")),
  out_root = normalizePath(
    get_arg("out-root", file.path("0608result", "scaledgd_initial_convergence_test")),
    winslash = "/", mustWork = FALSE
  )
)

if (!cfg$design %in% c("isotropic", "kronecker")) stop("--design must be isotropic or kronecker.")
if (cfg$p < 1L || cfg$q < 1L) stop("p and q must be positive integers.")
if (cfg$rank < 1L || cfg$rank > min(cfg$p, cfg$q)) stop("Invalid rank.")
if (length(cfg$signal) != cfg$rank) stop("--signal length must equal --rank.")
if (!is.finite(cfg$sigma) || cfg$sigma < 0) stop("--sigma must be nonnegative.")
if (!is.finite(cfg$cov_scale) || cfg$cov_scale <= 0) stop("--cov-scale must be positive.")
if (cfg$sample_size < 2L) stop("--sample-size must be at least 2.")
if (!is.finite(cfg$scaledgd_eta) || cfg$scaledgd_eta <= 0) stop("--scaledgd-eta must be positive.")
if (cfg$scaledgd_maxit < 1L) stop("--scaledgd-maxit must be positive.")
if (!is.finite(cfg$scaledgd_tol) || cfg$scaledgd_tol <= 0) stop("--scaledgd-tol must be positive.")
if (!is.finite(cfg$scaledgd_frob_tol) || cfg$scaledgd_frob_tol <= 0) {
  stop("--scaledgd-frob-tol must be positive.")
}
if (!cfg$scaledgd_stop_rule %in% c("objective", "frob", "both", "either", "window")) {
  stop("--scaledgd-stop-rule must be one of: objective, frob, both, either, window.")
}
if (!is.finite(cfg$scaledgd_stop_window) || cfg$scaledgd_stop_window < 1L) {
  stop("--scaledgd-stop-window must be positive.")
}
if (!is.finite(cfg$scaledgd_ridge) || cfg$scaledgd_ridge < 0) stop("--scaledgd-ridge must be nonnegative.")
if (cfg$reference_maxit < cfg$scaledgd_maxit) {
  stop("--reference-maxit must be at least --scaledgd-maxit.")
}
if (!is.finite(cfg$reference_tol) || cfg$reference_tol < 0) {
  stop("--reference-tol must be nonnegative.")
}

dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$out_root, "figures"), recursive = TRUE, showWarnings = FALSE)

message("Generating one test data set...")
set.seed(cfg$seed)
truth <- make_truth(cfg$p, cfg$q, cfg$rank, cfg$signal)
design <- make_design_objects(cfg)
dat <- generate_data(cfg$sample_size, truth, cfg, design$Sigmap, design$Sigmaq)
pilot_cov <- pilot_covariance_estimates(dat$xmat, cfg$p, cfg$q, ridge = cfg$scaledgd_ridge)

message("Tracing scaledgd_initial with maxit=", cfg$scaledgd_maxit, "...")
fit <- scaledgd_initial_trace(
  dat$xmat, dat$y, p = cfg$p, q = cfg$q, rank0 = cfg$rank,
  M_truth = truth$M0,
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

ref_fit <- NULL
if (cfg$reference_maxit > cfg$scaledgd_maxit) {
  message("Computing reference run with maxit=", cfg$reference_maxit,
          " and forced tol=0 objective stopping...")
  ref_fit <- scaledgd_initial_trace(
    dat$xmat, dat$y, p = cfg$p, q = cfg$q, rank0 = cfg$rank,
    M_truth = truth$M0,
    S1_raw = pilot_cov$S1_raw, S2_raw = pilot_cov$S2_raw, avefro = pilot_cov$avefro,
    eta = cfg$scaledgd_eta,
    maxit = cfg$reference_maxit,
    tol = 0,
    frob_tol = cfg$scaledgd_frob_tol,
    stop_rule = "objective",
    stop_window = 1L,
    ridge = cfg$scaledgd_ridge,
    backtrack = cfg$scaledgd_backtrack
  )
}

history <- fit$history
final <- history[nrow(history), , drop = FALSE]
init <- history[history$iter == 0L, , drop = FALSE]
truth_norm <- max(mat_fnorm(truth$M0), .Machine$double.eps)
minimax_rate <- cfg$sigma * sqrt(cfg$rank * (cfg$p + cfg$q) / cfg$sample_size)

summary_df <- data.frame(
  seed = cfg$seed,
  sample_size = cfg$sample_size,
  p = cfg$p,
  q = cfg$q,
  rank = cfg$rank,
  sigma = cfg$sigma,
  signal = paste(cfg$signal, collapse = ","),
  design = cfg$design,
  cov_scale = cfg$cov_scale,
  scaledgd_eta = cfg$scaledgd_eta,
  scaledgd_maxit = cfg$scaledgd_maxit,
  scaledgd_tol = cfg$scaledgd_tol,
  scaledgd_frob_tol = cfg$scaledgd_frob_tol,
  scaledgd_stop_rule = cfg$scaledgd_stop_rule,
  scaledgd_stop_window = cfg$scaledgd_stop_window,
  scaledgd_ridge = cfg$scaledgd_ridge,
  scaledgd_backtrack = cfg$scaledgd_backtrack,
  iter = fit$iter,
  converged = fit$converged,
  hit_maxit = fit$iter >= cfg$scaledgd_maxit && !fit$converged,
  init_objective = fit$init_objective,
  final_objective = fit$objective,
  objective_drop = fit$init_objective - fit$objective,
  final_rel_objective_change = final$rel_objective_change,
  window_rel_objective = fit$window_rel_objective,
  window_rel_frob = fit$window_rel_frob,
  init_truth_frob = init$truth_frob,
  init_rel_truth_frob = init$rel_truth_frob,
  final_step_frob = final$step_frob,
  final_rel_step_frob = final$rel_step_frob,
  final_truth_frob = mat_fnorm(fit$M - truth$M0),
  final_rel_truth_frob = mat_fnorm(fit$M - truth$M0) / truth_norm,
  minimax_rate = minimax_rate,
  final_over_minimax_rate = mat_fnorm(fit$M - truth$M0) / minimax_rate,
  final_sq_over_minimax_sq = mat_fnorm(fit$M - truth$M0)^2 / minimax_rate^2,
  eta_last = fit$eta_last,
  backtrack_total = fit$backtrack_total,
  stop_streak = fit$stop_streak,
  rank_1e8 = fit$rank_1e8,
  rank_1e4 = fit$rank_1e4,
  reference_maxit = cfg$reference_maxit,
  reference_tol = cfg$reference_tol,
  reference_iter = if (is.null(ref_fit)) NA_integer_ else ref_fit$iter,
  reference_converged = if (is.null(ref_fit)) NA else ref_fit$converged,
  reference_truth_frob = if (is.null(ref_fit)) NA_real_ else mat_fnorm(ref_fit$M - truth$M0),
  reference_rel_truth_frob = if (is.null(ref_fit)) NA_real_ else mat_fnorm(ref_fit$M - truth$M0) / truth_norm,
  reference_over_minimax_rate = if (is.null(ref_fit)) NA_real_ else mat_fnorm(ref_fit$M - truth$M0) / minimax_rate,
  frob_diff_vs_reference = if (is.null(ref_fit)) NA_real_ else mat_fnorm(fit$M - ref_fit$M),
  rel_frob_diff_vs_reference = if (is.null(ref_fit)) NA_real_ else mat_fnorm(fit$M - ref_fit$M) / truth_norm,
  objective_gap_vs_reference = if (is.null(ref_fit)) NA_real_ else fit$objective - ref_fit$objective
)

history_file <- file.path(cfg$out_root, "scaledgd_initial_convergence_history.csv")
summary_file <- file.path(cfg$out_root, "scaledgd_initial_convergence_summary.csv")
config_file <- file.path(cfg$out_root, "config.txt")
fig_pdf <- file.path(cfg$out_root, "figures", "fig_scaledgd_initial_convergence.pdf")
fig_png <- file.path(cfg$out_root, "figures", "fig_scaledgd_initial_convergence.png")

write.csv(history, history_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)
config_lines <- vapply(
  names(cfg),
  function(nm) paste(nm, paste(cfg[[nm]], collapse = ","), sep = " = "),
  character(1)
)
writeLines(config_lines, config_file)
plot_convergence(history, cfg, fit, fig_pdf, fig_png)

message("Done.")
message("Summary: ", summary_file)
message("History: ", history_file)
message("Figure: ", fig_pdf)
print(summary_df)
