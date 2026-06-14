base_dir <- "C:/Users/Ke Baofang/Desktop/matrix-subsampling-simulations/0608result/p30_q25_covscale02_eta1_sample_size_sweep_n300_2000"
sample_sizes <- c(300, 400, 500, 700, 1000, 1500, 2000)

read_one <- function(n) {
  file <- file.path(base_dir, paste0("m", n), "scaledgd_initial_convergence_summary.csv")
  z <- read.csv(file)
  data.frame(
    sample_size = z$sample_size,
    p = z$p,
    q = z$q,
    rank = z$rank,
    sigma = z$sigma,
    cov_scale = z$cov_scale,
    eta = z$scaledgd_eta,
    maxit = z$scaledgd_maxit,
    iter = z$iter,
    converged = z$converged,
    hit_maxit = z$hit_maxit,
    final_objective = z$final_objective,
    final_truth_frob = z$final_truth_frob,
    minimax_rate_standard = z$minimax_rate,
    ratio_standard = z$final_over_minimax_rate,
    minimax_rate_covscale_adjusted = z$minimax_rate / z$cov_scale,
    ratio_covscale_adjusted = z$final_truth_frob / (z$minimax_rate / z$cov_scale),
    window_rel_objective = z$window_rel_objective,
    window_rel_frob = z$window_rel_frob,
    eta_last = z$eta_last,
    backtrack_total = z$backtrack_total
  )
}

summary_df <- do.call(rbind, lapply(sample_sizes, read_one))
summary_file <- file.path(base_dir, "summary_scaledgd_eta1_sample_size_sweep_n300_2000.csv")
write.csv(summary_df, summary_file, row.names = FALSE)

fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png_file <- file.path(fig_dir, "fig_scaledgd_frob_error_vs_sample_size_n300_2000.png")
pdf_file <- file.path(fig_dir, "fig_scaledgd_frob_error_vs_sample_size_n300_2000.pdf")

fit <- lm(log(final_truth_frob) ~ log(sample_size), data = summary_df)
slope <- unname(coef(fit)[2])

anchor_n <- 500
anchor_err <- summary_df$final_truth_frob[summary_df$sample_size == anchor_n]
ref_err <- anchor_err * sqrt(anchor_n / summary_df$sample_size)

draw <- function() {
  op <- par(mfrow = c(1, 2), mar = c(4.6, 4.8, 3.2, 1.0))
  on.exit(par(op), add = TRUE)

  plot(summary_df$sample_size, summary_df$final_truth_frob,
       type = "b", lwd = 2, pch = 19, col = "#0072B2",
       xlab = "sample size n", ylab = expression("||" * hat(M) - M[0] * "||"[F]),
       main = "Frobenius error, n=300-2000")
  lines(summary_df$sample_size, ref_err, lty = 2, lwd = 2, col = "#666666")
  grid()
  legend("topright", bty = "n", cex = 0.82,
         legend = c("observed", paste0("C/sqrt(n), anchored at n=", anchor_n)),
         col = c("#0072B2", "#666666"),
         pch = c(19, NA), lty = c(1, 2), lwd = c(2, 2))

  plot(summary_df$sample_size, summary_df$final_truth_frob,
       log = "xy", type = "b", lwd = 2, pch = 19, col = "#0072B2",
       xlab = "sample size n (log)", ylab = expression("||" * hat(M) - M[0] * "||"[F] * " (log)"),
       main = paste0("Log-log slope: ", sprintf("%.3f", slope)))
  lines(summary_df$sample_size, ref_err, lty = 2, lwd = 2, col = "#666666")
  grid()
  legend("topright", bty = "n", cex = 0.82,
         legend = c("observed", "slope -1/2 reference"),
         col = c("#0072B2", "#666666"),
         pch = c(19, NA), lty = c(1, 2), lwd = c(2, 2))
}

png(png_file, width = 1650, height = 750, res = 150)
draw()
dev.off()

pdf(pdf_file, width = 11, height = 5)
draw()
dev.off()

message("Wrote: ", summary_file)
message("Wrote: ", png_file)
message("Wrote: ", pdf_file)
