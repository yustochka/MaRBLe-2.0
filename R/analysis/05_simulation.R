# R/analysis/05_simulation.R
#
# Gene-set level bootstrap stability simulation for threshold validation.
# Uses phyper()-based hypergeometric ORA, mirroring Pipeline A step06 logic.
# No external packages required.
#
# Simulation A — Pathway coverage vs bootstrap stability
#   Question: At what pathway observation fraction does TAU=0.70 break down?
#   Varies: fraction of full pathway genes present in the dataset (COVERAGE_GRID).
#   Fixed:  n_universe=5000, DE_FRAC=0.10, SIGNAL_FRAC=0.30 (of observed pathway).
#   Mechanism: at low coverage only a small fragment of the pathway is observed by
#              ORA; the effective pathway size shrinks, reducing test power even
#              when the true signal exists.
#   Connects to: LOW_COVERAGE classification (med_frac_in_matrix < 0.20) and
#                dataset_2 where frac_in_matrix ~ 0.10-0.15 for all final pathways.
#
# Simulation B — Signal strength vs bootstrap stability
#   Question: How do DE fraction and pathway signal overlap jointly affect stability?
#   Varies: DE fraction × signal overlap fraction (2D grid).
#   Fixed:  n_universe=5000 (representative of normal datasets).
#   Connects to: dataset_3 (DE frac ~9.2%, BORDERLINE) vs dataset_0 (DE frac ~19.5%, HEALTHY).
#
# Real dataset values are overlaid on plots when threshold_policy.csv is available.
#
# Run from project root:
#   Rscript R/analysis/05_simulation.R

ROOT <- getwd()

# =============================================================================
# Parameters — all tunable here
# =============================================================================

B_BOOT    <- 200L    # bootstrap iterations per stability estimate (matches Pipeline A)
SUB_FRAC  <- 0.85    # subsample fraction (matches Pipeline A)
P_CUT     <- 0.05    # enrichment significance threshold
N_REPS    <- 20L     # simulation replicates per parameter combo (for variance estimation)
BASE_SEED <- 42L     # reproducibility seed

# Simulation A: pathway coverage sweep
SIM_A_N_UNIVERSE    <- 5000L   # fixed healthy universe
SIM_A_N_FULL_PATHWAY <- 200L   # full annotated pathway size (typical WP/GO-BP)
SIM_A_DE_FRAC       <- 0.10    # DE fraction (10% of universe)
SIM_A_SIGNAL_FRAC   <- 0.30    # fraction of OBSERVED pathway genes that are DE
SIM_A_COVERAGE_GRID <- c(0.05, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00)

# Simulation B: signal strength sweep (unchanged)
SIM_B_N_UNIVERSE       <- 5000L
SIM_B_DE_FRAC_GRID     <- c(0.05, 0.08, 0.10, 0.15, 0.20, 0.25)
SIM_B_SIGNAL_FRAC_GRID <- c(0.10, 0.20, 0.30, 0.40)

# TAU reference values and LOW_COVERAGE fraction boundary
TAU70         <- 0.70
TAU50         <- 0.50
LOW_COV_FRAC  <- 0.20   # classification boundary: frac_in_matrix < 0.20 = LOW_COVERAGE

# Paths
OUT_DIR    <- file.path(ROOT, "results", "analysis")
FIG_DIR    <- file.path(ROOT, "results", "analysis", "figures")
POLICY_CSV <- file.path(OUT_DIR, "threshold_policy.csv")

# =============================================================================
# Core simulation function
# =============================================================================
#
# Universe layout (integer indices 1..n_universe):
#   1 .. n_signal                     : signal genes  (pathway AND DE)
#   n_signal+1 .. n_pathway           : pathway-only  (not DE)
#   n_pathway+1 .. n_pathway + n_bg   : background DE (not in pathway)
#   n_pathway + n_bg + 1 .. n_universe: inert genes
#
# n_pathway here = n_observed_pathway (what ORA sees, not the full annotation).
# Bootstrap: sample 85% of universe; recompute hypergeometric ORA; collect hit.
# Stability = fraction of B bootstraps where p < P_CUT.

simulate_one_stability <- function(n_universe, n_de, n_pathway, n_signal,
                                    B, sub_frac, p_cut, seed) {
  set.seed(seed)
  n_bg     <- n_de - n_signal
  sub_size <- floor(sub_frac * n_universe)
  n_pw_bg  <- n_pathway + n_bg   # boundary: pathway + background DE genes

  hits <- vapply(seq_len(B), function(b) {
    idx <- sample.int(n_universe, size = sub_size)

    n_sub_path <- sum(idx <= n_pathway)
    n_sub_both <- sum(idx <= n_signal)                                # pathway AND DE
    n_sub_de   <- n_sub_both + sum(idx > n_pathway & idx <= n_pw_bg) # all DE in sub

    if (n_sub_path == 0L || n_sub_de == 0L) return(0L)

    # P(X >= n_sub_both): one-sided hypergeometric enrichment test
    p_val <- phyper(n_sub_both - 1L,
                    n_sub_path, sub_size - n_sub_path,
                    n_sub_de, lower.tail = FALSE)
    as.integer(p_val < p_cut)
  }, integer(1L))

  mean(hits)
}

# =============================================================================
# Simulation A: pathway coverage sweep
# =============================================================================

run_sim_A_coverage <- function(n_universe, n_full_pathway, de_frac, signal_frac,
                                coverage_grid, B, sub_frac, p_cut, n_reps, base_seed) {
  cat("\n=== Simulation A: pathway coverage vs stability ===\n")
  cat(sprintf("  n_universe=%d | n_full_pathway=%d | DE_frac=%.2f | signal_frac=%.2f\n",
              n_universe, n_full_pathway, de_frac, signal_frac))
  cat(sprintf("  B=%d | reps=%d\n\n", B, n_reps))

  n_de_total <- round(de_frac * n_universe)  # total DE genes; fixed across coverage grid

  rows <- vector("list", length(coverage_grid) * n_reps)
  idx  <- 0L

  for (cov in coverage_grid) {
    n_obs_pathway <- max(1L, round(n_full_pathway * cov))   # observed (measured) pathway size
    n_signal      <- max(0L, round(n_obs_pathway * signal_frac))
    n_bg          <- n_de_total - n_signal

    # Sanity guards
    if (n_signal > n_de_total) {
      message(sprintf("  SKIP cov=%.2f: n_signal (%d) > n_de (%d)", cov, n_signal, n_de_total))
      next
    }
    if (n_bg > n_universe - n_obs_pathway) {
      message(sprintf("  SKIP cov=%.2f: n_bg (%d) > available non-pathway (%d)",
                      cov, n_bg, n_universe - n_obs_pathway))
      next
    }

    cat(sprintf("  coverage=%.2f  n_obs_pathway=%3d  n_signal=%2d  n_bg_de=%d\n",
                cov, n_obs_pathway, n_signal, n_bg))

    for (r in seq_len(n_reps)) {
      idx  <- idx + 1L
      # Seed stable per (coverage, rep): adding grid values won't reseed old combos
      seed <- base_seed + round(cov * 1000L) * 10L + r
      stab <- simulate_one_stability(n_universe, n_de_total,
                                      n_obs_pathway, n_signal,
                                      B, sub_frac, p_cut, seed)
      rows[[idx]] <- data.frame(
        n_universe     = n_universe,
        n_full_pathway = n_full_pathway,
        coverage       = cov,
        n_obs_pathway  = n_obs_pathway,
        de_fraction    = de_frac,
        n_de           = n_de_total,
        signal_frac    = signal_frac,
        n_signal       = n_signal,
        replicate      = r,
        stability      = stab,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows[seq_len(idx)])
}

# =============================================================================
# Simulation B: signal strength sweep (unchanged)
# =============================================================================

run_sim_B <- function(n_universe, de_frac_grid, signal_frac_grid, n_pathway,
                       B, sub_frac, p_cut, n_reps, base_seed) {
  cat("\n=== Simulation B: signal strength vs stability ===\n")
  cat(sprintf("  n_universe=%d | B=%d | reps=%d\n\n", n_universe, B, n_reps))

  n_combos <- length(de_frac_grid) * length(signal_frac_grid)
  rows     <- vector("list", n_combos * n_reps)
  idx      <- 0L

  for (de_frac in de_frac_grid) {
    n_de <- round(de_frac * n_universe)

    for (sig_frac in signal_frac_grid) {
      n_signal <- round(sig_frac * n_pathway)

      if (n_signal > n_de) {
        message(sprintf("  SKIP de=%.2f sig=%.2f: n_signal (%d) > n_de (%d)",
                        de_frac, sig_frac, n_signal, n_de))
        next
      }
      n_bg <- n_de - n_signal
      if (n_bg > n_universe - n_pathway) {
        message(sprintf("  SKIP de=%.2f sig=%.2f: n_bg exceeds available",
                        de_frac, sig_frac))
        next
      }

      cat(sprintf("  de_frac=%.2f  n_de=%4d  sig_frac=%.2f  n_signal=%2d\n",
                  de_frac, n_de, sig_frac, n_signal))

      for (r in seq_len(n_reps)) {
        idx  <- idx + 1L
        seed <- base_seed + 50000L + round(de_frac * 1000) * 100L +
                round(sig_frac * 100) * 10L + r
        stab <- simulate_one_stability(n_universe, n_de, n_pathway, n_signal,
                                        B, sub_frac, p_cut, seed)
        rows[[idx]] <- data.frame(
          n_universe  = n_universe,
          de_fraction = de_frac,
          n_de        = n_de,
          n_pathway   = n_pathway,
          signal_frac = sig_frac,
          n_signal    = n_signal,
          replicate   = r,
          stability   = stab,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows[seq_len(idx)])
}

# =============================================================================
# Summary: mean / SD / pass-rate per parameter group
# =============================================================================

summarise_sim <- function(df, group_vars) {
  key_df <- unique(df[, group_vars, drop = FALSE])
  key_df <- key_df[do.call(order, lapply(group_vars, function(v) key_df[[v]])), ,
                    drop = FALSE]
  rownames(key_df) <- NULL

  result <- lapply(seq_len(nrow(key_df)), function(i) {
    mask <- rep(TRUE, nrow(df))
    for (v in group_vars)
      mask <- mask & abs(df[[v]] - key_df[i, v]) < 1e-9
    sub      <- df$stability[mask]
    row      <- key_df[i, , drop = FALSE]
    row$n_reps            <- length(sub)
    row$stability_mean    <- round(mean(sub), 4)
    row$stability_sd      <- round(sd(sub),  4)
    row$frac_passes_tau70 <- round(mean(sub >= TAU70), 4)
    row$frac_passes_tau50 <- round(mean(sub >= TAU50), 4)
    row
  })
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  out
}

# =============================================================================
# Validation checks
# =============================================================================

validate_output <- function(df, label) {
  cat(sprintf("\nValidating %s ...\n", label))
  n_na  <- sum(is.na(df$stability))
  n_oor <- sum(df$stability < 0 | df$stability > 1, na.rm = TRUE)
  cat(sprintf("  rows=%d  NA=%d  out-of-range=%d\n", nrow(df), n_na, n_oor))
  if (n_na > 0 || n_oor > 0) warning(label, ": unexpected values detected")
}

# =============================================================================
# Plot A: stability vs pathway coverage fraction
# =============================================================================

plot_sim_A_coverage <- function(summ, real_data, out_file) {
  summ <- summ[order(summ$coverage), ]
  x    <- summ$coverage
  ymn  <- summ$stability_mean
  ylo  <- pmax(0, ymn - summ$stability_sd)
  yhi  <- pmin(1, ymn + summ$stability_sd)

  n_full   <- unique(summ$n_full_pathway)[1]
  n_obs_at <- round(n_full * x)   # observed pathway sizes along x-axis

  png(out_file, width = 860, height = 560, res = 96)
  par(mar = c(5.5, 5, 4.5, 2))

  plot(x, ymn, type = "n",
       xlab = "Pathway coverage (fraction of annotated members measured)",
       ylab = "Bootstrap stability (mean \u00b1 SD)",
       main = sprintf(
         "Sim A: Stability vs Pathway Coverage\n(N_universe=%d, N_full_pathway=%d, DE_frac=%.2f, signal_frac=%.2f, B=%d)",
         unique(summ$n_universe)[1], n_full,
         unique(summ$de_fraction)[1], unique(summ$signal_frac)[1], B_BOOT),
       xlim = c(0, 1.05), ylim = c(0, 1.1), las = 1,
       xaxt = "n")

  # Custom x-axis showing both fraction and observed pathway size
  axis(1L, at = x,
       labels = sprintf("%.2f\n(n=%d)", x, n_obs_at),
       cex.axis = 0.72, padj = 0.5)

  # SD band
  polygon(c(x, rev(x)), c(yhi, rev(ylo)),
          col = adjustcolor("steelblue", alpha.f = 0.20), border = NA)

  # Simulated curve
  lines(x, ymn, col = "steelblue", lwd = 2.5)
  points(x, ymn, pch = 16, col = "steelblue", cex = 1.5)

  # TAU reference lines
  abline(h = TAU70, lty = 2, col = "red3",      lwd = 2)
  abline(h = TAU50, lty = 3, col = "darkorange", lwd = 2)
  text(1.03, TAU70 + 0.02, "TAU=0.70", cex = 0.7, col = "red3",      adj = c(1, 0))
  text(1.03, TAU50 + 0.02, "TAU=0.50", cex = 0.7, col = "darkorange", adj = c(1, 0))

  # LOW_COVERAGE boundary at frac=0.20
  abline(v = LOW_COV_FRAC, lty = 4, col = "gray40", lwd = 1.8)
  text(LOW_COV_FRAC, 1.07,
       sprintf(" cov=%.2f\n LOW_COVERAGE\n boundary", LOW_COV_FRAC),
       adj = c(0, 1), cex = 0.68, col = "gray40")

  # Real dataset overlay using med_frac_in_matrix
  if (!is.null(real_data) && "med_frac_in_matrix" %in% colnames(real_data)) {
    rd <- real_data[!is.na(real_data$med_frac_in_matrix), ]
    if (nrow(rd) > 0) {
      point_cols <- c("black", "gray30", "gray50", "gray65")
      pchs_rd    <- c(17L, 18L, 15L, 16L)
      for (i in seq_len(nrow(rd))) {
        xr  <- rd$med_frac_in_matrix[i]
        yr  <- rd$mean_stability_go[i]
        pc  <- if (i <= length(point_cols)) point_cols[i] else "black"
        ph  <- if (i <= length(pchs_rd)) pchs_rd[i] else 17L
        points(xr, yr, pch = ph, col = pc, cex = 2.0)
        adj <- if (xr > 0.5) c(1.15, 0.5) else c(-0.15, 0.5)
        text(xr, yr,
             labels = sprintf("%s\n(%s)", rd$dataset[i], rd$classification[i]),
             cex = 0.62, col = pc, adj = adj)
      }
    }
  }

  legend("topleft",
         c("Simulated (mean \u00b1 SD)", "TAU = 0.70", "TAU = 0.50",
           "Coverage = 0.20 (LOW_COVERAGE)", "Real datasets (mean GO stability)"),
         col = c("steelblue", "red3", "darkorange", "gray40", "black"),
         lty = c(1, 2, 3, 4, NA), pch = c(16, NA, NA, NA, 17),
         lwd = c(2.5, 2, 2, 1.8, NA), cex = 0.78, bty = "n")

  dev.off()
  cat(sprintf("  Wrote: %s\n", out_file))
}

# =============================================================================
# Plot B: stability vs DE fraction (line plot) — unchanged
# =============================================================================

plot_sim_B_lines <- function(summ, real_data, out_file) {
  sig_fracs <- sort(unique(summ$signal_frac))
  de_fracs  <- sort(unique(summ$de_fraction))
  n_pathway <- unique(summ$n_pathway)[1]

  line_cols <- c("steelblue", "darkorange", "seagreen4", "firebrick3")
  pchs      <- c(16L, 17L, 15L, 18L)
  leg_labs  <- sprintf("Signal %.0f%% (S=%d)", sig_fracs * 100,
                       round(sig_fracs * n_pathway))

  png(out_file, width = 820, height = 540, res = 96)
  par(mar = c(5, 5, 4, 2))

  plot(range(de_fracs), c(0, 1.05), type = "n",
       xlab = "DE fraction (proportion of universe that is DE)",
       ylab = "Bootstrap stability (mean \u00b1 SD)",
       main = sprintf("Sim B: Stability vs DE Fraction\n(n_universe=%d, pathway size=%d, B=200)",
                      unique(summ$n_universe)[1], n_pathway),
       las = 1)

  for (j in seq_along(sig_fracs)) {
    sub <- summ[abs(summ$signal_frac - sig_fracs[j]) < 1e-9, ]
    sub <- sub[order(sub$de_fraction), ]
    ylo <- pmax(0, sub$stability_mean - sub$stability_sd)
    yhi <- pmin(1, sub$stability_mean + sub$stability_sd)

    polygon(c(sub$de_fraction, rev(sub$de_fraction)), c(yhi, rev(ylo)),
            col = adjustcolor(line_cols[j], 0.15), border = NA)
    lines(sub$de_fraction, sub$stability_mean, col = line_cols[j], lwd = 2)
    points(sub$de_fraction, sub$stability_mean, pch = pchs[j],
           col = line_cols[j], cex = 1.3)
  }

  abline(h = TAU70, lty = 2, col = "red3", lwd = 2)

  # Real dataset vertical markers
  if (!is.null(real_data)) {
    vline_cols <- c("gray25", "gray50", "gray65", "gray40")
    for (i in seq_len(nrow(real_data))) {
      rd  <- real_data[i, ]
      vc  <- if (i <= length(vline_cols)) vline_cols[i] else "gray40"
      abline(v = rd$de_fraction, lty = 3, col = vc, lwd = 1.5)
      y_label <- 0.03 + (i - 1) * 0.08
      text(rd$de_fraction, y_label,
           labels = sprintf(" %s\n (%s)", rd$dataset, rd$classification),
           adj = c(0, 0), cex = 0.62, col = vc)
    }
  }

  legend("topleft",
         c(leg_labs, "TAU = 0.70"),
         col = c(line_cols[seq_along(sig_fracs)], "red3"),
         lty = c(rep(1L, length(sig_fracs)), 2L),
         pch = c(pchs[seq_along(sig_fracs)], NA),
         lwd = c(rep(2L, length(sig_fracs)), 2L),
         cex = 0.8, bty = "n")

  dev.off()
  cat(sprintf("  Wrote: %s\n", out_file))
}

# =============================================================================
# Plot B heatmap: stability[de_fraction, signal_frac] — unchanged
# =============================================================================

plot_sim_B_heatmap <- function(summ, out_file) {
  de_fracs  <- sort(unique(summ$de_fraction))
  sig_fracs <- sort(unique(summ$signal_frac))

  mat <- matrix(NA_real_, nrow = length(de_fracs), ncol = length(sig_fracs))
  for (i in seq_along(de_fracs)) {
    for (j in seq_along(sig_fracs)) {
      sub <- summ[abs(summ$de_fraction - de_fracs[i]) < 1e-9 &
                   abs(summ$signal_frac  - sig_fracs[j]) < 1e-9, ]
      if (nrow(sub) == 1L) mat[i, j] <- sub$stability_mean
    }
  }

  pal <- colorRampPalette(c("#f7f7f7", "#d1e5f0", "#4393c3", "#2166ac", "#053061"))(100)

  png(out_file, width = 660, height = 580, res = 96)
  par(mar = c(5, 5.5, 5, 2))

  image(seq_along(de_fracs), seq_along(sig_fracs), mat,
        col = pal, zlim = c(0, 1),
        xaxt = "n", yaxt = "n",
        xlab = "DE fraction",
        ylab = "Signal fraction (prop. of pathway genes that are DE)",
        main = sprintf("Sim B: Stability Heatmap (n_universe=%d)\n",
                       unique(summ$n_universe)[1]))

  mtext("Light = low stability  \u2192  Dark blue = high stability | Red contour = TAU = 0.70",
        side = 3, cex = 0.75, line = 0.2, col = "gray30")

  axis(1L, at = seq_along(de_fracs),  labels = sprintf("%.2f", de_fracs))
  axis(2L, at = seq_along(sig_fracs), labels = sprintf("%.2f", sig_fracs), las = 1L)

  for (i in seq_along(de_fracs)) {
    for (j in seq_along(sig_fracs)) {
      v <- mat[i, j]
      if (!is.na(v)) {
        tc <- if (v > 0.55) "white" else "gray20"
        text(i, j, sprintf("%.2f", v), cex = 0.88, col = tc, font = 2L)
      }
    }
  }

  if (any(mat >= TAU70, na.rm = TRUE) && any(mat < TAU70, na.rm = TRUE)) {
    contour(seq_along(de_fracs), seq_along(sig_fracs), mat,
            levels = TAU70, add = TRUE, col = "red3", lwd = 2.5,
            labcex = 1.0, drawlabels = TRUE)
  } else {
    cat("  [heatmap] All cells on one side of TAU=0.70 — contour omitted\n")
  }

  dev.off()
  cat(sprintf("  Wrote: %s\n", out_file))
}

# =============================================================================
# Main
# =============================================================================

cat("=== 05_simulation.R ===\n")
cat(sprintf("ROOT: %s\n", ROOT))
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load real data for overlay -----------------------------------------------
real_data <- NULL
if (file.exists(POLICY_CSV)) {
  rd <- read.csv(POLICY_CSV, stringsAsFactors = FALSE)
  keep_cols <- c("dataset", "n_universe", "de_fraction",
                 "mean_stability_go", "med_frac_in_matrix", "classification")
  keep_cols <- keep_cols[keep_cols %in% colnames(rd)]
  real_data <- rd[, keep_cols]
  real_data <- real_data[!is.na(real_data$n_universe), ]
  cat(sprintf("\nLoaded real data for overlay: %d datasets\n", nrow(real_data)))
} else {
  cat("\n[info] threshold_policy.csv not found — plots will not include real dataset overlay\n")
}

# --- Simulation A: pathway coverage -------------------------------------------
raw_A  <- run_sim_A_coverage(
  SIM_A_N_UNIVERSE, SIM_A_N_FULL_PATHWAY,
  SIM_A_DE_FRAC, SIM_A_SIGNAL_FRAC,
  SIM_A_COVERAGE_GRID,
  B_BOOT, SUB_FRAC, P_CUT, N_REPS, BASE_SEED)

validate_output(raw_A, "Simulation A")
write.csv(raw_A, file.path(OUT_DIR, "sim_A_coverage_stability.csv"), row.names = FALSE)
cat(sprintf("\nWrote: results/analysis/sim_A_coverage_stability.csv (%d rows)\n", nrow(raw_A)))

summ_A <- summarise_sim(raw_A, "coverage")
write.csv(summ_A, file.path(OUT_DIR, "sim_A_coverage_summary.csv"), row.names = FALSE)
cat(sprintf("Wrote: results/analysis/sim_A_coverage_summary.csv (%d rows)\n", nrow(summ_A)))

# --- Simulation B: signal strength --------------------------------------------
raw_B  <- run_sim_B(SIM_B_N_UNIVERSE, SIM_B_DE_FRAC_GRID, SIM_B_SIGNAL_FRAC_GRID,
                     50L, B_BOOT, SUB_FRAC, P_CUT, N_REPS, BASE_SEED)

validate_output(raw_B, "Simulation B")
write.csv(raw_B, file.path(OUT_DIR, "sim_B_signal_stability.csv"), row.names = FALSE)
cat(sprintf("\nWrote: results/analysis/sim_B_signal_stability.csv (%d rows)\n", nrow(raw_B)))

summ_B <- summarise_sim(raw_B, c("de_fraction", "signal_frac"))
write.csv(summ_B, file.path(OUT_DIR, "sim_B_summary.csv"), row.names = FALSE)
cat(sprintf("Wrote: results/analysis/sim_B_summary.csv (%d rows)\n", nrow(summ_B)))

# --- Figures ------------------------------------------------------------------
cat("\nGenerating figures...\n")
plot_sim_A_coverage(summ_A, real_data,
                    file.path(FIG_DIR, "fig_sim_A_stability_vs_coverage.png"))
plot_sim_B_lines(summ_B, real_data,
                 file.path(FIG_DIR, "fig_sim_B_stability_vs_de_fraction.png"))
plot_sim_B_heatmap(summ_B,
                    file.path(FIG_DIR, "fig_sim_B_heatmap.png"))

# --- Console summary: Simulation A -------------------------------------------
cat("\n=== SUMMARY: Simulation A (pathway coverage) ===\n")
cat(sprintf("  n_universe=%d | n_full_pathway=%d | DE_frac=%.2f | signal_frac=%.2f\n\n",
            SIM_A_N_UNIVERSE, SIM_A_N_FULL_PATHWAY, SIM_A_DE_FRAC, SIM_A_SIGNAL_FRAC))
cat(sprintf("  %-10s  %-12s  %-14s  %-12s  %s\n",
            "coverage", "n_obs_path", "stability_mean", "passes_tau70", "passes_tau50"))

for (i in seq_len(nrow(summ_A))) {
  r         <- summ_A[i, ]
  n_obs     <- round(SIM_A_N_FULL_PATHWAY * r$coverage)
  marker    <- if (abs(r$coverage - LOW_COV_FRAC) < 1e-9) " <- LOW_COVERAGE boundary" else ""
  cat(sprintf("  %-10.2f  %-12d  %-14.3f  %-12s  %s%s\n",
              r$coverage, n_obs,
              r$stability_mean,
              sprintf("%.0f%%", r$frac_passes_tau70 * 100),
              sprintf("%.0f%%", r$frac_passes_tau50 * 100),
              marker))
}

# Where does TAU=0.70 first cross 50% pass rate (rising from low coverage)?
above50 <- summ_A[summ_A$frac_passes_tau70 >= 0.50, ]
if (nrow(above50) > 0) {
  cat(sprintf("\n  TAU=0.70 first reaches >=50%% pass rate at coverage = %.2f\n",
              min(above50$coverage)))
}
above90 <- summ_A[summ_A$frac_passes_tau70 >= 0.90, ]
if (nrow(above90) > 0) {
  cat(sprintf("  TAU=0.70 first reaches >=90%% pass rate at coverage = %.2f\n",
              min(above90$coverage)))
}

# TAU=0.50 recovery: where does it first reach >=90%?
rec50 <- summ_A[summ_A$frac_passes_tau50 >= 0.90, ]
if (nrow(rec50) > 0) {
  cat(sprintf("  TAU=0.50 first reaches >=90%% pass rate at coverage = %.2f\n",
              min(rec50$coverage)))
}

# --- Console summary: Simulation B -------------------------------------------
cat("\n=== SUMMARY: Simulation B ===\n")
cat(sprintf("  %-10s  %-11s  %-14s  %s\n",
            "de_frac", "signal_frac", "stability_mean", "passes_tau70"))
for (i in seq_len(nrow(summ_B))) {
  r <- summ_B[i, ]
  cat(sprintf("  %-10.2f  %-11.2f  %-14.3f  %.0f%%\n",
              r$de_fraction, r$signal_frac, r$stability_mean,
              r$frac_passes_tau70 * 100))
}

cross_B <- summ_B[summ_B$frac_passes_tau70 > 0 & summ_B$frac_passes_tau70 < 1, ]
if (nrow(cross_B) > 0) {
  cat("\n  TAU=0.70 boundary zone (0% < pass_rate < 100%):\n")
  for (i in seq_len(nrow(cross_B))) {
    r <- cross_B[i, ]
    cat(sprintf("    de_frac=%.2f  signal_frac=%.2f  stability=%.3f\n",
                r$de_fraction, r$signal_frac, r$stability_mean))
  }
}

cat("\nDone.\n")
