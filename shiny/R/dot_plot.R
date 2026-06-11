# shiny/R/dot_plot.R
# Dot plot for the Final Shortlist.
# Requires ggplot2 (available in the project renv).
#
# Usage: make_dot_plot(df, pipeline, selected_idx = 1L)
#   df           — data frame from load_shortlist()$data (standardised columns)
#   pipeline     — "A" or "B"
#   selected_idx — 1-based rank of the highlighted term (matches table row)

if (!requireNamespace("ggplot2", quietly = TRUE))
  stop("ggplot2 is required for the dot plot. Install it with install.packages('ggplot2').")

# ── Colour palettes ─────────────────────────────────────────────────────────
.PAL_A <- c(               # Pipeline A: stability tiers
  high    = "#1e7d58",     # green  ≥ 0.85
  medium  = "#3560c8",     # blue   ≥ 0.70
  low     = "#9aaabb"      # gray   < 0.70 or NA
)
.PAL_B <- c(               # Pipeline B: direction
  UP      = "#3560c8",     # blue
  DOWN    = "#c04848",     # red
  unknown = "#9aaabb"      # gray
)

# ── Label truncation ─────────────────────────────────────────────────────────
.trunc <- function(x, n = 44L) {
  long <- nchar(x) > n
  x[long] <- paste0(substr(x[long], 1L, n - 1L), "…")
  x
}

# ── Main function ────────────────────────────────────────────────────────────
make_dot_plot <- function(df, pipeline, selected_idx = 1L) {

  # ── Guard: need at least one plottable row ────────────────────────────────
  df <- df[!is.na(df$q_value) & as.numeric(df$q_value) > 0, , drop = FALSE]
  if (nrow(df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = "No plottable data",
                          colour = "#9aaabb", fontface = "italic", size = 3.5) +
        ggplot2::theme_void()
    )
  }

  # ── Prepare data ──────────────────────────────────────────────────────────
  # Floor to avoid -log10(0) = Inf
  df$q_safe <- pmax(as.numeric(df$q_value), 1e-300)
  df$neg_log10_q <- -log10(df$q_safe)

  # Sort ascending q (most significant first), cap at 25 terms
  df <- df[order(df$q_safe), , drop = FALSE]
  df <- utils::head(df, 25L)

  # Truncated labels (deduplicate if needed)
  df$label <- .trunc(as.character(df$term_name), 40L)
  dup <- anyDuplicated(df$label) > 0
  if (dup) {
    is_dup <- duplicated(df$label) | duplicated(df$label, fromLast = TRUE)
    df$label[is_dup] <- paste0(df$label[is_dup], " (", df$collection[is_dup], ")")
  }

  # y-axis factor: level order puts most significant (row 1) at the TOP
  df$y_fct <- factor(df$label, levels = rev(df$label))

  # ── Colour groups ─────────────────────────────────────────────────────────
  if (pipeline == "A") {
    s <- as.numeric(df$stability)
    df$cgrp <- ifelse(is.na(s), "low",
               ifelse(s >= 0.85, "high",
               ifelse(s >= 0.70, "medium", "low")))
    pal    <- .PAL_A
    leg_title <- "Stability"
    leg_labels <- c(high   = "Stab. ≥0.85",
                    medium = "Stab. ≥0.70",
                    low    = "Stab. <0.70")
  } else {
    d <- toupper(trimws(as.character(df$direction)))
    d[is.na(df$direction) | !nzchar(d) | !d %in% c("UP", "DOWN")] <- "unknown"
    df$cgrp <- d
    pal    <- .PAL_B
    leg_title <- "Direction"
    leg_labels <- c(UP = "↑ Up", DOWN = "↓ Down", unknown = "—")
  }
  df$cgrp <- factor(df$cgrp, levels = names(pal))

  # ── Dot size from gene_count ───────────────────────────────────────────────
  gc <- as.numeric(df$gene_count)
  gc[is.na(gc)] <- stats::median(gc, na.rm = TRUE)
  if (all(is.na(gc))) gc[] <- 20
  gc <- pmax(gc, 1)
  gc_rng <- max(gc) - min(gc)
  df$dot_sz <- if (gc_rng < 1) rep(2.5, nrow(df))
               else 1.5 + 3.5 * (gc - min(gc)) / gc_rng

  # ── Selected term ─────────────────────────────────────────────────────────
  # selected_idx matches the rank column (1-based, most significant = 1)
  if ("rank" %in% names(df)) {
    df$is_sel <- df$rank == as.integer(selected_idx)
  } else {
    df$is_sel <- seq_len(nrow(df)) == as.integer(selected_idx)
  }
  df_unsel <- df[!df$is_sel, , drop = FALSE]
  df_sel   <- df[ df$is_sel, , drop = FALSE]

  # ── Build plot ────────────────────────────────────────────────────────────
  p <- ggplot2::ggplot(df, ggplot2::aes(x = neg_log10_q, y = y_fct)) +

    # Subtle horizontal guide lines
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = as.numeric(y_fct)),
      colour = "#e2eaf4", linewidth = 0.35
    ) +

    # Non-selected dots
    ggplot2::geom_point(
      data = df_unsel,
      ggplot2::aes(size = dot_sz, colour = cgrp),
      shape = 19L, alpha = 0.72, stroke = 0
    )

  # Selected dot — navy ring, then filled centre
  if (nrow(df_sel) > 0) {
    sel_col <- as.character(pal[as.character(df_sel$cgrp[1])])
    if (is.na(sel_col)) sel_col <- "#9aaabb"

    p <- p +
      ggplot2::geom_point(          # outer ring
        data = df_sel,
        ggplot2::aes(x = neg_log10_q, y = y_fct, size = dot_sz + 1.4),
        shape = 21L, fill = NA, colour = "#18243a", stroke = 1.6, alpha = 1
      ) +
      ggplot2::geom_point(          # filled inner dot
        data = df_sel,
        ggplot2::aes(x = neg_log10_q, y = y_fct, size = dot_sz),
        shape = 19L, colour = sel_col, alpha = 1, stroke = 0
      )
  }

  p <- p +
    ggplot2::scale_size_identity() +
    ggplot2::scale_colour_manual(
      values = pal,
      labels = leg_labels[names(pal)],
      name   = leg_title,
      drop   = TRUE,
      guide  = ggplot2::guide_legend(
        override.aes = list(size = 2.5, alpha = 1, shape = 19L)
      )
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, 0.10)),
      limits = c(0, NA)
    ) +
    ggplot2::labs(
      x = expression(-log[10](italic(q)-value)),
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y  = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_line(colour = "#e2eaf4", linewidth = 0.4),
      axis.text.y         = ggplot2::element_text(size = 10.5, colour = "#18202e",
                                                  hjust = 1,
                                                  margin = ggplot2::margin(r = 3)),
      axis.text.x         = ggplot2::element_text(size = 9.5, colour = "#546070"),
      axis.title.x        = ggplot2::element_text(size = 10.5, colour = "#546070",
                                                  margin = ggplot2::margin(t = 4)),
      legend.title        = ggplot2::element_text(size = 9.5, face = "bold",
                                                  colour = "#9aaabb"),
      legend.text         = ggplot2::element_text(size = 9.5, colour = "#546070"),
      legend.position     = "right",
      legend.margin       = ggplot2::margin(l = 2, r = 0),
      legend.key.size     = grid::unit(0.6, "lines"),
      plot.margin         = ggplot2::margin(t = 6, r = 4, b = 6, l = 2),
      plot.background     = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background    = ggplot2::element_rect(fill = "white", colour = NA)
    )

  p
}
