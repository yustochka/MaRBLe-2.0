# shiny/R/funnel_plot.R
# Funnel plot for the Pipeline Reduction chart tab.
# Requires ggplot2 (available in the project renv).

#' Build and return a ggplot2 horizontal bar funnel chart.
#'
#' @param funnel_df  data.frame with columns: stage, count, source
#'                   Rows should be ordered first→last pipeline stage.
#' @param pipeline   "A" or "B"
make_funnel_plot <- function(funnel_df, pipeline) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for the funnel plot. Install it with install.packages('ggplot2').")

  df <- funnel_df[!is.na(funnel_df$count) & funnel_df$count >= 0L, , drop = FALSE]

  if (nrow(df) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = .5, y = .5,
                          label = "No funnel data available",
                          colour = "#9aaabb", face = "italic", size = 3.5) +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
    )
  }

  # Factor: first stage at top, last at bottom
  df$stage_f <- factor(df$stage, levels = rev(df$stage))

  # Colour gradient: sky (first/top) → navy (last/bottom)
  n <- nrow(df)
  bar_cols <- grDevices::colorRampPalette(c("#2a4270", "#6899e2"))(n)
  names(bar_cols) <- levels(df$stage_f)   # level 1 (bottom/final) → navy; level n (top/first) → sky

  is_partial <- any(!df$source %in% c("stage file", "final file"))
  subtitle_txt <- if (is_partial)
    "Note: some stage counts derived from summary data — full per-step files not found."
  else NULL

  max_count <- max(df$count, na.rm = TRUE)

  ggplot2::ggplot(df, ggplot2::aes(x = count, y = stage_f, fill = stage_f)) +
    ggplot2::geom_col(width = 0.55, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = format(count, big.mark = ",")),
      hjust = -0.12, size = 3.8, colour = "#18243a", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = bar_cols) +
    ggplot2::scale_x_continuous(
      limits = c(0, max_count * 1.22),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title    = paste0("Pipeline ", pipeline, " — reduction funnel"),
      subtitle = subtitle_txt,
      x        = "Number of terms / pathways",
      y        = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y  = ggplot2::element_blank(),
      panel.grid.minor    = ggplot2::element_blank(),
      panel.grid.major.x  = ggplot2::element_line(colour = "#e2eaf4", linewidth = 0.4),
      axis.text.y         = ggplot2::element_text(size = 11, colour = "#18202e", hjust = 1),
      axis.text.x         = ggplot2::element_text(size = 9.5, colour = "#546070"),
      axis.title.x        = ggplot2::element_text(size = 10, colour = "#546070",
                                                   margin = ggplot2::margin(t = 4)),
      plot.title          = ggplot2::element_text(size = 11.5, face = "bold",
                                                   colour = "#18243a"),
      plot.subtitle       = ggplot2::element_text(size = 9, colour = "#9aaabb",
                                                   face = "italic"),
      plot.margin         = ggplot2::margin(t = 8, r = 14, b = 6, l = 4),
      plot.background     = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background    = ggplot2::element_rect(fill = "white", colour = NA)
    )
}
