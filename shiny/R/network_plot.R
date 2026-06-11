# shiny/R/network_plot.R
# Redundancy network visualization for the chart tab.
# Requires ggplot2. No igraph needed — uses a circle layout computed in R.

# Internal: empty / no-data state plot
.net_empty <- function(msg = "No data available.") {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    return(invisible(NULL))
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                      colour = "#9aaabb", fontface = "italic",
                      size = 3.2, hjust = 0.5, vjust = 0.5) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

#' Render the redundancy network chart.
#'
#' @param network_data  list returned by load_network_edges()
#' @param selected_id   term_id of currently highlighted row, or NULL
make_network_plot <- function(network_data, selected_id = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for the network plot.")

  # Empty states
  if (is.null(network_data) || !network_data$found || nrow(network_data$nodes) == 0L) {
    msg <- if (!is.null(network_data) && nzchar(network_data$message))
             network_data$message
           else
             "No shortlist terms available."
    return(.net_empty(msg))
  }

  nodes <- network_data$nodes
  edges <- network_data$edges
  n     <- nrow(nodes)

  # ── Layout ────────────────────────────────────────────────────────────────
  # Group: direction (UP/DOWN) for Pipeline B; collection (GO/WP) for Pipeline A
  has_dir <- any(!is.na(nodes$direction) & nzchar(as.character(nodes$direction)))
  sort_by <- if (has_dir) {
    ifelse(is.na(nodes$direction) | !nzchar(nodes$direction), "z",
           toupper(as.character(nodes$direction)))
  } else {
    as.character(nodes$collection)
  }
  nodes <- nodes[order(sort_by, na.last = TRUE), , drop = FALSE]
  rownames(nodes) <- NULL

  angle <- if (n == 1L) 0 else seq(0, 2 * pi * (1 - 1 / n), length.out = n)

  nodes$x <- cos(angle)
  nodes$y <- sin(angle)

  # ── Node size ────────────────────────────────────────────────────────────
  gc  <- as.numeric(nodes$gene_count)
  gc[is.na(gc)] <- 20
  gc  <- pmax(gc, 1)
  grng <- max(gc) - min(gc)
  nodes$sz <- if (grng < 1) rep(2.8, n) else 1.8 + 3.2 * (gc - min(gc)) / grng

  # ── Node colour ───────────────────────────────────────────────────────────
  nodes$col <- ifelse(
    grepl("GO", nodes$collection, ignore.case = FALSE), "#3560c8",
    ifelse(grepl("WP", nodes$collection, ignore.case = FALSE), "#1e7d58",
           "#9aaabb")
  )

  # ── Selected highlight ────────────────────────────────────────────────────
  sel_id       <- if (!is.null(selected_id) && length(selected_id) > 0L &&
                      !is.na(selected_id[1L]) && nzchar(selected_id[1L]))
                    as.character(selected_id[1L]) else NA_character_
  nodes$is_sel <- !is.na(sel_id) & as.character(nodes$id) == sel_id

  # ── Prepare edge coordinates ──────────────────────────────────────────────
  pos <- nodes[, c("id", "x", "y"), drop = FALSE]
  has_edges <- !is.null(edges) && nrow(edges) > 0L

  if (has_edges) {
    edges <- merge(edges, pos, by.x = "from", by.y = "id", all.x = TRUE)
    names(edges)[names(edges) == "x"] <- "x1"
    names(edges)[names(edges) == "y"] <- "y1"
    edges <- merge(edges, pos, by.x = "to", by.y = "id", all.x = TRUE)
    names(edges)[names(edges) == "x"] <- "x2"
    names(edges)[names(edges) == "y"] <- "y2"
    edges <- edges[!is.na(edges$x1) & !is.na(edges$x2), , drop = FALSE]
    has_edges <- nrow(edges) > 0L
  }

  # ── Caption note ──────────────────────────────────────────────────────────
  src <- if (!is.null(network_data$source)) network_data$source else ""
  note <- if (!has_edges && nrow(network_data$edges) == 0L) {
    "No redundancy/consensus edges found — showing shortlist terms."
  } else if (grepl("cross_collection", src, fixed = TRUE)) {
    "Edges show GO–WikiPathways consensus or overlap relationships."
  } else if (grepl("cluster", src, fixed = TRUE)) {
    "Edges approximated from shared cluster IDs."
  } else {
    "Edges show pairwise overlap (Jaccard) between shortlisted terms."
  }

  # ── Selected-node label text (truncated, used below) ─────────────────────
  sel_label <- if (!is.na(sel_id)) {
    raw <- as.character(nodes$label[as.character(nodes$id) == sel_id][1L])
    if (!is.na(raw) && nzchar(raw)) {
      if (nchar(raw) > 28L) paste0(substr(raw, 1L, 27L), "…") else raw
    } else NA_character_
  } else NA_character_

  # ── Build ggplot ──────────────────────────────────────────────────────────
  p <- ggplot2::ggplot()

  # Edges: overlap/cluster (light) drawn first, then consensus (blue) on top
  if (has_edges) {
    e_other <- edges[edges$edge_type != "consensus", , drop = FALSE]
    e_cons  <- edges[edges$edge_type == "consensus", , drop = FALSE]

    if (nrow(e_other) > 0L)
      p <- p + ggplot2::geom_segment(
        data = e_other,
        ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
        colour = "#ccd6e8", linewidth = 0.65, alpha = 0.85, lineend = "round"
      )

    if (nrow(e_cons) > 0L)
      p <- p + ggplot2::geom_segment(
        data = e_cons,
        ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2),
        colour = "#6899e2", linewidth = 1.4, alpha = 0.9, lineend = "round"
      )
  }

  # Non-selected nodes
  df_u <- nodes[!nodes$is_sel, , drop = FALSE]
  if (nrow(df_u) > 0L)
    p <- p + ggplot2::geom_point(
      data = df_u,
      ggplot2::aes(x = x, y = y, size = sz, colour = col),
      alpha = 0.85, stroke = 0, show.legend = FALSE
    )

  # Selected node: navy ring + refilled dot
  df_s <- nodes[nodes$is_sel, , drop = FALSE]
  if (nrow(df_s) > 0L)
    p <- p +
      ggplot2::geom_point(
        data = df_s,
        ggplot2::aes(x = x, y = y, size = sz + 1.8),
        shape = 21L, fill = NA, colour = "#18243a", stroke = 1.8
      ) +
      ggplot2::geom_point(
        data = df_s,
        ggplot2::aes(x = x, y = y, size = sz, colour = col),
        alpha = 1, stroke = 0, show.legend = FALSE
      )

  # Label: only the selected node, placed inside the circle at 52% of radius
  # so it can never clip at the coord boundary (always within ±0.52).
  if (!is.na(sel_id) && !is.na(sel_label) && any(nodes$is_sel)) {
    df_s_lbl <- nodes[nodes$is_sel, , drop = FALSE]
    df_s_lbl$lbl_x <- df_s_lbl$x * 0.52
    df_s_lbl$lbl_y <- df_s_lbl$y * 0.52
    df_s_lbl$disp  <- sel_label
    p <- p + ggplot2::geom_label(
      data = df_s_lbl,
      ggplot2::aes(x = lbl_x, y = lbl_y, label = disp),
      fill          = "white",
      colour        = "#18243a",
      size          = 3.4,
      fontface      = "bold",
      label.padding = ggplot2::unit(0.18, "lines"),
      label.size    = 0.25,
      hjust         = 0.5,
      vjust         = 0.5
    )
  }

  p +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_size_identity() +
    ggplot2::coord_fixed(ratio = 1,
                         xlim = c(-2.6, 2.6), ylim = c(-2.3, 2.3),
                         expand = FALSE) +
    ggplot2::labs(title = NULL, caption = note) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.caption    = ggplot2::element_text(size = 9.5, colour = "#9aaabb",
                                              face = "italic", hjust = 0.5,
                                              margin = ggplot2::margin(t = 4)),
      plot.margin     = ggplot2::margin(t = 8, r = 24, b = 6, l = 24),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}
