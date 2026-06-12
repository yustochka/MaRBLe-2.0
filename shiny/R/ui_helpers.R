# shiny/R/ui_helpers.R
# UI helper functions extracted from app.R to keep it manageable.
# Sourced by app.R at startup.
#
# Contents:
#   Advanced Settings modal (adv_param_row, adv_section, make_adv_modal)
#   Detail panel (make_detail_html, make_gene_preview)

# ── Advanced Settings modal helpers ─────────────────────────────────────────

adv_param_row <- function(label, value, desc = NULL) {
  tags$div(
    class = "adv-param-row",
    tags$div(
      class = "adv-param-left",
      tags$span(class = "adv-param-label", label),
      if (!is.null(desc)) tags$div(class = "adv-param-desc", desc)
    ),
    tags$span(class = "adv-param-value", value)
  )
}

adv_section <- function(title, ...) {
  tags$div(
    class = "adv-section",
    tags$div(class = "adv-section-title", title),
    ...
  )
}

# ── Pipeline-specific tab content ────────────────────────────────────────────

.adv_b_content <- function(b_top_n) {
  sel <- as.character(b_top_n)
  # Custom HTML replaces radioButtons() to avoid Bootstrap wrapper artifacts.
  # The shiny-input-radiogroup class is required for Shiny's radio binding.
  make_pill <- function(value, label) {
    tags$label(
      class = "top-n-pill",
      tags$input(type  = "radio",  name  = "adv_b_top_n", value = value,
                 checked = if (sel == value) NA else NULL),
      label
    )
  }

  tagList(
    tags$div(
      class = "adv-tab-content",
      adv_section(
        "Shortlist options",
        tags$div(
          class = "adv-radio-group",
          tags$div(class = "setup-label", "Top terms per direction"),
          tags$div(
            id    = "adv_b_top_n",
            class = "shiny-input-radiogroup top-n-group",
            make_pill("3",   "3"),
            make_pill("5",   "5"),
            make_pill("All", "All available")
          ),
          tags$div(class = "adv-param-desc",
                   "Filters the currently available precomputed shortlist. ",
                   tags$em("All available"), " shows all rows present in the loaded result file; ",
                   "it does not recover terms removed by the original pipeline cap.")
        )
      ),
      adv_section(
        "Pipeline B method (read-only)",
        adv_param_row("CAMERA / fgsea agreement", "Enabled",
          "Both tools must agree on direction and FDR for a term to advance."),
        adv_param_row("Semantic collapse",        "Enabled",
          "GO:BP terms collapsed using Wang semantic similarity (τ = 0.60)."),
        adv_param_row("Overlap clustering",       "Enabled",
          "Jaccard-based Louvain clustering removes redundant overlapping terms."),
        adv_param_row("Default cap",              "5 per direction",
          "Precomputed results are capped at 5 terms per direction (UP / DOWN).")
      )
    )
  )
}

.adv_a_content <- function(ds, a_source) {
  # Custom HTML replaces radioButtons() to avoid Bootstrap absolute-positioning
  # (Bootstrap 3 sets .radio label { padding-left: 20px } and
  # .radio input[type=radio] { position: absolute; margin-left: -20px }
  # which causes the circle to overlap the text).
  make_radio_row <- function(value, label) {
    tags$label(
      style = paste0(
        "display:flex;align-items:center;gap:8px;cursor:pointer;",
        "font-size:13px;font-weight:500;color:#18243a;margin:0;"
      ),
      tags$input(
        type    = "radio",
        name    = "adv_a_source",
        value   = value,
        checked = if (a_source == value) NA else NULL,
        style   = paste0(
          "position:static;margin:0;padding:0;",
          "width:14px;height:14px;flex-shrink:0;",
          "cursor:pointer;accent-color:#3560c8;"
        )
      ),
      label
    )
  }

  tagList(
    tags$div(
      class = "adv-tab-content",
      adv_section(
        "Result source",
        tags$div(
          class = "adv-radio-group",
          tags$div(class = "setup-label", "Load result from"),
          tags$div(
            id    = "adv_a_source",
            class = "shiny-input-radiogroup",
            style = "display:flex;flex-direction:column;gap:10px;margin-top:8px;margin-bottom:4px;",
            make_radio_row("Baseline",     "Baseline"),
            make_radio_row("Recommended",  "Recommended policy")
          ),
          tags$div(class = "adv-param-desc",
                   tags$strong("Baseline:"), " loads the standard precomputed result. ",
                   tags$strong("Recommended policy:"),
                   " loads results/analysis/policy_runs/", ds, "/FINAL.csv if available; ",
                   "falls back to Baseline if missing.")
        )
      ),
      adv_section(
        "Pipeline A method (read-only)",
        adv_param_row("ORA enrichment",        "Fisher's exact test",
          "Overrepresentation analysis for WikiPathways and GO:BP."),
        adv_param_row("GO semantic collapse",   "τ = 0.60",
          "GO:BP terms collapsed using Wang semantic similarity."),
        adv_param_row("Bootstrap stability",    "200 iterations, τ = 0.70",
          "Subsample-based stability filtering."),
        adv_param_row("Consensus Jaccard",      "0.15",
          "Minimum overlap for GO–WP cross-collection consensus pairs.")
      )
    )
  )
}

# ── Modal builder (pipeline-aware) ────────────────────────────────────────────
make_adv_modal <- function(ds = "", pipeline = "B",
                           settings = list(b_top_n = 5L, a_source = "Baseline")) {
  b_top_n  <- if (!is.null(settings$b_top_n))  settings$b_top_n  else 5L
  a_source <- if (!is.null(settings$a_source)) settings$a_source else "Baseline"

  modal_title <- tags$div(
    tags$div(class = "adv-modal-title", "Advanced Settings"),
    tags$div(class = "adv-modal-subtitle",
             "Settings control how precomputed results are displayed.")
  )

  body_content <- if (pipeline == "B") .adv_b_content(b_top_n)
                  else                  .adv_a_content(ds, a_source)

  modalDialog(
    title     = modal_title,
    size      = "l",
    easyClose = TRUE,
    footer    = tags$div(
      class = "adv-footer adv-footer-actions",
      modalButton("Cancel"),
      actionButton("adv_save_settings", "Save settings",
                   class = "marble-btn-primary adv-save-btn")
    ),

    tags$div(
      class = "adv-modal-body",

      tags$div(
        class = "adv-info-box",
        tags$span(class = "adv-info-icon", "ℹ"),
        tags$span(
          "These settings adjust the displayed shortlist using available precomputed results."
        )
      ),

      body_content
    )
  )
}

# ── Detail panel helpers ─────────────────────────────────────────────────────

make_detail_html <- function(row, pipeline, genes = NULL) {
  safe_chr <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x[1])) "—" else as.character(x[1])
  }
  fmt_q    <- function(q) if (!is.na(q)) formatC(as.numeric(q), format="e", digits=2) else "—"
  fmt_stab <- function(s) if (!is.na(s)) formatC(as.numeric(s), digits=2, format="f") else "—"

  tid        <- safe_chr(row$term_id)
  gene_src   <- if (grepl("^GO:", tid)) "org.Hs.eg.db" else "Pathway2Gene.csv"

  coll      <- safe_chr(row$collection)
  coll_lbl  <- if (grepl("GO", coll, fixed=TRUE)) "GO Biological Process" else "WikiPathways"
  dir_raw   <- row$direction
  dir_lbl   <- if (!is.na(dir_raw) && nzchar(as.character(dir_raw))) {
    if (toupper(dir_raw) == "UP") "↑ Upregulated" else "↓ Downregulated"
  } else NULL
  coll_line <- paste0(coll_lbl, if (!is.null(dir_lbl)) paste0(" · ", dir_lbl) else "")

  stat_row <- function(key, val, cls = "") {
    tags$div(
      class = "detail-stat-row",
      tags$span(class = "detail-stat-key", key),
      tags$span(class = paste("detail-stat-val", cls), val)
    )
  }

  if (pipeline == "B") {
    reason_val <- if ("final_reason" %in% names(row) &&
                      !is.na(row$final_reason) &&
                      nzchar(as.character(row$final_reason))) {
      as.character(row$final_reason)
    } else "Selected"

    tagList(
      tags$div(class = "detail-coll-line", coll_line),
      tags$div(
        class = "detail-stats-grid",
        stat_row("q-value",    fmt_q(row$q_value),       "qval-cell"),
        stat_row("Gene count", safe_chr(row$gene_count), ""),
        if (reason_val != "Selected") stat_row("Reason", reason_val, "") else NULL
      ),
      make_gene_preview(genes, gene_src)
    )
  } else {
    stab_raw <- row$stability
    stab_cls <- if (!is.na(stab_raw)) {
      s <- as.numeric(stab_raw)
      if (s >= 0.85) "stab-high" else if (s >= 0.75) "stab-mid" else "stab-low"
    } else "stab-na"
    cons_is_true <- !is.na(row$consensus) && as.logical(row$consensus)

    why <- paste0("Ranked #", safe_chr(row$rank), " by q-value. ",
                  if (!is.na(stab_raw))
                    paste0("Stability ", fmt_stab(stab_raw), " meets threshold. ")
                  else "",
                  if (cons_is_true) "Cross-collection consensus pair identified."
                  else "No consensus partner — kept on stability/rank alone.")

    tagList(
      tags$div(class = "detail-coll-line", coll_line),
      tags$div(
        class = "detail-stats-grid",
        stat_row("q-value",    fmt_q(row$q_value),       "qval-cell"),
        stat_row("Stability",  fmt_stab(stab_raw),       stab_cls),
        stat_row("Gene count", safe_chr(row$gene_count), ""),
        stat_row("Consensus",  if (cons_is_true) "Yes" else "No",
                               if (cons_is_true) "cons-yes" else "cons-no")
      ),
      tags$div(
        class = "why-kept-box",
        tags$div(class = "why-kept-label", "Why kept"),
        tags$div(class = "why-kept-text",  why)
      ),
      make_gene_preview(genes, gene_src)
    )
  }
}

# genes: character vector from gene_map(), or NULL
# source_label: string shown in the "Source:" note
make_gene_preview <- function(genes, source_label = "Pathway2Gene.csv") {
  max_show <- 12L

  if (is.null(genes)) {
    return(tagList(
      tags$div(class = "detail-gene-label", "Gene list preview"),
      tags$div(class = "gene-empty-note",
               "No gene mapping available for this term.")
    ))
  }

  if (length(genes) == 0L) {
    return(tagList(
      tags$div(class = "detail-gene-label", "Gene list preview"),
      tags$div(class = "gene-empty-note", "No genes annotated to this term.")
    ))
  }

  show <- head(genes, max_show)
  rest <- if (length(genes) > max_show)
            genes[(max_show + 1L):length(genes)] else character(0)

  tagList(
    tags$div(class = "detail-gene-label",
             paste0("Gene list preview (", length(genes), " mapped)")),
    tags$div(
      class = "gene-chips",
      lapply(show, function(g) tags$span(class = "gene-chip", g)),
      # Overflow genes live in a <details> toggle — click "+N more" to expand.
      if (length(rest) > 0L) tags$details(
        class = "gene-more",
        tags$summary(
          class = "gene-chip-more",
          tags$span(class = "more-collapsed", paste0("+", length(rest), " more")),
          tags$span(class = "more-expanded", "Show fewer")
        ),
        tags$div(
          class = "gene-chips gene-more-list",
          lapply(rest, function(g) tags$span(class = "gene-chip", g))
        )
      )
    ),
    tags$div(class = "gene-source-note", paste0("Source: ", source_label))
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  DATASET INFO PANEL  (setup screen)
# ══════════════════════════════════════════════════════════════════════════════

# ── Internal helpers ─────────────────────────────────────────────────────────

.ds_group_label <- function(dataset) {
  if (grepl("^dataset_", dataset, ignore.case = FALSE)) "Reference"
  else "Biochemistry application"
}

.fmt_final <- function(n, source) {
  if (is.null(n) || is.na(n)) return("—")
  n_int <- suppressWarnings(as.integer(n))
  if (is.na(n_int)) return("—")
  lbl <- switch(as.character(source),
    tier1                     = "Tier 1",
    tier2                     = "Tier 2",
    candidates                = "cand.",
    solo                      = "solo",
    fallback_topn_last_resort = "top-N",
    ""
  )
  if (nzchar(lbl)) paste0(n_int, " (", lbl, ")") else as.character(n_int)
}

.cls_badge <- function(cls) {
  if (is.null(cls) || is.na(cls) || !nzchar(cls))
    return(tags$span(class = "ds-cls-badge ds-cls-na", "—"))
  badge_cls <- switch(as.character(cls),
    HEALTHY      = "ds-cls-badge ds-cls-healthy",
    BORDERLINE   = "ds-cls-badge ds-cls-borderline",
    NO_OVERLAP   = "ds-cls-badge ds-cls-no-overlap",
    LOW_COVERAGE = "ds-cls-badge ds-cls-low-coverage",
    "ds-cls-badge ds-cls-na"
  )
  tags$span(class = badge_cls, cls)
}

.ds_metric <- function(label, value, sub = NULL) {
  tags$div(
    class = "ds-metric-item",
    tags$span(class = "ds-metric-label", label),
    tags$span(class = "ds-metric-value", value),
    if (!is.null(sub) && nzchar(sub))
      tags$span(class = "ds-metric-sub", sub)
    else NULL
  )
}

# ── Dataset interpretation sentence ─────────────────────────────────────────

#' Returns a short interpretation string, or NULL if none applicable.
ds_interpretation <- function(meta) {
  cls <- if (is.null(meta$classification) || is.na(meta$classification))
           NA_character_ else as.character(meta$classification)
  de  <- if (is.null(meta$de_fraction) || is.na(meta$de_fraction))
           NA_real_ else as.numeric(meta$de_fraction)
  ds  <- meta$dataset

  if (!is.na(de) && de < 0.01) {
    pct <- formatC(de * 100, digits = 2, format = "f")
    return(paste0(
      "Very low DE fraction (", pct, "%) — ",
      "Pipeline B is likely more informative than ORA-based Pipeline A."
    ))
  }

  # Classification-based fallback
  if (is.na(cls)) return(NULL)
  switch(cls,
    HEALTHY      = "Strong DE signal; Pipeline A results are expected to be reliable.",
    BORDERLINE   = paste0(
      "Moderate signal; consensus threshold was relaxed to achieve Tier 1 — ",
      "treat results with medium confidence."
    ),
    NO_OVERLAP   = paste0(
      "No GO–WP consensus achievable; Pipeline A reports Tier 2 only. ",
      "Pipeline B adds complementary directional information."
    ),
    LOW_COVERAGE = paste0(
      "Focused panel with low pathway coverage; ",
      "enrichment results should be interpreted with caution."
    ),
    NULL
  )
}

# ── Main panel builder ───────────────────────────────────────────────────────

#' Build the compact dataset info panel for the setup screen.
#' @param meta  list returned by load_dataset_meta()
make_dataset_info_panel <- function(meta) {
  if (is.null(meta)) return(NULL)

  ds  <- meta$dataset
  cls <- meta$classification
  grp <- .ds_group_label(ds)

  grp_class <- if (grepl("^dataset_", ds)) "ds-group-badge ds-grp-reference"
               else                         "ds-group-badge ds-grp-bio"

  # Check whether we have any metadata at all
  has_data <- any(!is.na(c(meta$n_universe, meta$de_fraction,
                            meta$classification, meta$n_final_a)))

  if (!has_data) {
    note <- "No pipeline results found for this dataset."
    return(tags$div(
      class = "ds-info-panel",
      tags$div(
        class = "ds-info-header",
        tags$span(class = "ds-info-name", ds),
        tags$span(class = grp_class, grp)
      ),
      tags$div(class = "ds-info-no-data", note)
    ))
  }

  # Format values
  safe_num <- function(x, fmt = "as.is", mult = 1) {
    if (is.null(x) || is.na(x)) return("—")
    v <- as.numeric(x) * mult
    switch(fmt,
      pct  = paste0(formatC(v, digits = 1, format = "f"), "%"),
      int  = format(as.integer(round(v)), big.mark = ","),
      dec2 = formatC(v, digits = 2, format = "f"),
      as.character(x)
    )
  }

  n_a   <- .fmt_final(meta$n_final_a, meta$final_source_a)
  n_b   <- .fmt_final(meta$n_final_b, meta$final_source_b)
  interp <- ds_interpretation(meta)

  has_policy <- !is.null(cls) && !is.na(cls) && nzchar(cls)

  tags$div(
    class = "ds-info-panel",

    # Header: name + group badge
    tags$div(
      class = "ds-info-header",
      tags$span(class = "ds-info-name", ds),
      tags$span(class = grp_class, grp)
    ),

    # Metrics row
    tags$div(
      class = "ds-info-metrics",
      .ds_metric("Universe",     safe_num(meta$n_universe, "int"), "genes"),
      .ds_metric("DE",           safe_num(meta$de_fraction, "pct", 100)),
      tags$div(
        class = "ds-metric-item",
        tags$span(class = "ds-metric-label", "Classification"),
        .cls_badge(cls)
      ),
      .ds_metric("Pipeline A", n_a),
      .ds_metric("Pipeline B", n_b)
    ),

    # Thresholds row (only when policy data available)
    if (has_policy) tags$div(
      class = "ds-info-thresholds",
      tags$span(class = "ds-thr-label", "Recommended"),
      tags$span(class = "ds-thr-item",
                tags$span(class = "ds-thr-key", "TAU"),
                tags$span(class = "ds-thr-val",  safe_num(meta$tau, "dec2"))),
      tags$span(class = "ds-thr-item",
                tags$span(class = "ds-thr-key", "CONS_J"),
                tags$span(class = "ds-thr-val",  safe_num(meta$cons_j, "dec2"))),
      tags$span(class = "ds-thr-item",
                tags$span(class = "ds-thr-key", "Confidence"),
                tags$span(class = "ds-thr-val",
                          if (is.null(meta$confidence) || is.na(meta$confidence))
                            "—" else as.character(meta$confidence)))
    ) else NULL,

    # Interpretation
    if (!is.null(interp))
      tags$div(class = "ds-info-interp", interp)
    else NULL
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  PIPELINE A vs B COMPARISON VIEW
# ══════════════════════════════════════════════════════════════════════════════

# Internal: outcome badge <span>
.cmp_outcome_badge <- function(outcome_raw) {
  lbl <- switch(as.character(outcome_raw),
    "Complementary" = "Complementary",
    "B > A"         = "B outperforms A",
    "A > B"         = "A outperforms B",
    as.character(outcome_raw)
  )
  cls <- switch(as.character(outcome_raw),
    "Complementary" = "cmp-outcome-badge cmp-outcome-complementary",
    "B > A"         = "cmp-outcome-badge cmp-outcome-b-gt-a",
    "A > B"         = "cmp-outcome-badge cmp-outcome-a-gt-b",
    "cmp-outcome-badge cmp-outcome-complementary"
  )
  tags$span(class = cls, lbl)
}

#' Render the A ↔ B comparison tab panel for a given dataset.
#'
#' @param ds           Dataset name string.
#' @param comp_summary Named list from load_comparison_summary(), or NULL.
#' @param overlap_tbl  Data frame from load_overlap_table().
make_compare_tab_content <- function(ds, comp_summary, overlap_tbl) {

  # ── No data: friendly message ───────────────────────────────────────────────
  if (is.null(comp_summary)) {
    return(tags$div(
      class = "compare-view",
      tags$div(
        class = "cmp-not-available",
        "No precomputed A↔B comparison is available for this dataset."
      )
    ))
  }

  # ── Helpers for safe extraction ─────────────────────────────────────────────
  s2c <- function(x) {
    if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) return("")
    as.character(x)
  }
  s2i <- function(x) {
    v <- suppressWarnings(as.integer(x))
    if (is.na(v)) "—" else as.character(v)
  }
  s2n <- function(x, digits = 1) {
    v <- suppressWarnings(as.numeric(x))
    if (is.na(v)) "—" else formatC(v, digits = digits, format = "f")
  }

  # ── Summary card ────────────────────────────────────────────────────────────
  s <- comp_summary

  n_ov    <- s2i(s$n_exact_overlap)
  agr     <- s2c(s$agreement_level)
  de_pct  <- {v <- suppressWarnings(as.numeric(s$de_pct)); if (is.na(v)) "—" else paste0(v, "%")}
  n_de    <- s2i(s$n_de)
  obs_txt <- s2c(s$key_observation)

  a_tier <- s2c(s$a_tier)
  a_cls  <- s2c(s$a_class)
  a_sub  <- paste0(if (nzchar(a_tier)) paste0(a_tier, " · ") else "", a_cls)
  b_sub  <- {b <- s2c(s$b_agreement); if (nzchar(b)) paste0("agreement: ", b) else ""}

  header_right <- paste0(
    "overlap: ", n_ov,
    if (nzchar(agr)) paste0("  ·  ", agr) else ""
  )

  summ_card <- tags$div(
    class = "cmp-summary-card",
    tags$div(
      class = "cmp-summary-header",
      .cmp_outcome_badge(s$outcome),
      tags$span(class = "cmp-agreement-text", header_right)
    ),
    tags$div(
      class = "cmp-metrics-row",
      tags$div(class = "cmp-metric",
               tags$div(class = "cmp-metric-label", "Pipeline A"),
               tags$div(class = "cmp-metric-value", s2i(s$n_pipeline_a)),
               tags$div(class = "cmp-metric-sub", a_sub)),
      tags$div(class = "cmp-metric-div"),
      tags$div(class = "cmp-metric",
               tags$div(class = "cmp-metric-label", "Pipeline B"),
               tags$div(class = "cmp-metric-value", s2i(s$n_pipeline_b)),
               tags$div(class = "cmp-metric-sub", b_sub)),
      tags$div(class = "cmp-metric-div"),
      tags$div(class = "cmp-metric",
               tags$div(class = "cmp-metric-label", "DE fraction"),
               tags$div(class = "cmp-metric-value", de_pct),
               tags$div(class = "cmp-metric-sub", paste0("n = ", n_de, " genes")))
    ),
    if (nzchar(obs_txt)) tags$div(class = "cmp-observation", obs_txt) else NULL
  )

  # ── Overlap table ────────────────────────────────────────────────────────────
  overlap_section <- if (!is.null(overlap_tbl) && nrow(overlap_tbl) > 0L) {
    rows <- lapply(seq_len(nrow(overlap_tbl)), function(i) {
      r <- overlap_tbl[i, ]
      tid <- s2c(r$ID)
      dsc <- s2c(r$Description)

      stab_span <- {
        v <- suppressWarnings(as.numeric(r$stability))
        if (!is.na(v)) {
          cls <- if (v >= 0.90) "stab-high" else if (v >= 0.70) "stab-mid" else "stab-low"
          tags$span(class = cls, formatC(v, digits = 2, format = "f"))
        } else tags$span("—")
      }

      dir_raw  <- toupper(s2c(r$direction))
      dir_span <- if (dir_raw == "UP")
        tags$span(class = "dir-badge dir-badge-up",   HTML("&#8593; UP"))
      else if (dir_raw == "DOWN")
        tags$span(class = "dir-badge dir-badge-down", HTML("&#8595; DOWN"))
      else tags$span(dir_raw)

      nes_txt <- {
        v <- suppressWarnings(as.numeric(r$absNES))
        if (!is.na(v)) formatC(v, digits = 2, format = "f") else "—"
      }

      tags$tr(
        tags$td(class = "cmp-id-cell", tid),
        tags$td(dsc),
        tags$td(stab_span),
        tags$td(dir_span),
        tags$td(class = "qval-cell", nes_txt)
      )
    })

    tagList(
      tags$div(class = "cmp-section-label", "Shared pathway IDs"),
      tags$div(
        class = "cmp-overlap-table",
        tags$table(
          tags$thead(tags$tr(
            tags$th("ID"), tags$th("Description"),
            tags$th("A stability"), tags$th("B direction"), tags$th("B |NES|")
          )),
          tags$tbody(rows)
        )
      )
    )
  } else {
    tagList(
      tags$div(class = "cmp-section-label", "Shared pathway IDs"),
      tags$div(class = "cmp-no-data", "No exact pathway ID overlap for this dataset.")
    )
  }

  tags$div(
    class = "compare-view",
    summ_card,
    overlap_section
  )
}
