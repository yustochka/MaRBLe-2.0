library(shiny)

# ── Source data loader ──────────────────────────────────────────────────────
dl_path <- if (file.exists("R/data_loader.R"))        "R/data_loader.R" else
           if (file.exists("shiny/R/data_loader.R")) "shiny/R/data_loader.R" else
           stop("Cannot find data_loader.R. Run: shiny::runApp('shiny')")
source(dl_path, local = TRUE)

# ── Source dot plot helper ───────────────────────────────────────────────────
dp_path <- if (file.exists("R/dot_plot.R"))        "R/dot_plot.R" else
           if (file.exists("shiny/R/dot_plot.R")) "shiny/R/dot_plot.R" else
           stop("Cannot find dot_plot.R. Run: shiny::runApp('shiny')")
source(dp_path, local = TRUE)

# ── Source funnel plot helper ────────────────────────────────────────────────
fp_path <- if (file.exists("R/funnel_plot.R"))        "R/funnel_plot.R" else
           if (file.exists("shiny/R/funnel_plot.R")) "shiny/R/funnel_plot.R" else
           stop("Cannot find funnel_plot.R. Run: shiny::runApp('shiny')")
source(fp_path, local = TRUE)

# ── Source network plot helper ───────────────────────────────────────────────
np_path <- if (file.exists("R/network_plot.R"))        "R/network_plot.R" else
           if (file.exists("shiny/R/network_plot.R")) "shiny/R/network_plot.R" else
           stop("Cannot find network_plot.R. Run: shiny::runApp('shiny')")
source(np_path, local = TRUE)

# ── Source UI helpers ────────────────────────────────────────────────────────
uh_path <- if (file.exists("R/ui_helpers.R"))        "R/ui_helpers.R" else
           if (file.exists("shiny/R/ui_helpers.R")) "shiny/R/ui_helpers.R" else
           stop("Cannot find ui_helpers.R. Run: shiny::runApp('shiny')")
source(uh_path, local = TRUE)

# ── Bootstrap values ────────────────────────────────────────────────────────
datasets     <- detect_datasets()
n_datasets   <- length(datasets)

# Datasets that have at least one pipeline run recorded in dataset_summary.csv
.ds_summ             <- load_dataset_summary()
datasets_with_results <- if (nrow(.ds_summ) > 0 && "dataset" %in% names(.ds_summ)) {
  unique(.ds_summ$dataset)
} else {
  datasets   # fallback: treat all as available
}

# Default selection: first dataset that actually has results; else first overall
.avail_ds    <- datasets[datasets %in% datasets_with_results]
init_dataset <- if (length(.avail_ds) > 0) .avail_ds[1] else
                if (length(datasets)  > 0) datasets[1]  else "dataset_0"

# Pre-load policy CSV for fast panel rendering (tiny file; shared across users)
.policy_cache <- load_threshold_policy()

# ── NULL-safe helper ────────────────────────────────────────────────────────
`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (is.na(a[1])  || !nzchar(a[1])) return(b)
  a[1]
}

# ══════════════════════════════════════════════════════════════════════════════
#  SHARED UI: NAVIGATION BAR
# ══════════════════════════════════════════════════════════════════════════════
make_nav <- function(active_step   = 1,
                     context_sub   = NULL,    # NULL → default subtitle
                     right_actions = NULL) {
  make_step <- function(i, label) {
    state <- if (i < active_step) "done" else if (i == active_step) "active" else "pending"
    tags$div(
      class = "step-item",
      tags$div(class = paste0("step-circle step-", state),
               if (state == "done") "✓" else as.character(i)),
      tags$span(class = paste0("step-label step-label-", state), label)
    )
  }
  tags$nav(
    class = "marble-nav",
    tags$div(
      class = "marble-nav-inner",
      tags$div(
        class = "marble-brand",
        tags$img(
          src   = "logo-header-cropped.png?v=2",
          alt   = "logo",
          style = "width:28px;height:28px;display:block;flex-shrink:0;object-fit:contain;"
        ),
        tags$span(class = "marble-title", "Pathway Analysis Dashboard"),
        if (!is.null(context_sub)) tagList(
          tags$span(class = "marble-sep"),
          tags$span(class = "marble-subtitle", context_sub)
        )
      ),
      tags$div(
        class = "marble-stepper",
        make_step(1, "Setup"),
        tags$div(class = "step-line"),
        make_step(2, "Results")
      ),
      if (!is.null(right_actions))
        tags$div(class = "nav-right-actions", right_actions)
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  SETUP SCREEN COMPONENTS
# ══════════════════════════════════════════════════════════════════════════════
make_dataset_chips <- function(datasets, selected = datasets[1],
                               available = datasets) {
  make_chip <- function(d) {
    has_data   <- d %in% available
    chip_class <- if (!has_data) "dataset-chip chip-unavailable"
                  else if (identical(d, selected)) "dataset-chip chip-active"
                  else "dataset-chip"
    if (has_data) {
      tags$button(class = chip_class, `data-ds` = d,
                  onclick = sprintf("selectDataset(this, '%s')", d), d)
    } else {
      tags$button(class = chip_class, `data-ds` = d, disabled = NA,
                  title = "No pipeline results available for this dataset",
                  d, tags$span(class = "chip-no-data-badge", "no results"))
    }
  }

  make_group <- function(ds_group, label) {
    if (length(ds_group) == 0L) return(NULL)
    tags$div(
      class = "chip-group",
      tags$div(class = "chip-group-label", label),
      tags$div(class = "dataset-chips", lapply(ds_group, make_chip))
    )
  }

  ref_ds <- datasets[grepl("^dataset_", datasets)]
  bio_ds <- datasets[!grepl("^dataset_", datasets)]

  tags$div(
    tags$div(class = "setup-label", "Dataset"),
    make_group(ref_ds, "Reference"),
    make_group(bio_ds, "Biochemistry application"),
    tags$div(
      class = "dataset-note",
      tags$span(class = "note-dot"),
      tags$span(
        class = "note-text",
        as.character(n_datasets),
        if (n_datasets == 1) " dataset" else " datasets",
        " auto-detected from ",
        tags$code(class = "note-code", "data/")
      )
    )
  )
}

pipeline_defs <- list(
  list(key = "A", name = "Pipeline A",
       desc = "ORA + redundancy reduction + stability + consensus"),
  list(key = "B", name = "Pipeline B",
       desc = "CAMERA / fgsea + redundancy reduction + top-N shortlist")
)

make_pipeline_cards <- function(selected = "B") {
  tags$div(
    tags$div(class = "setup-label", "Pipeline"),
    tags$div(
      class = "pipeline-cards",
      lapply(pipeline_defs, function(p) {
        is_sel <- identical(p$key, selected)
        radio_input <- if (is_sel) {
          tags$input(type = "radio", name = "pipeline_sel",
                     value = p$key, checked = NA)
        } else {
          tags$input(type = "radio", name = "pipeline_sel", value = p$key)
        }
        tags$label(
          class           = if (is_sel) "pipeline-card card-active" else "pipeline-card",
          `data-pipeline` = p$key,
          radio_input,
          tags$div(class = "pipeline-radio-dot"),
          tags$div(
            class = "pipeline-card-text",
            tags$div(class = "pipeline-card-name", p$name),
            tags$div(class = "pipeline-card-desc", p$desc)
          )
        )
      })
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  ADVANCED SETTINGS MODAL — helpers in shiny/R/ui_helpers.R
# ══════════════════════════════════════════════════════════════════════════════

# (adv_param_row, adv_section, make_adv_modal are defined in shiny/R/ui_helpers.R)

# (adv_param_row, adv_section, make_adv_modal fully defined in shiny/R/ui_helpers.R)


# ── Setup screen ──────────────────────────────────────────────────────────────
make_setup_screen <- function(ds = init_dataset, pipeline = "B", preset = "Default") {
  tagList(
    make_nav(active_step = 1),
    tags$div(
      class = "marble-screen-bg",
      tags$div(
        class = "marble-setup-card",
        tags$div(
          class = "setup-card-header",
          tags$div(class = "setup-card-dot"),
          tags$span("Configure Analysis")
        ),
        tags$div(
          class = "setup-card-body",
          make_dataset_chips(datasets, selected = ds, available = datasets_with_results),
          uiOutput("dataset_info_panel"),
          tags$hr(class = "setup-divider"),
          make_pipeline_cards(selected = pipeline),
          tags$hr(class = "setup-divider"),
          actionButton(
            "open_advanced",
            label = tagList(
              tags$span(class = "adv-row-icon",  "⚙"),
              tags$span(class = "adv-row-label", "Advanced settings"),
              tags$span(class = "adv-row-arrow", "→")
            ),
            class = "advanced-row-btn"
          ),
          tags$div(
            class = "load-button-wrap",
            actionButton(
              inputId = "load_results",
              label   = HTML('<span class="btn-run-icon">&#9654;</span>&nbsp;Load shortlist'),
              class   = "marble-btn-primary"
            ),
            tags$div(class = "load-hint",
                     "Loads cached results automatically if already computed.")
          )
        )
      )
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  RESULTS SCREEN COMPONENTS
# ══════════════════════════════════════════════════════════════════════════════

# ── Summary card (shared) ────────────────────────────────────────────────────
make_summary_card <- function(label, value, sub = NULL,
                              border_class = "card-border-blue") {
  disp <- if (is.null(value) || (length(value) == 1 && is.na(value))) "—"
          else if (inherits(value, "html")) value
          else as.character(value)
  tags$div(
    class = paste("summary-card", border_class),
    tags$div(
      class = "summary-card-inner",
      tags$div(class = "summary-card-label", label),
      tags$div(class = "summary-card-value", disp),
      if (!is.null(sub)) tags$div(class = "summary-card-sub", sub)
    )
  )
}

# ── Pipeline B top-N filter ──────────────────────────────────────────────────
# n: integer (3, 5) or "All" — applied after loading the precomputed shortlist.
filter_pipelineB_top_n <- function(df, n) {
  if (is.null(n) || identical(n, "All") || identical(n, "All available")) return(df)
  n_int <- as.integer(n)
  if (is.na(n_int) || n_int <= 0L) return(df)

  has_dir <- "direction" %in% names(df) &&
             any(!is.na(df$direction) & nzchar(as.character(df$direction)))

  if (has_dir) {
    dir_key <- toupper(trimws(as.character(df$direction)))
    dir_key[is.na(df$direction) | !nzchar(as.character(df$direction))] <- "UNKNOWN"
    groups <- split(seq_len(nrow(df)), dir_key)
    keep_idx <- unlist(lapply(groups, function(idx) {
      sub_ord <- idx[order(df$q_value[idx], na.last = TRUE)]
      sub_ord[seq_len(min(n_int, length(sub_ord)))]
    }))
    filtered <- df[sort(keep_idx), , drop = FALSE]
  } else {
    ord      <- order(df$q_value, na.last = TRUE)
    filtered <- df[ord, , drop = FALSE][seq_len(min(n_int, nrow(df))), , drop = FALSE]
  }

  filtered <- filtered[order(filtered$q_value, na.last = TRUE), , drop = FALSE]
  filtered$rank <- seq_len(nrow(filtered))
  rownames(filtered) <- NULL
  filtered
}

# ── Build display data frame for DT ─────────────────────────────────────────
# Returns a pipeline-specific 6-column data.frame ready for DT::datatable().
# Original data (r$data) is NOT modified; this is display-only.
build_dt_display <- function(df, pipeline) {
  if (nrow(df) == 0) return(data.frame())

  q_vals <- ifelse(
    is.na(df$q_value), "",
    formatC(as.numeric(df$q_value), format = "e", digits = 2)
  )

  gene_vals <- ifelse(
    is.na(df$gene_count), "",
    as.character(as.integer(df$gene_count))
  )

  if (pipeline == "B") {
    dir_vals <- ifelse(
      is.na(df$direction) | !nzchar(as.character(df$direction)),
      "",
      ifelse(toupper(as.character(df$direction)) == "UP",
        '<span class="dir-badge dir-badge-up">&#8593; UP</span>',
        '<span class="dir-badge dir-badge-down">&#8595; DOWN</span>')
    )

    data.frame(
      "Term"       = as.character(df$term_name),
      "Collection" = as.character(df$collection),
      "Direction"  = dir_vals,
      "q-value"    = q_vals,
      "Genes"      = gene_vals,
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
  } else {
    stab_vals <- ifelse(
      is.na(df$stability), "",
      formatC(as.numeric(df$stability), digits = 2, format = "f")
    )
    cons_vals <- ifelse(
      !is.na(df$consensus) & as.logical(df$consensus), "✓", "–"
    )

    data.frame(
      "Term"       = as.character(df$term_name),
      "Collection" = as.character(df$collection),
      "q-value"    = q_vals,
      "Stability"  = stab_vals,
      "Genes"      = gene_vals,
      "Consensus"  = cons_vals,
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
  }
}

# ── Chart card (with tab-switching via JS → input$chart_tab) ─────────────────
.CHART_TABS <- list(
  list(id = "dot_plot",    label = "Dot plot"),
  list(id = "funnel",      label = "Selection funnel"),
  list(id = "network",     label = "Overlap network"),
  list(id = "compare_ab",  label = "A ↔ B")
)

make_chart_card <- function() {
  tags$div(
    class = "results-card chart-card",
    tags$div(
      class = "chart-tabs",
      lapply(seq_along(.CHART_TABS), function(i) {
        t <- .CHART_TABS[[i]]
        tags$button(
          class   = if (i == 1L) "chart-tab chart-tab-active" else "chart-tab",
          onclick = sprintf("switchChartTab(this, '%s')", t$id),
          t$label
        )
      })
    ),
    # Single output; server switches content based on chart_tab reactive
    tags$div(class = "chart-area", uiOutput("chart_area_content"))
  )
}

# (make_detail_html and make_gene_preview are defined in shiny/R/ui_helpers.R)

# ── Parameters strip ─────────────────────────────────────────────────────────
make_params_strip <- function(cards, result) {
  preset <- result$effective_preset %||% "Default"
  ds     <- cards$dataset %||% ""
  disp_preset <- sub("\\s*\\(.*\\)$", "", preset)   # strip parenthetical note

  # Default thresholds; override from policy CSV for Recommended
  tau <- 0.70; cons_j <- 0.15
  if (grepl("Recommended", preset, fixed = TRUE) && nzchar(ds)) {
    pol <- tryCatch(load_threshold_policy(), error = function(e) data.frame())
    if (nrow(pol) > 0 && "dataset" %in% names(pol)) {
      r <- pol[pol$dataset == ds, , drop = FALSE]
      if (nrow(r) > 0) {
        if (!is.na(r$recommended_TAU_STABILITY[1]))  tau    <- r$recommended_TAU_STABILITY[1]
        if (!is.na(r$recommended_CONS_JACCARD_MIN[1])) cons_j <- r$recommended_CONS_JACCARD_MIN[1]
      }
    }
  }

  params <- list(
    FDR       = "0.05",
    Stability = formatC(tau,    digits = 2, format = "f"),
    Jaccard   = formatC(cons_j, digits = 2, format = "f"),
    Bootstrap = "200",
    Preset    = disp_preset
  )

  tags$div(
    class = "params-strip",
    tags$span(class = "params-label", "Parameters"),
    lapply(names(params), function(k) {
      tags$span(
        class = "params-item",
        tags$span(class = "params-key", k, " "),
        tags$span(class = "params-val", params[[k]])
      )
    })
  )
}

# ── Assembled results screen ──────────────────────────────────────────────────
make_results_screen <- function(cards, result) {
  df       <- result$data
  pipeline <- cards$pipeline %||% "B"
  preset   <- result$effective_preset %||% "Default"
  disp_preset <- sub("\\s*\\(.*\\)$", "", preset)
  is_fallback <- grepl("not available", preset, fixed = TRUE)

  policy_border <- switch(
    as.character(cards$policy_class %||% ""),
    "HEALTHY"      = "card-border-green",
    "BORDERLINE"   = "card-border-sky",
    "NO_OVERLAP"   = "card-border-sky",
    "LOW_COVERAGE" = "card-border-slate",
    "card-border-slate"
  )

  context_sub <- paste0(cards$dataset %||% "", " · Pipeline ", pipeline,
                        " · ", disp_preset)

  right_acts <- tagList(
    actionButton("back_to_setup", "← Setup",
                 class = "marble-btn-ghost nav-btn-sm"),
    downloadButton("download_results",
                   label = "⬇ Download shortlist",
                   class = "nav-download-btn")
  )

  tagList(
    make_nav(active_step = 2, context_sub = context_sub, right_actions = right_acts),

    tags$div(
      class = "res-page",

      # ── Summary cards ────────────────────────────────────────────────
      tags$div(
        class = "summary-cards-row",
        make_summary_card("Dataset",     cards$dataset,       border_class = "card-border-slate"),
        make_summary_card("Pipeline",    paste("Pipeline", pipeline), border_class = "card-border-blue"),
        make_summary_card("Final terms", nrow(df),            border_class = "card-border-blue"),
        make_summary_card("Universe",    cards$universe_size, sub = "genes", border_class = "card-border-sky"),
        if (pipeline == "B") {
          n_up   <- sum(!is.na(df$direction) & toupper(df$direction) == "UP",   na.rm = TRUE)
          n_down <- sum(!is.na(df$direction) & toupper(df$direction) == "DOWN", na.rm = TRUE)
          dir_disp <- HTML(paste0(
            '<span class="dir-up-count">&#8593;&nbsp;', n_up, '</span>',
            '<span class="dir-sep"> UP</span>',
            '&ensp;',
            '<span class="dir-down-count">&#8595;&nbsp;', n_down, '</span>',
            '<span class="dir-sep"> DOWN</span>'
          ))
          make_summary_card("Direction", dir_disp, border_class = "card-border-sky")
        } else {
          make_summary_card("Policy", cards$policy_class, border_class = policy_border)
        }
      ),

      # ── Preset fallback notice ────────────────────────────────────────
      if (is_fallback) tags$div(
        class = "preset-fallback-note",
        tags$span("⚠️"),
        tags$span(result$message)
      ),

      # ── Main dashboard grid ───────────────────────────────────────────
      tags$div(
        class = "res-grid",

        # LEFT: shortlist table
        tags$div(
          class = "res-left results-card",
          tags$div(
            class = "results-card-header",
            tags$span(class = "results-card-dot"),
            tags$span(class = "results-card-header-title", "Final Shortlist"),
            tags$span(class = "shortlist-count-badge",
                      paste0(nrow(df), " term", if (nrow(df) != 1) "s"))
          ),
          if (pipeline == "B") tags$div(
            class = "dir-filter-bar",
            tags$span(class = "dir-filter-label", "Direction"),
            tags$div(
              class = "dir-filter-btns",
              tags$button(class = "dir-filter-btn dir-filter-active",
                          onclick = "setDirFilter(this, 'All')", "All"),
              tags$button(class = "dir-filter-btn dir-up-btn",
                          onclick = "setDirFilter(this, 'UP')", HTML("&#8593;&nbsp;UP")),
              tags$button(class = "dir-filter-btn dir-down-btn",
                          onclick = "setDirFilter(this, 'DOWN')", HTML("&#8595;&nbsp;DOWN"))
            )
          ),
          tags$div(
            class = "shortlist-table-wrap",
            DT::DTOutput("shortlist_table", width = "100%")
          )
        ),

        # RIGHT: chart + detail (server-rendered) + params
        tags$div(
          class = "res-right",
          make_chart_card(),
          uiOutput("detail_panel"),
          make_params_strip(cards, result),
          if (!is.null(result$filter_note) &&
              length(result$filter_note) > 0 &&
              nzchar(result$filter_note)) {
            tags$div(class = "filter-note", result$filter_note)
          }
        )
      )
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  NOT-FOUND SCREEN
# ══════════════════════════════════════════════════════════════════════════════
make_not_found_screen <- function(result, ds, pipeline, preset) {
  tagList(
    make_nav(active_step = 1),
    tags$div(
      class = "marble-screen-bg",
      tags$div(
        class = "not-found-card",
        tags$div(
          class = "not-found-header",
          tags$span(class = "not-found-icon", "⚠"),
          tags$span(class = "not-found-title", "Results not found")
        ),
        tags$div(
          class = "not-found-body",
          tags$div(class = "not-found-message", result$message),
          tags$div(
            class = "not-found-meta",
            tags$div(tags$strong("Dataset:  "), ds),
            tags$div(tags$strong("Pipeline: "), paste("Pipeline", pipeline)),
            tags$div(tags$strong("Preset:   "), preset)
          ),
          actionButton("back_to_setup", "←  Back to Setup",
                       class = "marble-btn-ghost")
        )
      )
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  JAVASCRIPT
# ══════════════════════════════════════════════════════════════════════════════
setup_js <- tags$script(HTML(paste0("
(function () {

  // ── Setup screen: dataset chip ──────────────────────────────────────────
  window.selectDataset = function (btn, dataset) {
    document.querySelectorAll('.dataset-chip').forEach(function (b) {
      b.classList.remove('chip-active');
    });
    btn.classList.add('chip-active');
    Shiny.setInputValue('dataset', dataset, { priority: 'event' });
  };

  // ── Setup screen: pipeline card ─────────────────────────────────────────
  document.addEventListener('change', function (e) {
    if (e.target && e.target.name === 'pipeline_sel') {
      document.querySelectorAll('.pipeline-card').forEach(function (card) {
        card.classList.remove('card-active');
      });
      e.target.closest('.pipeline-card').classList.add('card-active');
      Shiny.setInputValue('pipeline', e.target.value, { priority: 'event' });
    }
  });

  // ── Chart tab switching ─────────────────────────────────────────────────
  window.switchChartTab = function (btn, tab) {
    document.querySelectorAll('.chart-tab').forEach(function (b) {
      b.classList.remove('chart-tab-active');
    });
    btn.classList.add('chart-tab-active');
    Shiny.setInputValue('chart_tab', tab, { priority: 'event' });
  };

  // ── Direction filter (Pipeline B) ───────────────────────────────────────
  window.setDirFilter = function (btn, val) {
    document.querySelectorAll('.dir-filter-btn').forEach(function (b) {
      b.classList.remove('dir-filter-active');
    });
    btn.classList.add('dir-filter-active');
    Shiny.setInputValue('dir_filter', val, { priority: 'event' });
  };

  // ── Initial Shiny values (once connected) ───────────────────────────────
  document.addEventListener('shiny:connected', function () {
    Shiny.setInputValue('dataset',  '", init_dataset, "');
    Shiny.setInputValue('pipeline', 'B');
    Shiny.setInputValue('preset',   'Default');
  });

}());
")))

# ══════════════════════════════════════════════════════════════════════════════
#  UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- bootstrapPage(
  title = "Pathway Analysis Dashboard",
  tags$head(
    tags$meta(charset = "UTF-8"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1.0"),
    tags$link(rel = "stylesheet", href = "marble.css?v=6")
  ),
  tags$div(class = "marble-page", uiOutput("page_content")),
  setup_js
)

# ══════════════════════════════════════════════════════════════════════════════
#  SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  app_state    <- reactiveVal("setup")
  loaded_data  <- reactiveVal(NULL)
  loaded_cards <- reactiveVal(NULL)
  raw_result   <- reactiveVal(NULL)   # unfiltered Pipeline B result for re-filter on settings change
  selected_idx <- reactiveVal(1L)
  gene_map     <- reactiveVal(list())
  chart_tab    <- reactiveVal("dot_plot")
  funnel_data  <- reactiveVal(NULL)
  network_data <- reactiveVal(NULL)
  adv_settings <- reactiveVal(list(b_top_n = 5L, a_source = "Baseline"))
  dir_filter   <- reactiveVal("All")

  # filtered_data: applies direction filter on top of loaded_data (Pipeline B only).
  # All other display outputs consume this instead of loaded_data().
  filtered_data <- reactive({
    r   <- loaded_data()
    flt <- dir_filter()
    if (is.null(r) || !r$found) return(r)
    if (is.null(flt) || identical(flt, "All")) return(r)
    cards <- loaded_cards()
    if (is.null(cards) || (cards$pipeline %||% "B") != "B") return(r)
    df <- r$data
    if (!"direction" %in% names(df) || nrow(df) == 0L) return(r)
    dir_col <- toupper(trimws(as.character(df$direction)))
    keep <- !is.na(dir_col) & dir_col == toupper(flt)
    result <- r
    result$data <- df[keep, , drop = FALSE]
    rownames(result$data) <- NULL
    result
  })

  # ── Dataset info panel (setup screen) ────────────────────────────────────
  output$dataset_info_panel <- renderUI({
    ds <- input$dataset %||% init_dataset
    if (is.null(ds) || !nzchar(as.character(ds))) return(NULL)
    meta <- load_dataset_meta(ds,
                              policy_df  = .policy_cache,
                              summary_df = .ds_summ)
    make_dataset_info_panel(meta)
  })

  # ── Chart area: switches between dot plot, funnel, network, A↔B compare ──
  output$chart_area_content <- renderUI({
    tab <- chart_tab()

    # Comparison tab: reads precomputed files, doesn't depend on loaded_data()
    if (tab == "compare_ab") {
      ds <- loaded_cards()$dataset %||% ""
      if (!nzchar(as.character(ds))) return(NULL)
      comp_summ <- load_comparison_summary(ds)
      ovlp_tbl  <- load_overlap_table(ds)
      return(make_compare_tab_content(ds, comp_summ, ovlp_tbl))
    }

    r <- filtered_data()
    if (is.null(r) || !r$found) return(NULL)

    if (tab == "dot_plot") {
      n <- nrow(r$data)
      h <- paste0(max(170L, min(320L, as.integer(n * 19L + 30L))), "px")
      plotOutput("dot_plot", height = h, width = "100%")
    } else if (tab == "funnel") {
      plotOutput("funnel_plot", height = "290px", width = "100%")
    } else {
      tagList(
        plotOutput("network_plot", height = "268px", width = "100%"),
        tags$div(
          class = "network-note",
          "Edges show available overlap relationships. ",
          "For Pipeline B, some edges may be approximated from shared cluster IDs."
        )
      )
    }
  })

  # ── Dot plot ─────────────────────────────────────────────────────────────
  output$dot_plot <- renderPlot({
    r <- filtered_data()
    if (is.null(r) || !r$found || nrow(r$data) == 0) return(invisible(NULL))
    pipeline <- loaded_cards()$pipeline %||% "B"
    make_dot_plot(r$data, pipeline, selected_idx())
  }, res = 96, bg = "white")

  # ── Funnel plot ───────────────────────────────────────────────────────────
  output$funnel_plot <- renderPlot({
    fd <- funnel_data()
    if (is.null(fd)) return(invisible(NULL))
    pipeline <- loaded_cards()$pipeline %||% "B"
    make_funnel_plot(fd, pipeline)
  }, res = 96, bg = "white")

  # ── Network plot ──────────────────────────────────────────────────────────
  output$network_plot <- renderPlot({
    nd  <- network_data()
    sel <- {
      r <- filtered_data()
      if (!is.null(r) && r$found && nrow(r$data) > 0L) {
        idx <- min(max(1L, selected_idx()), nrow(r$data))
        as.character(r$data$term_id[idx])
      } else NULL
    }
    make_network_plot(nd, selected_id = sel)
  }, res = 96, bg = "white")

  # ── Chart tab switching ───────────────────────────────────────────────────
  observeEvent(input$chart_tab, {
    chart_tab(input$chart_tab)
  }, ignoreInit = TRUE)

  # ── Direction filter (Pipeline B only) ───────────────────────────────────
  observeEvent(input$dir_filter, {
    dir_filter(input$dir_filter)
    selected_idx(1L)
  }, ignoreInit = TRUE)

  # ── Page content ─────────────────────────────────────────────────────────
  output$page_content <- renderUI({
    state <- app_state()
    if (state == "setup") {
      ds       <- isolate(input$dataset)  %||% init_dataset
      pipeline <- isolate(input$pipeline) %||% "B"
      preset   <- isolate(input$preset)   %||% "Default"
      return(make_setup_screen(ds, pipeline, preset))
    }
    if (state == "results") {
      return(make_results_screen(loaded_cards(), loaded_data()))
    }
    if (state == "not_found") {
      ds       <- isolate(input$dataset)  %||% init_dataset
      pipeline <- isolate(input$pipeline) %||% "B"
      preset   <- isolate(input$preset)   %||% "Default"
      return(make_not_found_screen(loaded_data(), ds, pipeline, preset))
    }
    make_setup_screen()
  })

  # ── Detail panel (re-renders on row selection) ───────────────────────────
  output$detail_panel <- renderUI({
    r <- filtered_data()
    if (is.null(r) || !r$found || nrow(r$data) == 0) return(NULL)

    idx      <- min(max(1L, selected_idx()), nrow(r$data))
    row      <- r$data[idx, ]
    pipeline <- loaded_cards()$pipeline %||% "B"

    # Look up genes for the selected term from the pre-loaded map
    gm    <- gene_map()
    tid   <- as.character(row$term_id %||% "")
    genes <- if (nzchar(tid) && tid %in% names(gm)) gm[[tid]] else NULL

    tags$div(
      class = "results-card detail-card",
      tags$div(
        class = "results-card-header",
        tags$span(class = "detail-header-name",
                  as.character(row$term_name %||% "—")),
        tags$span(class = "detail-header-id",
                  as.character(row$term_id   %||% "—"))
      ),
      tags$div(
        class = "detail-body",
        make_detail_html(row, pipeline, genes)
      )
    )
  })

  # ── DT shortlist table ───────────────────────────────────────────────────
  output$shortlist_table <- DT::renderDT({
    r <- filtered_data()
    if (is.null(r) || !r$found || nrow(r$data) == 0)
      return(DT::datatable(data.frame()))

    pipeline <- loaded_cards()$pipeline %||% "B"
    disp     <- build_dt_display(r$data, pipeline)

    DT::datatable(
      disp,
      escape       = FALSE,
      selection    = list(mode = "single", selected = 1L),
      rownames     = FALSE,
      class        = "compact stripe hover",
      options      = list(
        paging       = FALSE,
        dom          = "ft",          # filter + table, no pagination/info
        scrollY      = "400px",
        scrollCollapse = TRUE,
        autoWidth    = FALSE,
        order        = list(),        # preserve data order (already ranked)
        columnDefs   = if (pipeline == "B") list(
          list(className = "dt-left",   targets = 0L),
          list(className = "dt-center", targets = c(1L, 2L, 3L, 4L)),
          list(width = "32%", targets = 0L),   # Term
          list(width = "13%", targets = 1L),   # Collection
          list(width = "10%", targets = 2L),   # Direction
          list(width = "14%", targets = 3L),   # q-value
          list(width =  "8%", targets = 4L)    # Genes
        ) else list(
          list(className = "dt-left",   targets = 0L),
          list(className = "dt-center", targets = c(1L, 2L, 3L, 4L, 5L)),
          list(width = "28%", targets = 0L),   # Term
          list(width = "12%", targets = 1L),   # Collection
          list(width = "14%", targets = 2L),   # q-value
          list(width = "12%", targets = 3L),   # Stability
          list(width =  "7%", targets = 4L)    # Genes
        ),
        language = list(search = "Search:")
      )
    )
  }, server = FALSE)

  # ── Row selection: DT → selected_idx → detail panel + dot plot ───────────
  observeEvent(input$shortlist_table_rows_selected, {
    sel <- input$shortlist_table_rows_selected
    if (!is.null(sel) && length(sel) > 0) {
      new_idx <- as.integer(sel[1])
      if (!identical(new_idx, selected_idx())) selected_idx(new_idx)
    }
  }, ignoreInit = FALSE)

  # ── Download Results ──────────────────────────────────────────────────────
  output$download_results <- downloadHandler(
    filename = function() {
      r <- loaded_data()
      c <- loaded_cards()
      ds          <- c$dataset          %||% "dataset"
      pipeline    <- c$pipeline         %||% "B"
      eff_preset  <- r$effective_preset %||% "Default"
      preset_safe <- gsub("[^A-Za-z0-9]", "_",
                          sub("\\s*\\(.*\\)$", "", eff_preset))
      paste0("pathway_shortlist_", ds, "_Pipeline", pipeline,
             "_", preset_safe, ".csv")
    },
    content = function(file) {
      r <- isolate(filtered_data())
      if (is.null(r) || !r$found) return(invisible(NULL))
      df <- r$data

      std_cols <- c("term_id", "term_name", "collection", "direction",
                    "rank", "q_value", "stability", "gene_count",
                    "cluster", "consensus", "source_pipeline")
      present_std <- intersect(std_cols, names(df))
      rest        <- setdiff(names(df), present_std)

      write.csv(df[, c(present_std, rest), drop = FALSE], file, row.names = FALSE)
    }
  )

  # ── Advanced Settings modal ───────────────────────────────────────────────
  observeEvent(input$open_advanced, {
    ds       <- input$dataset  %||% init_dataset
    pipeline <- input$pipeline %||% "B"
    showModal(make_adv_modal(ds, pipeline, adv_settings()))
  })

  observeEvent(input$adv_save_settings, {
    pipeline <- input$pipeline %||% "B"
    curr     <- adv_settings()
    if (pipeline == "B") {
      raw  <- input$adv_b_top_n %||% "5"
      new_n <- if (raw == "All") "All" else as.integer(raw)
      adv_settings(modifyList(curr, list(b_top_n = new_n)))
    } else {
      adv_settings(modifyList(curr, list(
        a_source = input$adv_a_source %||% "Baseline"
      )))
    }
    removeModal()
  }, ignoreInit = TRUE)

  # ── Load shortlist ────────────────────────────────────────────────────────
  observeEvent(input$load_results, {
    ds       <- input$dataset  %||% init_dataset
    pipeline <- input$pipeline %||% "B"
    settings <- adv_settings()

    # Determine which precomputed result to load
    preset <- if (pipeline == "A" && settings$a_source == "Recommended") {
      "Recommended"
    } else {
      "Default"
    }

    result <- load_shortlist(ds, pipeline, preset)

    # Pipeline B: apply top-N per direction filter to precomputed data
    if (pipeline == "B" && result$found && nrow(result$data) > 0) {
      raw_result(result)                             # store unfiltered for re-filter on settings change
      n <- settings$b_top_n
      result$data <- filter_pipelineB_top_n(result$data, n)
      result$filter_note <- if (identical(n, "All")) {
        "Showing all rows in the loaded result file (pipeline already capped at 5 per direction)."
      } else {
        paste0("Showing top ", n, " terms per direction from the precomputed shortlist.")
      }
    } else if (pipeline == "A" && result$found) {
      raw_result(NULL)
      result$filter_note <- if (preset == "Recommended" &&
                                grepl("Recommended", result$effective_preset, fixed = TRUE)) {
        "Loaded recommended (policy-optimised) result."
      } else NULL
    } else {
      raw_result(NULL)
      result$filter_note <- NULL
    }

    cards <- load_summary_cards(ds, pipeline, preset)

    loaded_data(result)
    loaded_cards(cards)
    selected_idx(1L)
    dir_filter("All")
    chart_tab("dot_plot")

    if (result$found && nrow(result$data) > 0) {
      wp_map <- load_pathway_gene_map(ds)
      go_map <- load_go_gene_map(result$data$term_id)
      gene_map(c(wp_map, go_map))
      funnel_data(tryCatch(
        load_funnel_counts(ds, pipeline, preset),
        error = function(e) { warning("load_funnel_counts: ", conditionMessage(e)); NULL }
      ))
      network_data(tryCatch(
        load_network_edges(ds, pipeline, preset, shortlist_df = result$data),
        error = function(e) { warning("load_network_edges: ", conditionMessage(e)); NULL }
      ))
    } else {
      gene_map(list())
      funnel_data(NULL)
      network_data(NULL)
    }

    if (result$found) app_state("results") else app_state("not_found")
  })

  # ── Auto-refilter Pipeline B when Top N setting changes ──────────────────
  # Fires when adv_settings() changes while results are displayed.
  # Uses the stored unfiltered raw_result() so no file re-read is needed.
  observeEvent(adv_settings(), {
    raw <- raw_result()
    if (is.null(raw) || app_state() != "results") return()
    cards <- loaded_cards()
    if (is.null(cards) || (cards$pipeline %||% "B") != "B") return()

    n      <- adv_settings()$b_top_n
    result <- raw
    result$data <- filter_pipelineB_top_n(raw$data, n)
    result$filter_note <- if (identical(n, "All")) {
      "Showing all available precomputed terms."
    } else {
      paste0("Showing top ", n, " terms per direction from the precomputed shortlist.")
    }

    loaded_data(result)
    dir_filter("All")
    selected_idx(1L)
  }, ignoreInit = TRUE)

  # ── Back to Setup ─────────────────────────────────────────────────────────
  observeEvent(input$back_to_setup, {
    raw_result(NULL)
    app_state("setup")
  })
}

shinyApp(ui, server)
