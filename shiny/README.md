# Pathway Analysis Dashboard

Interactive R/Shiny dashboard for browsing precomputed Pipeline A and Pipeline B pathway enrichment shortlists across multiple RNA-seq datasets.

## How to run

From the **project root** (`marble/`):

```bash
Rscript -e "shiny::runApp('shiny')"
```

The app auto-detects available datasets from `data/` and opens on the setup screen.

## Required R packages

| Package | Purpose |
|---------|---------|
| `shiny` | App framework |
| `DT` | Interactive shortlist table |
| `ggplot2` | Dot plot, selection funnel, overlap network |

Install if needed:

```r
install.packages(c("shiny", "DT", "ggplot2"))
```

## Expected project structure

```
marble/
├── data/
│   ├── dataset_0/processed/Pathway2Gene.csv   ← shared gene map
│   └── dataset_0/, dataset_1/, dataset_2/, dataset_3/
├── results/
│   ├── pipelineA/<RUN_DIR>/FINAL/FINAL.csv
│   ├── pipelineB/<RUN_DIR>/FINAL/FINAL.csv
│   ├── analysis/
│   │   ├── dataset_summary.csv
│   │   ├── threshold_policy.csv
│   │   └── policy_runs/<dataset>/FINAL.csv   ← Pipeline A Recommended source
│   └── run_registry.csv
└── shiny/
    ├── app.R
    ├── R/  (data_loader.R, dot_plot.R, funnel_plot.R, network_plot.R, ui_helpers.R)
    └── www/marble.css
```

## Dashboard features

### Setup screen
- Dataset selection — auto-detected from `data/` subdirectories
- Pipeline selection — Pipeline A (ORA + bootstrap) or Pipeline B (CAMERA + fgsea); defaults to Pipeline B
- Advanced Settings button — opens pipeline-specific settings modal:
  - **Pipeline B**: Top terms per direction — 3 / 5 / All available (default 5)
  - **Pipeline A**: Result source — Baseline or Recommended policy

### Results screen
- **Five summary cards**: Dataset, Pipeline, Final terms, Universe (genes), Direction balance (B) or Policy class (A)
- **Shortlist table** (DT, searchable):
  - Pipeline B columns: Term, Collection, Direction, q-value, Genes
  - Pipeline A columns: Term, Collection, q-value, Stability, Genes, Consensus
- **Dot plot** — terms by significance; colored by direction (B) or stability (A); selected row highlighted
- **Selection funnel** — term counts at each pipeline reduction stage
- **Overlap network** — shortlisted terms as nodes; overlap/consensus edges; labels shown for selected node, its neighbors, and top-3 by gene count
- **Detail panel** — shown on row selection; displays collection, direction, q-value, stability/reason, gene count, consensus, and gene list preview
- **Parameters strip** — active FDR, stability, Jaccard, and bootstrap thresholds
- **Download shortlist CSV** — exports full standardised data frame including all pipeline columns

## v1 limitations

| Limitation | Notes |
|---|---|
| Loads precomputed results only | Does not rerun Pipeline A or Pipeline B |
| Pipeline B Top N filters display only | Only trims from the precomputed cap (default: 5 per direction); does not rerun pipeline |
| Pipeline A Recommended source | Loads `results/analysis/policy_runs/<dataset>/FINAL.csv` if available; falls back to Baseline |
| GO gene preview | Requires `AnnotationDbi` + `org.Hs.eg.db` (Bioconductor); shows fallback note if unavailable |
| Overlap network edges | Depends on cluster/edge files in the results directory; shows nodes-only if no edges found |
| No cross-pipeline comparison | Pipeline A and B results are browsed separately |
