# MaRBLe 2.0 — Pathway Redundancy & False-Positive Reduction

**Yustyna Babichuk**

A reproducible R project that turns gene-level differential-expression statistics
into a concise, de-duplicated, and stable shortlist of enriched pathways. It
provides two complementary enrichment pipelines, a threshold-optimisation and
simulation analysis, and an interactive Shiny dashboard for exploring the
results.

The repository ships with precomputed results for every bundled dataset, so the
**dashboard runs straight from a fresh clone** — no need to re-run the pipelines
first.

---

## Quick start — view the dashboard

The dashboard runs from a fresh clone with **only three packages**
(`shiny`, `DT`, `ggplot2`) — you do **not** need `renv` or the heavy
Bioconductor pipeline stack just to view the results.

From the **project root**, one command:

```bash
Rscript run_dashboard.R
```

This installs the three packages if they are missing and opens the app in your
browser. Equivalently, from an R session:

```r
install.packages(c("shiny", "DT", "ggplot2"))
shiny::runApp("shiny")
```

Pick a dataset and a pipeline, click **Load shortlist**, and explore the dot
plot, selection funnel, overlap network, and the A ↔ B comparison tab.

> Optional: installing `org.Hs.eg.db` + `AnnotationDbi` (Bioconductor) also
> shows the gene members behind each GO term. The dashboard works fine without
> them — it just skips that lookup.

## Re-running the pipelines (optional)

Re-running Pipeline A/B needs the full analysis stack, pinned with
[renv](https://rstudio.github.io/renv/):

```r
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::restore()   # heavy: CRAN + Bioconductor (clusterProfiler, fgsea, limma, …)
```

> **Note:** if `renv::restore()` fails on `RcppArmadillo` (a gfortran compile
> error on newer macOS), install the affected CRAN packages as binaries:
> `install.packages(<pkg>, type = "binary")`.

Then edit `config/default.yml` to select a dataset (exactly one `dataset:`
block uncommented — `dataset_0` is the default) and run:

```bash
Rscript R/run_pipelineA.R   # ORA + redundancy reduction + bootstrap + tiers
Rscript R/run_pipelineB.R   # CAMERA + fgsea + semantic collapse + top-N shortlist
```

Results are written to `results/pipelineA/<run_dir>/` or
`results/pipelineB/<run_dir>/`.

---

## The two pipelines

| | Pipeline A | Pipeline B |
|---|---|---|
| **Method** | ORA (over-representation) | CAMERA + fgsea (competitive enrichment) |
| **Input** | DE gene list + universe | Full expression matrix |
| **Redundancy reduction** | Semantic collapse (GO) + Jaccard clustering | Semantic collapse + Jaccard clustering |
| **Stability** | Bootstrap resampling (B = 200) | — |
| **Output** | Tier 1 (stable + consensus) → Tier 2 (stable) → candidates | Top-N per direction (default 5 UP + 5 DOWN) |
| **Entry point** | `R/run_pipelineA.R` | `R/run_pipelineB.R` |

The two methods look at the data through different lenses and are meant to be
read together: Pipeline A gives a stable ORA consensus (it needs a robust DE
list), while Pipeline B works without a DE pre-filter and adds an UP/DOWN
direction.

---

## Repository structure

```
run_dashboard.R            # One-command dashboard launcher
R/
  run_pipelineA.R          # Pipeline A entry point
  run_pipelineB.R          # Pipeline B entry point
  pipelineA/steps/         # Step scripts 01–07 + step98 (final export)
  pipelineB/steps/         # Step scripts 01–06
  analysis/                # Post-run analysis: diagnostics, threshold sweep,
  │                        #   policy, simulations, plots
  ablation/                # Pipeline B ablation experiment scripts
  utils/                   # Shared utilities (config, gene-ID conversion, …)
config/
  default.yml              # Dataset selection, pipeline parameters, thresholds
  ablation/                # Ablation variant configs
data/
  dataset_0/ … dataset_3/  # Reference RNA-seq datasets
  ACTA2KO_iVSMC/           # ACTA2-knockout iVSMC application dataset
results/                   # Precomputed pipeline + analysis outputs (curated)
report_assets/             # Publication-ready figures (PNG) and source tables
shiny/                     # Interactive dashboard (see shiny/README.md)
renv.lock                  # Locked package versions
```

---

## Datasets

| Dataset | Description | Gene IDs |
|---|---|---|
| `dataset_0` | Reference cohort — full transcriptome | SYMBOL |
| `dataset_1` | GSE199939 — EV / KD groups | ENTREZ |
| `dataset_2` | GSE243836 — focused panel (~2 200 genes) | ENTREZ |
| `dataset_3` | GSE247345 | ENTREZ |
| `ACTA2KO_iVSMC` | ACTA2 knockout vs wild-type in iVSMC | ENSEMBL → SYMBOL |

All datasets share a single pathway annotation,
`data/dataset_0/processed/Pathway2Gene.csv` (GO-BP + WikiPathways).

> Additional vascular smooth-muscle-cell datasets used during development are
> **withheld pending publication** and are intentionally not included here.

---

## Analysis scripts (`R/analysis/`)

| Script | Purpose |
|---|---|
| `00_diagnostics.R` | Per-dataset / per-pipeline summary table |
| `01_sweep_pipelineA.R` | Threshold-sensitivity sweep for Pipeline A |
| `02_expressed_fraction.R` | Expression-coverage QC for final pathways |
| `03_threshold_policy.R` | Dataset classification + recommended thresholds |
| `04_run_with_policy.R` | Run Pipeline A with auto-selected thresholds |
| `05_simulation.R` | Simulation studies A (coverage) and B (signal) |
| `06_plot_stability_vs_significance.R` | Stability-vs-significance figure |

`R/ablation/` contains the Pipeline B component-ablation experiment
(`run_ablation_B.R` + comparison scripts) with its variant configs under
`config/ablation/`.

---

## Configuration (`config/default.yml`)

Key parameters under the `pipelineA:` / `pipelineB:` blocks and shared `thresholds:`:

- `dataset:` — the active dataset (expression + metadata + statistics files)
- `gene_id_type:` — `SYMBOL` | `ENTREZ` | `ENSEMBL` (Pipeline B normalises to SYMBOL)
- `TAU_STABILITY:` — bootstrap stability threshold (Pipeline A, default 0.70)
- `CONS_JACCARD_MIN:` — cross-collection consensus Jaccard threshold (default 0.15)
- `agreement_mode:` — `intersection` | `relax_fdr_then_topn` (Pipeline B step 03)

---

## Environment

Managed with renv. After installing or updating packages, refresh the lockfile:

```r
renv::snapshot()
```

---

## License

Academic project. Code is shared for review and reference; please get in touch
before reusing the code or data.
