<div align="center">

# MaRBLe 2.0 — Pathway Analysis of RNA-seq Data
### Reducing Redundancy and False Positives

**Yustyna Babichuk**
BSc Data Science & Artificial Intelligence · Maastricht University
Honours research project — *Maastricht Research Based Learning (MaRBLe 2.0)*

Supervisors: Rachel Cavill · Jarno Koetsier · Pepijn Saraber

![R](https://img.shields.io/badge/R-%E2%89%A5%204.5-276DC3?logo=r&logoColor=white)
![Dashboard](https://img.shields.io/badge/dashboard-Shiny-1BA39C?logo=rstudio&logoColor=white)
![Reproducible](https://img.shields.io/badge/environment-renv-F46800)
![Programme](https://img.shields.io/badge/Maastricht-MaRBLe%202.0%20Honours-002b5c)

<br>

<img src="report_assets/dashboard/dashboard_overview.png" alt="MaRBLe pathway dashboard — Pipeline A results view" width="920">

<sub><em>The Shiny dashboard — Pipeline A results for <code>dataset_0</code>: final shortlist, dot plot, and per-pathway detail.</em></sub>

</div>

---

## Overview

Pathway enrichment analysis is widely used to interpret RNA-seq experiments, but
it has two well-known weaknesses: enriched-pathway lists are often **long and
redundant** (overlapping gene sets, hierarchical GO terms), and many methods
**ignore gene–gene correlation**, which inflates pathway-level significance and
produces false positives.

This project develops and compares **two complementary pathway-analysis
pipelines** that reduce redundancy and false positives while preserving
biologically meaningful signal, and wraps the final outputs in an **interactive
Shiny dashboard** for inspection and comparison.

It was carried out over two semesters of the MaRBLe honours programme:

- **Semester 1** — designed and validated the two pipelines on a reference RNA-seq dataset (proof of concept).
- **Semester 2** — generalised the framework across multiple datasets, added a dataset-aware threshold policy, simulation studies, a Pipeline B ablation study, a cross-pipeline comparison, and the dashboard.

> 📄 The work is written up in two reports (*Semester 1* and *Semester 2*). This
> repository contains the code, the publicly shareable datasets, and the
> precomputed results that back those reports.

---

## ▶️ Run the dashboard

The dashboard runs **straight from a fresh clone** — the precomputed results are
included, so there is nothing to recompute. You only need three R packages
(`shiny`, `DT`, `ggplot2`); the heavy Bioconductor pipeline stack is **not**
required just to view results.

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

In the app: pick a **dataset** and a **pipeline**, click **Load shortlist**, and
explore the dot plot, selection funnel, overlap network, per-pathway detail
panel, and the **A ↔ B comparison** tab. The displayed shortlist can be exported
as CSV.

> Optional: installing `org.Hs.eg.db` + `AnnotationDbi` (Bioconductor) additionally
> shows the gene members behind each GO term. The dashboard works fine without
> them — it just skips that lookup.

<div align="center">

<img src="report_assets/dashboard/dashboard_comparison.png" alt="Pipeline A ↔ B comparison tab" width="920">

<sub><em>The <strong>A ↔ B</strong> tab (<code>dataset_0</code>): the complementarity summary, shared pathway IDs across both pipelines, and per-pathway detail.</em></sub>

</div>

---

## Methods at a glance

**Pipeline A — trust-first over-representation analysis**

```
DE gene list ─▶ ORA (GO-BP + WikiPathways) ─▶ semantic collapse (GO, Wang)
            ─▶ overlap clustering (Jaccard) ─▶ bootstrap stability (B = 200)
            ─▶ GO↔WP consensus ─▶ Tier 1 / Tier 2 / candidates
```

**Pipeline B — correlation-aware, ranked, cutoff-free**

```
expression matrix ─▶ CAMERA  +  fgsea ─▶ agreement filter (both, same direction)
                  ─▶ semantic collapse ─▶ overlap clustering
                  ─▶ top-N per direction (UP / DOWN shortlist)
```

| | **Pipeline A** | **Pipeline B** |
|---|---|---|
| **Strategy** | ORA (over-representation) | CAMERA + fgsea (competitive / ranked) |
| **Input** | DE gene list + measured-gene universe | full expression matrix |
| **Handles gene–gene correlation** | indirectly | **yes** (CAMERA) |
| **Redundancy reduction** | semantic collapse + Jaccard clustering | semantic collapse + Jaccard clustering |
| **Stability** | bootstrap resampling (B = 200) | — |
| **Direction (UP/DOWN)** | — | **yes** |
| **Output** | Tier 1 (stable + consensus) → Tier 2 (stable) → candidates | top-N per direction (default 5 UP + 5 DOWN) |
| **Entry point** | `R/run_pipelineA.R` | `R/run_pipelineB.R` |

The two pipelines are **complementary, not competing**: Pipeline A is strongest
when a reliable DE gene list is available, while Pipeline B adds directionality
and stays informative when gene-level differential expression is weak.

---

## Key findings

- **The framework generalises across datasets**, but the reliability of the final shortlist tracks dataset properties — DE signal strength, pathway coverage, and GO↔WP overlap — rather than failing in one uniform way.
- **Dataset-aware thresholds recover borderline cases.** Relaxing the consensus Jaccard threshold while keeping the stability requirement strict raised `dataset_3` from 2 → 8 Tier-1 pathways, and `ACTA2KO` from 4 → 9, without weakening stability.
- **Redundancy reduction is not cosmetic.** A Pipeline B ablation showed that *overlap clustering* most strongly shapes the candidate pool, while the final top-N cap is mainly a presentation choice.
- **A and B agree rarely but usefully.** Exact pathway-ID overlap is low (0–2 per dataset); where it occurs it carries extra confidence because two different designs converge.
- **A dashboard ties it together**, turning multiple output files into one interactive inspection and comparison layer.

---

## Repository layout

```
run_dashboard.R            One-command dashboard launcher
R/
  run_pipelineA.R          Pipeline A entry point
  run_pipelineB.R          Pipeline B entry point
  pipelineA/steps/         Step scripts 01–07 + step98 (final export)
  pipelineB/steps/         Step scripts 01–06
  analysis/                Diagnostics, threshold sweep & policy, simulations, plots
  ablation/                Pipeline B component-ablation experiment
  utils/                   Shared utilities (config, gene-ID conversion, …)
config/
  default.yml              Active dataset, pipeline parameters, thresholds
  ablation/                Ablation variant configs
data/                      Publicly shareable datasets (expression + metadata + DE stats)
results/                   Precomputed pipeline + analysis outputs (curated, dashboard-ready)
report_assets/             Publication-ready figures (PNG) and source tables
shiny/                     Interactive dashboard  (see shiny/README.md)
renv.lock                  Locked package versions for full reproducibility
```

---

## Datasets

| Dataset | Description | Gene IDs |
|---|---|---|
| `dataset_0` | Reference cohort — full transcriptome (Semester-1 baseline) | SYMBOL |
| `dataset_1` | GSE199939 | ENTREZ |
| `dataset_2` | GSE243836 — focused panel (~2,200 genes; low-coverage case) | ENTREZ |
| `dataset_3` | GSE247345 | ENTREZ |
| `ACTA2KO_iVSMC` | ACTA2 knockout vs wild-type in iVSMC (Biochemistry collaboration) | ENSEMBL → SYMBOL |

All datasets share one pathway annotation,
`data/dataset_0/processed/Pathway2Gene.csv` (GO Biological Process + WikiPathways).
The enrichment universe is always the intersection of measured genes and
annotated pathway genes, so enrichment is only tested against genes that could
have been observed.

> Additional vascular smooth-muscle-cell datasets from the Biochemistry
> collaboration (`Aneurysm_pVSMC`, `PAR1KO_iVSMC`, `Lineages_iVSMC`) were part of
> the study but are **withheld pending publication** and are intentionally not
> included in this public repository.

---

## Re-running the pipelines (optional)

Re-running Pipeline A/B requires the full analysis stack, pinned with
[renv](https://rstudio.github.io/renv/):

```r
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::restore()   # heavy: CRAN + Bioconductor (clusterProfiler, fgsea, limma, …)
```

> **Note:** if `renv::restore()` fails on `RcppArmadillo` (a gfortran compile
> error on newer macOS), install the affected CRAN packages as binaries:
> `install.packages(<pkg>, type = "binary")`.

Then edit `config/default.yml` to select a dataset (exactly one `dataset:` block
uncommented — `dataset_0` is the default) and run:

```bash
Rscript R/run_pipelineA.R   # ORA + redundancy reduction + bootstrap + tiers
Rscript R/run_pipelineB.R   # CAMERA + fgsea + semantic collapse + top-N shortlist
```

Outputs are written to `results/pipelineA/<run_dir>/` or `results/pipelineB/<run_dir>/`.

### Post-run analysis (`R/analysis/`)

| Script | Purpose |
|---|---|
| `00_diagnostics.R` | Per-dataset / per-pipeline summary table |
| `01_sweep_pipelineA.R` | Threshold-sensitivity sweep for Pipeline A |
| `02_expressed_fraction.R` | Expression-coverage QC for final pathways |
| `03_threshold_policy.R` | Dataset classification + recommended thresholds |
| `04_run_with_policy.R` | Run Pipeline A with auto-selected thresholds |
| `05_simulation.R` | Simulation studies (coverage & signal) for threshold behaviour |
| `06_plot_stability_vs_significance.R` | Stability-vs-significance figure |

`R/ablation/` holds the Pipeline B component-ablation experiment
(`run_ablation_B.R` + comparison scripts), configured by `config/ablation/`.

### Key configuration (`config/default.yml`)

- `dataset:` — the active dataset (expression + metadata + statistics files)
- `gene_id_type:` — `SYMBOL` | `ENTREZ` | `ENSEMBL` (Pipeline B normalises to SYMBOL)
- `TAU_STABILITY:` — bootstrap stability threshold (Pipeline A, default `0.70`)
- `CONS_JACCARD_MIN:` — cross-collection consensus Jaccard threshold (default `0.15`)
- `agreement_mode:` — `intersection` | `relax_fdr_then_topn` (Pipeline B step 03)

---

## Acknowledgements

Carried out within the **MaRBLe 2.0** honours programme at Maastricht
University, in collaboration with the Department of Biochemistry. Thanks to
supervisors **Rachel Cavill**, **Jarno Koetsier**, and **Pepijn Saraber** for
their guidance throughout the project.

## License

Academic project shared for review and as a portfolio reference. Please get in
touch before reusing the code or data.
