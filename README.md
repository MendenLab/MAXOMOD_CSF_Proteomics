# MAXOMOD_CSF_Proteomics

This repository provides demo data and scripts to reproduce the bioinformatic pipelines and figures from the manuscript *Proteomic Profiling of Cerebrospinal Fluid Identifies Immune and Synapto-Axonal ALS Subtypes*.

---

# Project structure

| Path | Description |
|------|-------------|
| `Script/` | R scripts for the full analysis pipeline (00–20) |
| `Analysis.sh` | Master script: Discovery → Validation → between-cohort → External → Brain |
| `vignettes/` | R Markdown files to reproduce selected main/supplementary figures |
| `demo/` | Demo inputs (copy of pipeline outputs) used by vignettes |
| `Discovery/`, `Validation/` | Pipeline outputs for each cohort |
| `Input/` | Reference data (e.g. ID mapping, annotation files) |
| `Plots/` | Default figure output directory for vignettes |

---

# System requirements

## R and RStudio
- R v.4.2.3 (or compatible). Install CRAN packages with `install.packages()`. For Bioconductor packages (e.g. `DEP`, `SummarizedExperiment`, `clusterProfiler`, `org.Hs.eg.db`) use `BiocManager::install("PackageName")`.

### Packages

- AnnotationDbi
- biomaRt
- boot
- caret
- circlize
- cluster
- clusterProfiler
- ComplexHeatmap
- cowplot
- data.table
- dendextend
- DEP
- dplyr
- dunn.test
- e1071
- GOSemSim
- ggalluvial
- ggpubr
- ggrepel
- ggsci
- ggplot2
- ggthemes
- glmnet
- grid
- gridExtra
- GSVA
- limma
- mclust
- msigdbr
- NMF
- nnet
- nortest
- openxlsx
- optparse
- org.Hs.eg.db
- patchwork
- pheatmap
- pROC
- purrr
- randomForest
- readr
- readxl
- reshape2
- RColorBrewer
- scales
- stats
- stringr
- SummarizedExperiment
- tibble
- tidyr
- tidyverse
- umap
- visdat
- WGCNA
- writexl
- xgboost

## Operating System:
macOS Tahoe 26.2 (Apple Silicon, M3 Pro)

## Any required non-standard hardware
None for the custom codes.

## Instruction for use
Codes to generate all main and supplementary figures are in the `Script/` folder. Run the full pipeline from the **project root**:
```bash
bash Analysis.sh
```
Or in the background: `nohup bash Analysis.sh > Analysis.log 2>&1 &`  
Paths in scripts are relative to the project root; do not run from another directory.

## Demo
Demo data are under `demo/` (e.g. `demo/Discovery/`, `demo/Validation/`). Populate them by copying from pipeline outputs, or run the pipeline first. Then open and knit the R Markdown files from the project root (or from `vignettes/`; the documents set the project root automatically).

| Figure | Vignette | Description |
|--------|----------|-------------|
| Fig. 1b | `vignettes/Fig1b_clinical_heatmap.Rmd` | Heatmap of the 100 most variable CSF protein abundances in the discovery cohort (run time: &lt; 1 min). |
| Fig. 1c, 1d, 1e | `vignettes/Fig1cde_volcano_scatter.Rmd` | Discovery and validation volcano plots (1c, 1d) and Discovery vs Validation signed FDR scatter (1e). Uses `demo/Discovery/` and `demo/Validation/` DE results (run time: < 1 min).|
| Fig. 1f | `vignettes/Fig1f_GESA_IC.Rmd` | GSEA IC heatmap (GO terms clustered by semantic similarity) for Discovery and Validation. Uses `demo/Discovery/06_GSEA/` and `demo/Validation/06_GSEA/` (run time: < 1 min).|
| Fig. 2b, 2c | `vignettes/Fig2bc_Clustering.Rmd` | AIC by clustering method and number of clusters (2b); Sankey plot of cluster assignment from k=2 to k=3 (2c). Exports `cluster_assignments_2.csv`. Uses `demo/Discovery/02_Missing_Inspection_subclusters/` (run time: < 1 min).|
| Fig. 2d, 2e, 2f | `vignettes/Fig2def_volcano_scatter.Rmd` | Subclusters (alpha vs beta): Discovery volcano (2d), Validation volcano (2e), Discovery vs Validation signed FDR scatter (2f). Uses `demo/Discovery/` and `demo/Validation/` subcluster DE results (`03_Differential_expression_analysis_subclusters/`) (run time: < 1 min).|
| Fig. 3a–3g | `vignettes/Fig3_WGCNA.Rmd` | WGCNA: dendrogram and module–trait bars (3a), eigengene networks (3b), top GO terms per module IC heatmap (3c), eigenprotein boxplots for turquoise and blue modules in Discovery/Validation/External (3d), module–trait heatmaps Discovery (3e) and Validation (3f), MEs vs clinical variables (3g). Saves PDFs to `Plots/Fig3_WGCNA/`. Uses `demo/Discovery/02_Missing_Inspection_subclusters/`, `08_Clustering_als/`, and optionally Validation/External inputs (run time: ~ 8 min).|

## About issues:
To report bugs, ask questions, or provide feedback, please use the GitHub Issues page.