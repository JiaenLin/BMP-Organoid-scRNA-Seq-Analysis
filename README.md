# BMP-Organoid-scRNA-Seq-Analysis
Reproducibility repository for the Paper: BMP Signaling Modulates Human Cortical Expansion and Maturation
The repository is organized based on the main figures in the manuscript.
## System Requirements

### Software Dependencies

| Software | Version | Purpose |
|---|---|---|
| R | 4.4.0 | All R-based analyses |
| Seurat | v5 | scRNA-seq clustering, annotation, DEG |
| Monocle3 | 3.0 | Pseudotime trajectory inference |
| SCENIC / pySCENIC | 0.12.1 | Gene regulatory network (GRN) analysis |
| clusterProfiler | v4.0.5 | Gene set enrichment analysis (GSEA) |
| scANVI | 1.0 | Reference mapping / label transfer |
| Python | ≥ 3.8 | Fig 8 label transfer pipeline |
| scanpy | 1.12.1| AnnData I/O and preprocessing (Fig 8) |

> **Tested on:** Linux (HPC cluster), x86-64. No non-standard hardware required.

### Databases Required
- **MSigDB** — for GSEA pathway enrichment (`clusterProfiler`)
- **cisTarget databases** — for RcisTarget motif enrichment (pySCENIC)
- **Human developing neocortex atlas** — reference AnnData for scANVI label transfer (Fig 8); see Wang, L., Wang, C., Moriano, J.A. et al. Molecular and cellular dynamics of the developing human neocortex. Nature 647, 169–178 (2025).
## Installation Guide

1. **Clone this repository:**
   ```bash
   git clone <repo-url>
   cd BMP4_organoid_analysis
   ```

2. **Install R packages** (in R ≥ 4.4.0):
   ```r
   install.packages(c("Seurat", "ggplot2", "dplyr", "pheatmap", "ggpubr"))
   BiocManager::install(c("clusterProfiler", "AUCell", "RcisTarget", "GENIE3", "monocle3"))
   ```

3. **Install Python packages** (for Fig 8):
   ```bash
   pip install scanpy scvi-tools scipy pandas seaborn matplotlib
   ```

4. **Typical install time:** ~20–40 minutes on a standard desktop computer.

---

## Data Availability

Raw and processed scRNA-seq data are deposited in NCBI Gene Expression Omnibus (GEO):

**Accession: [GSE306803]**
## Repository Structure

```
BMP4_paper_code/
├── Fig2_scRNAseq_clustering/
│   └── Fig2_scRNAseq_clustering.R        # QC, merge, normalise, cluster, annotate, UMAP,
│                                          # cell composition, whole-organoid SCENIC heatmap
├── Fig3_RGoRG_trajectory/
│   └── Fig3_RGoRG_trajectory.R           # RG/oRG pseudotime, TF correlation,
│                                          # RSS heatmap, differential TF regulons
├── Fig6_NPC_GABA_trajectory/
│   └── Fig6_NPC_GABA_trajectory.R        # GABA lineage Monocle3 trajectory,
│                                          # gene/TF-pseudotime correlation, ARX/DLX plots
├── Fig8_reference_mapping/
   └── Fig8_reference_mapping.py         # scANVI label transfer, prediction score dist.,
                                         # fraction-overlap heatmaps, pseudobulk Pearson-R,
                                          # paired Wilcoxon boxplot
```
## Methods Summary

Detailed methods are described in the manuscript (Methods section).
## Demo
For a quick demo, download the processed matrix files from GEO (GSE306803)

## Reproduction Instructions

To reproduce all quantitative results in the manuscript:

1. Download raw count matrices from GEO (GSE306803).
2. Run Cell Ranger (human GRCh38 reference) to generate per-sample filtered matrix files.
3. Run scripts in figure order (Fig 2 → Fig 8); each script loads the processed object from the previous step where applicable.
4. All random seeds are set within each script (`set.seed()` / `np.random.seed(42)`) for reproducibility.
5. Statistical thresholds: p < 0.05 (*), p < 0.01 (**), p < 0.001 (***), p < 0.0001 (****) throughout.

---

## License

MIT License. See `LICENSE` file for details.

---
---
