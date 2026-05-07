# BMP-Organoid-scRNA-Seq-Analysis
Reproducibility repository for the Paper: BMP Signaling Modulates Human Cortical Expansion and Maturation
The repository is organized based on the main figures in the manuscript.
## Directory Structure

BMP4_paper_code/
├── Fig2_scRNAseq_clustering/
│   └── Fig2_scRNAseq_clustering.R
├── Fig3_RGoRG_trajectory/
│   └── Fig3_RGoRG_trajectory.R
├── Fig4_NPC_GABA_trajectory/
│   └── Fig4_NPC_GABA_trajectory.R
├── Fig5_fetal_brain/
│   └── Fig5_fetal_reference_mapping.R
├── Fig6_NPC_glutamatergic/
│   └── Fig6_NPC_Glu_trajectory.R
├── Fig7_GABA_regulons/
│   └── Fig7_GABA_regulon_analysis.R
├── utils/
│   └── pySCENIC_loom_prep.R
└── README.md

## Figure-to-Script Mapping

| Figure   | Script                          | Key analysis                               |
|----------|---------------------------------|--------------------------------------------|
| Fig 2A–C | Fig2_scRNAseq_clustering.R      | UMAP, annotation, marker DotPlot           |
| Fig 2D–F | Fig2_scRNAseq_clustering.R      | Cell composition bar charts                |
| Fig 2G   | Fig2_scRNAseq_clustering.R      | Whole-organoid regulon heatmap (SCENIC)    |
| Fig 3A–B | Fig3_RGoRG_trajectory.R         | RG/oRG UMAP + composition                 |
| Fig 3C   | Fig3_RGoRG_trajectory.R         | Monocle3 pseudotime                        |
| Fig 3D   | Fig3_RGoRG_trajectory.R         | Gene-pseudotime correlation scatter        |
| Fig 3E–F | Fig3_RGoRG_trajectory.R         | Differential TF regulons + RSS heatmap     |
| Fig 4A–B | Fig4_NPC_GABA_trajectory.R      | GABA UMAP + composition + pseudotime       |
| Fig 4C   | Fig4_NPC_GABA_trajectory.R      | Pseudotime-correlated genes DotPlot        |
| Fig 4D   | Fig4_NPC_GABA_trajectory.R      | DEGs in early GABA trajectory              |
| Fig 4E–F | Fig4_NPC_GABA_trajectory.R      | ARX/DLX regulon DotPlot + FeaturePlots     |
| Fig 5A   | Fig5_fetal_reference_mapping.R  | Reference vs organoid UMAP                 |
| Fig 5B–E | Fig5_fetal_reference_mapping.R  | Fetal markers, oRG validation              |
| Fig 6A–C | Fig6_NPC_Glu_trajectory.R       | NPC-Glu pseudotime + correlated genes      |
| Fig 6D–E | Fig6_NPC_Glu_trajectory.R       | NPC marker DotPlots                        |
| Fig 7A   | Fig7_GABA_regulon_analysis.R    | RSS heatmap (all cell types)               |
| Fig 7B   | Fig7_GABA_regulon_analysis.R    | Pseudo-bulk AUC heatmap                    |
| Fig 7C–E | Fig7_GABA_regulon_analysis.R    | Per-cell-type DEG + dot-plots   
