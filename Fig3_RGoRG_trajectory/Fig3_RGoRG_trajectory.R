# ==============================================================================
# Figure 3: RG/oRG Pseudotime, DEG, and TF Regulon Analysis
# Corresponds to: Figure 3A-F
# Input:  subRGoRG.RDS
# Author: Jiaen Lin
# Date: 08-May-2026
# ==============================================================================

library(Seurat); library(SeuratWrappers); library(monocle3)
library(SCENIC); library(SCopeLoomR); library(AUCell)
library(dplyr); library(ggplot2); library(ComplexHeatmap)

path_deposit <- "data/processed/"; path_output <- "results/Fig3/"
dir.create(path_output, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load and prepare ───────────────────────────────────────────────────────
subRGoRG <- LoadSeuratRds(file.path(path_deposit, "subRGoRG.RDS"))
Idents(subRGoRG) <- subRGoRG$CellType
levels(subRGoRG)  <- c("RG","oRG")

# ── 2. Fig 3A: UMAP split by condition ───────────────────────────────────────
p_umap <- DimPlot(subRGoRG, split.by = "orig.ident", cols = c("#1FA195","#F68282"))
ggsave(file.path(path_output, "Fig3A_UMAP_split.pdf"), p_umap, width = 10, height = 5)

# ── 3. Fig 3B: Composition bar ────────────────────────────────────────────────
prop_data <- subRGoRG@meta.data %>%
  group_by(orig.ident, CellType) %>% summarise(n=n(),.groups="drop") %>%
  mutate(percent = n/sum(n)*100)
p_bar <- ggplot(prop_data, aes(x=orig.ident, y=percent, fill=CellType)) +
  geom_col(width=0.7) + coord_flip() + theme_classic() +
  scale_fill_manual(values=c("#1FA195","tomato2"))
ggsave(file.path(path_output, "Fig3B_composition_bar.pdf"), p_bar, width=6, height=4)

# 4. Fig 3C:DEG analysis ───────────────────────────────────────────
DefaultAssay(sub_RGoRG)<-"RNA"
Idents(sub_RGoRG)<-"orig.ident"
DEG<-FindMarkers(sub_RGoRG,ident.1 = "BMP4",ident.2 = "control")
DEG$combined_score<- DEG$avg_log2FC * (-log10(DEG$p_val_adj))
DEG<-DEG[order(DEG$combined_score,decreasing = TRUE),]
label<-append(rownames(head(DEG, 15)),rownames(tail(DEG, 15)))
p_DEG<-EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 10e-5,FCcutoff = 0.5,selectLab = label,labSize = 4.0)
ggsave(file.path(path_output,"Fig3C_pseudotime_featureplot.pdf"), p_DEG, width=5, height=5)

# ── 5. Fig 3D: Monocle3 pseudotime ───────────────────────────────────────────
subRGoRG_rna <- subRGoRG
for (a in c("AUC","AUCBinary","Treat.RSS","RSS")) subRGoRG_rna[[a]] <- NULL

cds <- as.cell_data_set(subRGoRG_rna)
cds <- cluster_cells(cds, resolution = 2e-3)
cds <- learn_graph(cds, use_partition = FALSE)
root_cells <- colnames(cds)[clusters(cds) == 10]   # verify visually
cds <- order_cells(cds, root_cells = root_cells)

p_pseudo <- plot_cells(cds, color_cells_by="pseudotime",
                       label_cell_groups=FALSE, label_branch_points=FALSE,
                       label_roots=FALSE, label_leaves=FALSE, trajectory_graph_color="grey60")
ggsave(file.path(path_output,"Fig3C_pseudotime.pdf"), p_pseudo, width=6, height=5)

subRGoRG$pseudotime_m3 <- pseudotime(cds)[colnames(subRGoRG)]
p_ft <- FeaturePlot(subRGoRG, features="pseudotime_m3", cols=c("grey","red"), order=TRUE)
ggsave(file.path(path_output,"Fig3D_pseudotime_featureplot.pdf"), p_ft, width=5, height=5)

# ── 6. Fig 3E,G: Gene/regulon-pseudotime correlation ───────────────────────────────────
DefaultAssay(subRGoRG) <- "RNA"
expr_scaled <- ScaleData(subRGoRG,
                         features = rownames(subRGoRG)[rowSums(subRGoRG@assays$RNA@data) > 1])
cor_gene <- cor(t(expr_scaled@assays$RNA@scale.data), subRGoRG$pseudotime_m3, method="pearson")
cor_gene <- tibble::rownames_to_column(as.data.frame(cor_gene), "gene") %>%
  rename(cor=1) %>% arrange(desc(cor))
write.csv(cor_gene, file.path(path_output,"Fig3E_RGoRG_pseudotime_cor_genes.csv"), row.names=FALSE)

DefaultAssay(subRGoRG) <- "AUC"
expr_scaled <- ScaleData(subRGoRG,
                         features = rownames(subRGoRG)[rowSums(subRGoRG@assays$RNA@data) > 1])
cor_TF <- cor(t(expr_scaled@assays$RNA@scale.data), subRGoRG$pseudotime_m3, method="pearson")
cor_TF <- tibble::rownames_to_column(as.data.frame(cor_TF), "gene") %>%
  rename(cor=1) %>% arrange(desc(cor))
write.csv(cor_TF, file.path(path_output,"Fig3G_RGoRG_pseudotime_cor_TFs.csv"), row.names=FALSE)

# ── 7. Fig 3F, H: Driver genes and TF regulons identification ────────────────────────────────
DefaultAssay(subRGoRG) <- "RNA"
Idents(subRGoRG) <- subRGoRG$orig.ident
markers <- FindMarkers(subRGoRG, ident.1="BMP4", ident.2="control",
                          only.pos=FALSE, logfc.threshold=0.1)
markers$combined_score <- markers$avg_log2FC * (-log10(markers$p_val))
markers <- tibble::rownames_to_column(markers, "GENE")
write.csv(markers, file.path(path_output,"RGoRG_DEG_BMP4vsCtrl.csv"), row.names=FALSE)

DefaultAssay(subRGoRG) <- "AUC"
Idents(subRGoRG) <- subRGoRG$orig.ident
TF_markers <- FindMarkers(subRGoRG, ident.1="BMP4", ident.2="control",
                          only.pos=FALSE, logfc.threshold=0.1)
TF_markers$combined_score <- TF_markers$avg_log2FC * (-log10(TF_markers$p_val))
TF_markers <- tibble::rownames_to_column(TF_markers, "GENE")
write.csv(TF_markers, file.path(path_output,"RGoRG_TF_BMP4vsCtrl.csv"), row.names=FALSE)

driver_genes <- intersect(
  cor_gene$gene[cor_tf$cor > 0.1 & cor_gene$pvalue<0.05],
  markers$GENE[markers$avg_log2FC > 0.5 & markers$FDR<0.05])

driver_TFs <- intersect(
  cor_tf$gene[cor_tf$cor > 0.1 & cor_tf$pvalue<0.05],
  TF_markers$GENE[TF_markers$avg_log2FC > 0.5 & TF_markers$FDR<0.05])

# ── 8. Fig 3I: TF network building ────────────────────────────────
#Detail code in TF_gene_network_building.r


