# ==============================================================================
# Figure 6: NPC-to-GABA Neuron Trajectory, Pseudotime Correlations, TF Regulons
# Corresponds to: Figure 6A-F  |  Input: GABAtrac.RDS
# Author: Jiaen Lin
# ==============================================================================

library(Seurat); library(SeuratWrappers); library(monocle3)
library(SCENIC); library(SCopeLoomR); library(AUCell)
library(dplyr); library(ggplot2)

path_deposit <- "data/processed/"; path_output <- "results/Fig6/"
dir.create(path_output, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load GABA lineage subset ───────────────────────────────────────────────
GABAtrac <- LoadSeuratRds(file.path(path_deposit, "GABAtrac.RDS"))
Idents(GABAtrac) <- GABAtrac$CellType
levels(GABAtrac) <- c("NPC","GABA.NPC","Imm.GABA.N","Cajal","GABA.N")

# ── 2. Fig 6A-B: UMAP + cell composition ───────────────────────────────────────
p_umap_split <- DimPlot(GABAtrac, label = TRUE, split.by = "orig.ident")
ggsave(file.path(path_output,"Fig6A_UMAP_split.pdf"), p_umap_split, width = 10, height = 5)

b_num <- sum(GABAtrac$orig.ident == "BMP4")
c_num <- sum(GABAtrac$orig.ident == "control")
counts <- GABAtrac@meta.data %>%
  group_by(orig.ident, Idents = CellType) %>% summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = ifelse(orig.ident == "BMP4", count/b_num*100, count/c_num*100))

p_bar <- ggplot(counts, aes(Idents, percentage, fill = orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = round(percentage, 2)), vjust = -1,
            position = position_dodge(0.8), colour = "black", size = 5) +
  theme(axis.text.x  = element_text(angle = 10, hjust = 1, size = 15),
        axis.text.y  = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.text  = element_text(size = 20),
        axis.line    = element_line(linewidth = 0))
ggsave(file.path(path_output,"Fig6B_composition.pdf"), p_bar, width = 8, height = 5)

# ── 3. Fig 6A: Monocle3 pseudotime ───────────────────────────────────────────
GABAtrac_rna <- GABAtrac
GABAtrac_rna[["AUC"]] <- NULL; GABAtrac_rna[["AUCBinary"]] <- NULL

cds <- as.cell_data_set(GABAtrac_rna)
cds <- cluster_cells(cds, resolution = 2e-3)
cds <- learn_graph(cds, use_partition = FALSE)
root_cells <- colnames(cds)[clusters(cds) == 16]   # NPC root; need to verify visually
cds <- order_cells(cds, root_cells = root_cells)

p_pseudo <- plot_cells(cds, color_cells_by = "pseudotime",
                       label_cell_groups = FALSE, label_branch_points = FALSE,
                       label_roots = FALSE, label_leaves = FALSE, trajectory_graph_color = "grey60")
ggsave(file.path(path_output,"Fig6A_pseudotime.pdf"), p_pseudo, width = 6, height = 5)

# ── 4. Fig 6C: heatmap of regulon activity ───────

regulon_bulk<-AverageExpression(GABAtrac,group.by = 'Treat.CellType',return.seurat = TRUE)
mat1<-regulon_bulk@assays$AUC$scale.data[,rev(c( "BMP4.Cajal",'control.Cajal','BMP4.Imm.GABA.N',"control.Imm.GABA.N",
                                                 'BMP4.GABA.NPC',"control.GABA.NPC",'BMP4.NPC','control.NPC'))]
library(circlize)
circos.clear()
circos.par(gap.after = c(2, 2, 2, 2, 25))
split = sample(letters[1:5], 112, replace = TRUE)
split = factor(split, levels = letters[1:5])
col_fun1 = colorRamp2(c(-2, 0, 2), c("darkgreen", "white", "darkred"))
circos.heatmap(mat1,split = split,col =col_fun1,rownames.side = "outside", rownames.cex = 0.5,
               cluster = TRUE)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 5) { # the last sector
    cn = rev(colnames(mat1))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                c(1:n*(1.0)), cn, 
                cex =0.4, adj = c(0,1), facing = "bending.inside")
  }
}, bg.border = NA)
circos.clear()
lgd = Legend(title = "Regulon activity", col_fun = col_fun1, direction = "horizontal", 
             title_position = "topcenter", at = c( -2, 0, 2))
draw(lgd, x = unit(1, "cm"), y = unit(1, "cm"), just = c("left", "bottom"))
dev.off()

# ── 5. Fig 6D-E: Differential Regulon-pseudotime correlation (BMP4 vs control separately) ───────


DefaultAssay(GABAtrac) <- "AUC"
expr_scaled <- ScaleData(GABAtrac,
                         features = rownames(GABAtrac)[rowSums(GABAtrac@assays$RNA@data) > 1])
cor_TF <- cor(t(expr_scaled@assays$RNA@scale.data), subRGoRG$pseudotime_m3, method="pearson")
cor_TF <- tibble::rownames_to_column(as.data.frame(cor_TF), "gene") %>%
  rename(cor=1) %>% arrange(desc(cor))
write.csv(cor_TF, file.path(path_output,"GABAtrac_pseudotime_cor_TFs.csv"), row.names=FALSE)


Idents(GABAtrac) <- GABAtrac$orig.ident
TF_markers <- FindMarkers(GABAtrac, ident.1="BMP4", ident.2="control",
                          only.pos=FALSE, logfc.threshold=0.1)
TF_markers$combined_score <- TF_markers$avg_log2FC * (-log10(TF_markers$p_val))
TF_markers <- tibble::rownames_to_column(TF_markers, "GENE")
write.csv(TF_markers, file.path(path_output,"GABAtrac_TF_BMP4vsCtrl.csv"), row.names=FALSE)

driver_TFs <- intersect(
  cor_tf$gene[cor_tf$cor > 0.4 & cor_tf$pvalue<0.05],
  TF_markers$GENE[TF_markers$avg_log2FC > 0.5 & TF_markers$FDR<0.05])


# ── 6. Fig 6F: TF network building ───────
#see TF_gene_network_building.r, run with NPC_to_GABA trajectory driver_TF and driver_genes.
