# ==============================================================================
# Figure 2: scRNA-seq Clustering, Cell Type Annotation, and Regulon Analysis
# of BMP4-treated Human Cerebral Organoids
#
# Corresponds to: Figure 2A-G
# Input: raw 10X Chromium count matrices (Control, BMP4)
# Output: Annotated Seurat object (bmp4recluster2.RDS), UMAP plots,
#         cell composition bar charts, regulon heatmap
#Author: Jiaen Lin
#Date: 08-May-2026 
# ==============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)
library(scCustomize)
library(SCopeLoomR)
library(SCENIC)
library(AUCell)
library(ComplexHeatmap)

path_input   <- "data/raw_counts/"
path_output  <- "results/Fig2/"
path_deposit <- "data/processed/"
dir.create(path_output,  showWarnings = FALSE, recursive = TRUE)
dir.create(path_deposit, showWarnings = FALSE, recursive = TRUE)

# ── 1. Load 10X data ──────────────────────────────────────────────────────────
ctrl.data <- Read10X(data.dir = file.path(path_input, "Ctrl/filtered_feature_bc_matrix"))
bmp4.data <- Read10X(data.dir = file.path(path_input, "BMP4/filtered_feature_bc_matrix"))

ctrl <- CreateSeuratObject(counts = ctrl.data, project = "control",  min.cells = 3, min.features = 200)
bmp4 <- CreateSeuratObject(counts = bmp4.data, project = "BMP4",    min.cells = 3, min.features = 200)

# ── 2. QC ─────────────────────────────────────────────────────────────────────
for (obj_name in c("ctrl", "bmp4")) {
  obj <- get(obj_name)
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$percent.hb <- PercentageFeatureSet(obj, pattern = "^HBA|HBB")
  assign(obj_name, obj)
}
ctrl <- subset(ctrl, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 20)
bmp4 <- subset(bmp4, subset = nFeature_RNA > 200 & nFeature_RNA < 12000 & percent.mt < 20)

# ── 3. Merge, normalise, cluster ──────────────────────────────────────────────
sample.combined <- merge(ctrl, bmp4, project = "combined", merge.data = TRUE)
sample.combined <- NormalizeData(sample.combined, normalization.method = "LogNormalize", scale.factor = 1e4)
sample.combined <- FindVariableFeatures(sample.combined, selection.method = "vst", nfeatures = 5000)
sample.combined <- ScaleData(sample.combined, features = rownames(sample.combined))
sample.combined <- RunPCA(sample.combined, features = VariableFeatures(sample.combined), seed.use = 42)
sample.combined <- RunUMAP(sample.combined, dims = 1:30)
sample.combined <- FindNeighbors(sample.combined, dims = 1:30)
sample.combined <- FindClusters(sample.combined, resolution = 0.8)

# ── 4. Cell type annotation ───────────────────────────────────────────────────
marker.list <- list(
  RG          = c("SOX2","PAX6","HES1","VIM"),
  oRG         = c("HOPX","PTPRZ1","FAM107A","TNC","MOXD1","SOX2","PTN","LIFR"),
  IPC         = c("EOMES","NEUROD1","NEUROD4","NEUROG2","ASCL1"),
  Mature.N    = c("RBFOX3","NEUROD2","SYP","DLG4","STMN2"),
  Glu.N       = c("GRIN2B","GRIN1","SLC17A7","SLC17A6"),
  GABA.N      = c("GAD1","DLX1","DLX5","ERBB4","GAD2"),
  Astrocytes  = c("GFAP","S100B","AQP4","SLC1A2","ALDH1L1"),
  Cajal       = c("RELN","LHX1"),
  Mesenchymal = c("DCN","LUM","COL1A2","COL5A1","COL3A1")
)
sample.combined <- AddModuleScore(sample.combined, features = marker.list, ctrl = 10, name = "marker.")

Idents(sample.combined) <- "RNA_snn_res.0.8"
#manual cell type annotation based on marker expressions
sample.combined <- RenameIdents(sample.combined,
  "0"="RG","2"="RG","7"="RG","9"="RG","13"="RG","19"="RG",
  "3"="IPC","11"="IPC","15"="IPC",
  "4"="Mature.N","8"="Mature.N","12"="Mature.N",
  "5"="Glu.N","10"="Glu.N","14"="Glu.N",
  "6"="Immature.N","17"="Immature.N",
  "16"="Astrocytes","18"="GABA.N","20"="Cajal","1"="oRG"
)
sample.combined$annotation <- Idents(sample.combined)
levels(sample.combined) <- c("RG","oRG","IPC","Immature.N","Mature.N",
                              "Glu.N","GABA.N","Astrocytes","Cajal","Mesenchymal")

saveRDS(sample.combined, file.path(path_deposit, "bmp4.RDS"))

sample.combined<-readRDS("bmp4.RDS")
# ── 5. Fig 2A–C: UMAP + marker plots ─────────────────────────────────────────
color_palette <- c("#F68282","#B95FBB","firebrick","#31C53F","yellow3",
                   "#1FA195","darkgray","tomato2","#AC8F14","orange")

p_umap_all   <- DimPlot(sample.combined, label = TRUE, pt.size = 0.5, cols = color_palette)
p_umap_split <- DimPlot(sample.combined, split.by = "orig.ident", label = TRUE, pt.size = 0.5, cols = color_palette)
p_markers    <- FeaturePlotScCustom(sample.combined,
                  features = c("SOX2","HOPX","EOMES","GRIN2B","DLX5","GFAP"),
                  num_columns = 3, na_cutoff = 1, colors_use = viridis::magma(100, direction = -1))
p_dotplot    <- DotPlot(sample.combined,
                  features = c("SOX2","GLI3","PAX6","HOPX","GFAP","AQP4","S100B","EOMES",
                               "NEUROD4","DCX","STMN2","RBFOX3","SYP","GRIN2B","GRIN1",
                               "RELN","LHX1","DLX1","DLX5","ERBB4","GAD2","LUM","COL1A2"),
                  cols = c("cyan","red"), col.min = -2, col.max = 2, dot.scale = 6) +
                RotatedAxis() + coord_flip()

ggsave(file.path(path_output, "Fig2A_UMAP_all.pdf"),     p_umap_all,   width = 8,  height = 6)
ggsave(file.path(path_output, "Fig2B_dotplot.pdf"),      p_dotplot,    width = 10, height = 8)
ggsave(file.path(path_output, "Fig2C_markers_UMAP.pdf"), p_markers,    width = 12, height = 8)
ggsave(file.path(path_output, "Fig2D_UMAP_split.pdf"),   p_umap_split, width = 14, height = 6)

# ── 6. Fig 2E–F: Cell composition ────────────────────────────────────────────
sample.combined$celltype   <- Idents(sample.combined)
sample.combined$orig.ident <- factor(sample.combined$orig.ident, levels = c("control","BMP4"))

prop_data <- sample.combined@meta.data %>%
  group_by(orig.ident, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percent = n / sum(n) * 100)

p_bar <- ggplot(prop_data, aes(x = orig.ident, y = percent, fill = celltype)) +
  geom_col(width = 0.7) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() +
  scale_fill_manual(values = color_palette) +
  coord_flip()
ggsave(file.path(path_output, "Fig2E_composition_bar.pdf"), p_bar, width = 5, height = 8)

p_bar2<-SCpubr::do_BarPlot(sample.combined, 
                   group.by = "orig.ident",
                   split.by = "celltype",
                   position = "fill",
                   flip = FALSE,
                   add.n = TRUE,
                   add.n.size = 3,
                   return_data = TRUE,colors.use=c("control"='#1f77b4',"BMP4"= '#fee8c8'))
ggsave(file.path(path_output, "Fig2F_composition_bar.pdf"), p_bar2, width = 8, height = 5)

# ── 7. Fig 2G: Regulon heatmap ────────────────────────────────────────────────
loom <- open_loom(file.path("results/SCENIC","bmprecluster2_pyscenic.loom"))
regulon_AUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
close_loom(loom)

AUCmat <- AUCell::getAUC(regulon_AUC)
rownames(AUCmat) <- gsub("\\(.*\\)", "", rownames(AUCmat))
AUCassay <- ScaleData(CreateAssayObject(counts = AUCmat))
sample.combined[["AUC"]] <- AUCassay

rss <- SCENIC::calcRSS(AUC = AUCmat, cellAnnotation = sample.combined$celltype)
top_labels <- character(0)
for (cn in colnames(rss))
  top_labels <- unique(c(top_labels, rownames(rss)[order(rss[, cn], decreasing = TRUE)[1:10]]))

DefaultAssay(sample.combined) <- "AUC"
sample.combined$class<-Idents(sample.combined)
sample.combined$treat.class<-paste0(sample.combined$class,".",sample.combined$orig.ident)
seurat_bulk<-AverageExpression(sample.combined,group.by = 'treat.class',return.seurat = TRUE)
seurat_bulk$group<-str_split_i(seurat_bulk$treat.class, "\\.", -1) 
seurat_bulk$celltype<-gsub("\\.[^.]+$", "", seurat_bulk$treat.class)
splitlevel<-c()
list<-c('RG','oRG','IPC','Immature.N','Mature.N','Glu.N','Cajal',"GABA.N","Astrocytes")
seurat_bulk$group <- factor(seurat_bulk$group, levels = c("control", "BMP4"))
seurat_bulk$celltype <- factor(seurat_bulk$celltype, levels = c('RG','oRG','NPC','Immature.N','Mature.N','Glu.N','Cajal',"Astrocytes"))
for (id in list){
  for (time in rev(unique(seurat_bulk$group))){
    splitlevel<-c(splitlevel,paste0(id,'.',time))
  }
}
splitlevel<-intersect(splitlevel,seurat_bulk$treat.class)
table(bmp4v2$celltype)
labels<-unique(labels)
ggData <- GetAssayData(object = seurat_bulk, slot = "data")
ggData <- ggData[labels,]
ggData <- t(scale(t(ggData)))
ggData[ggData > 2] <- 2; ggData[ggData < -2] <- -2
ggData<-ggData[,splitlevel]
tmp <- data.frame(group = rev(seurat_bulk$group),
                  celltype = rev(seurat_bulk$celltype),row.names = splitlevel)
colormap_l1 = c("control"='#fee8c8',"BMP4"='#1f77b4',"RG"='#31a354',"oRG"='#31a354',
                "NPC"='#31a354',"Immature.N"='#31a354',"Mature.N"='#31a354',"Glu.N"='#31a354',"Cajal"='#31a354',"Astrocytes"='#31a354')
ann_colors = list(
  group = c("control"='#4292c6',"BMP4"='#fcbba1'),
  celltype = c("RG"='#F68282',"oRG"='#B95FBB',
               "NPC"='firebrick',"Immature.N"='#31C53F',"Mature.N"='yellow3',"Glu.N"="#1FA195","Cajal"='#31a354',"Astrocytes"='#AC8F14')
)
p_regulon <- pheatmap(as.matrix(ggData),annotation_colors= ann_colors, cluster_cols = FALSE, show_colnames =TRUE,cluster_rows=TRUE,
         annotation_col = tmp, fontsize_row = 10, fontsize_col =15,angle_col = c("45"),fontsize=15,
         width = 20, height = 8,
         cellwidth = 40, cellheight = 10)

ggsave(file.path(path_output, "Fig2G_regulon_heatmap.pdf"), p_regulon, width = 14, height = 10)










