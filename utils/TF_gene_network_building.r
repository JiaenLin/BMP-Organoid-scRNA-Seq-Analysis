# ==============================================================================
# TF–Gene Regulatory Network: Gene–Regulon Binary Matrix
# Input:  regulons — named list from pySCENIC (TF → target genes)
# Output: gene_reg_map — binary data.frame (genes × TF regulons)
# Author: Jiaen Lin
# Date: 08-May-2026
# ==============================================================================

library(Matrix)

all_genes     <- unique(unlist(regulons))
regulon_names <- names(regulons)

gene_regulon_mat <- Matrix(
  0,
  nrow     = length(all_genes),
  ncol     = length(regulon_names),
  dimnames = list(all_genes, regulon_names),
  sparse   = FALSE
)

for (i in seq_along(regulon_names)) {
  matched_rows <- match(regulons[[regulon_names[i]]], all_genes)
  gene_regulon_mat[matched_rows, i] <- 1
}

gene_reg_map <- as.data.frame(gene_regulon_mat)


# ==============================================================================
# QUERY: Subset to driver genes × driver TF regulons
# ==============================================================================
### get driver_TFs list from Fig3_RGoRG_trajectory.R driver_TF calculation.
driver_TFs        <- c("SOX2", "CEBPG", "SOX11", "PAX6", "NFIB",
                       "BCL11A", "KLF7", "TCF7L1", "ZEB1", "ZBTB20")
driver_TFs_scenic <- paste0(driver_TFs, "(+)")  # match pySCENIC naming
###get DEG list from DEG analysis
network_query <- gene_reg_map[
  intersect(DEG, rownames(gene_reg_map)),
  intersect(driver_TFs_scenic, colnames(gene_reg_map))
]
network_query <- network_query[rowSums(network_query) > 0, ]
network_query <- network_query[, order(colSums(network_query), decreasing = TRUE)]
