# =============================================================================
# Main analysis: GSE298023
# Includes module scores, overlap, marker enrichment, GO,
# sensitivity analysis, and fibroblast regional analysis
# =============================================================================

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

source("gene_lists.R")

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
seu <- readRDS("GSE298023_HH_dataset_seurat_object.rds")

dcm_genes <- intersect(dcm_genes, rownames(seu))
cad_genes <- intersect(cad_genes, rownames(seu))

# -----------------------------------------------------------------------------
# 2. Module scores
# -----------------------------------------------------------------------------
seu <- AddModuleScore(seu, features = list(dcm_genes), name = "DCM_exp")
seu <- AddModuleScore(seu, features = list(cad_genes), name = "CAD_exp")

score_summary <- seu@meta.data %>%
  group_by(ident) %>%
  summarise(
    DCM = mean(DCM_exp1, na.rm = TRUE),
    CAD = mean(CAD_exp1, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(DCM))

write.csv(score_summary, "results_celltype_module_scores.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. Gene set overlap
# -----------------------------------------------------------------------------
overlap_genes <- intersect(dcm_genes, cad_genes)
N <- nrow(seu)
p_overlap <- phyper(length(overlap_genes), length(dcm_genes), 
                    N - length(dcm_genes), length(cad_genes), lower.tail = TRUE)

# -----------------------------------------------------------------------------
# 4. Sensitivity analysis (exclude highly expressed genes)
# -----------------------------------------------------------------------------
# Remove genes with very high expression to test robustness
avg_expr <- AverageExpression(seu, features = c(dcm_genes, cad_genes), 
                              group.by = "ident")$RNA

high_expr_genes <- names(which(rowMeans(avg_expr) > quantile(rowMeans(avg_expr), 0.9)))

dcm_sens <- setdiff(dcm_genes, high_expr_genes)
cad_sens <- setdiff(cad_genes, high_expr_genes)

seu <- AddModuleScore(seu, features = list(dcm_sens), name = "DCM_sens")
seu <- AddModuleScore(seu, features = list(cad_sens), name = "CAD_sens")

sens_summary <- seu@meta.data %>%
  group_by(ident) %>%
  summarise(
    DCM_sens = mean(DCM_sens1, na.rm = TRUE),
    CAD_sens = mean(CAD_sens1, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(DCM_sens))

write.csv(sens_summary, "results_sensitivity_scores.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Fibroblast regional analysis
# -----------------------------------------------------------------------------
fb <- subset(seu, idents = "Fibroblasts")

fb_region <- fb@meta.data %>%
  group_by(region) %>%
  summarise(
    DCM = mean(DCM_exp1, na.rm = TRUE),
    CAD = mean(CAD_exp1, na.rm = TRUE),
    n = n()
  )

write.csv(fb_region, "results_fibroblast_region_scores.csv", row.names = FALSE)

# Optional: fibroblast state (if available)
# if ("fb_state" %in% colnames(fb@meta.data)) {
#   fb_state_summary <- fb@meta.data %>%
#     group_by(fb_state) %>%
#     summarise(DCM = mean(DCM_exp1), CAD = mean(CAD_exp1), n = n())
#   write.csv(fb_state_summary, "results_fibroblast_state_scores.csv", row.names = FALSE)
# }

message("Main analysis + sensitivity + fibroblast analysis completed.")