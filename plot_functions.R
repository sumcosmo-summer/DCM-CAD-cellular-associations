# =============================================================================
# Plotting functions
# =============================================================================

library(Seurat)
library(ggplot2)

# Cell type module score violin plot
plot_celltype_scores <- function(seu_obj, dcm_feature = "DCM_exp1", 
                                 cad_feature = "CAD_exp1") {
  p <- VlnPlot(seu_obj, features = c(dcm_feature, cad_feature), 
               group.by = "ident", pt.size = 0) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(p)
}

# Save main figure
save_main_figure <- function(seu_obj, filename = "figures/Main_VlnPlot.pdf") {
  dir.create("figures", showWarnings = FALSE)
  p <- plot_celltype_scores(seu_obj)
  ggsave(filename, p, width = 10, height = 6)
}

# Fibroblast region violin plot
plot_fb_region <- function(fb_obj) {
  p <- VlnPlot(fb_obj, features = c("DCM_exp1", "CAD_exp1"), 
               group.by = "region", pt.size = 0)
  return(p)
}

# GO dotplot helper
save_go_dotplot <- function(go_result, title, filename) {
  if (!is.null(go_result)) {
    p <- dotplot(go_result, showCategory = 10) + ggtitle(title)
    ggsave(filename, p, width = 8, height = 6)
  }
}