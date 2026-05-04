suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(ggplot2)
})

config <- yaml::read_yaml("configs/config.yaml")

obj <- readRDS(config$paths$seurat_rds)

processed_dir <- config$paths$processed_dir
tables_dir <- config$paths$tables_dir
figures_dir <- config$paths$figures_dir

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

condition_col <- config$metadata$condition_column
if (!condition_col %in% colnames(obj@meta.data)) {
  stop(paste("Missing condition column:", condition_col))
}

target_candidates <- unlist(config$perturbation$target_gene_candidates)
target_gene <- target_candidates[target_candidates %in% rownames(obj)][1]
if (is.na(target_gene)) {
  stop("Could not find BCL6/Bcl6/bcl6 in rownames(obj).")
}

obj$cell_id <- colnames(obj)
meta <- obj@meta.data
meta$cell_id <- rownames(meta)
write.csv(meta, config$paths$exported_meta_csv, row.names = TRUE, quote = FALSE)

if ("SCT" %in% Assays(obj)) {
  DefaultAssay(obj) <- "SCT"
  expr <- GetAssayData(obj, assay = "SCT", slot = "scale.data")
  if (nrow(expr) == 0 || ncol(expr) == 0) {
    warning("SCT scale.data is empty. Falling back to RNA data slot.")
    DefaultAssay(obj) <- "RNA"
    expr <- GetAssayData(obj, assay = "RNA", slot = "data")
  }
} else {
  DefaultAssay(obj) <- "RNA"
  expr <- GetAssayData(obj, assay = "RNA", slot = "data")
}

write.csv(as.matrix(expr), config$paths$exported_expression_csv, quote = FALSE)

bcl6_df <- FetchData(obj, vars = c(target_gene, condition_col))
colnames(bcl6_df)[1] <- "BCL6_expression"
bcl6_df$cell_id <- rownames(bcl6_df)
write.csv(bcl6_df, file.path(tables_dir, "bcl6_expression_by_cell.csv"), row.names = FALSE)

summary_df <- aggregate(
  BCL6_expression ~ get(condition_col),
  data = bcl6_df,
  FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x), n = length(x))
)
write.csv(summary_df, file.path(tables_dir, "bcl6_expression_by_condition.csv"), row.names = FALSE)

p <- VlnPlot(obj, features = target_gene, group.by = condition_col, pt.size = 0) +
  ggtitle(paste(target_gene, "expression by condition"))

ggsave(file.path(figures_dir, "bcl6_expression_by_condition.png"), p, width = 7, height = 5, dpi = 300)

if ("umap" %in% names(obj@reductions)) {
  p2 <- FeaturePlot(obj, features = target_gene, reduction = "umap") +
    ggtitle(paste(target_gene, "expression on UMAP"))
  ggsave(file.path(figures_dir, "bcl6_expression_umap.png"), p2, width = 6, height = 5, dpi = 300)
}

cat("\nExport complete.\n")
cat("Target gene used:", target_gene, "\n")
