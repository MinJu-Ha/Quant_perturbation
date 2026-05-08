suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
  library(ggplot2)
  library(Matrix)
})

config <- yaml::read_yaml("configs/config.yaml")

obj <- readRDS(config$paths$seurat_rds)

get_assay_matrix <- function(object, assay, layer) {
  tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(layer_error) {
      GetAssayData(object, assay = assay, slot = layer)
    }
  )
}

pick_features <- function(available_features, requested_features) {
  exact <- requested_features[requested_features %in% available_features]
  lower_map <- stats::setNames(available_features, tolower(available_features))
  case_insensitive <- unname(lower_map[tolower(requested_features)])
  unique(c(exact, case_insensitive[!is.na(case_insensitive)]))
}

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

chip <- read.csv(config$paths$chip_targets_csv, check.names = FALSE)
chip_targets <- unique(gsub("\ufeff", "", as.character(chip[[1]]), fixed = TRUE))
chip_targets <- chip_targets[!is.na(chip_targets) & chip_targets != ""]
min_target_genes <- as.integer(config$features$min_target_genes)

obj$cell_id <- colnames(obj)
meta <- obj@meta.data
meta$cell_id <- rownames(meta)
write.csv(meta, config$paths$exported_meta_csv, row.names = TRUE, quote = FALSE)

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- ScaleData(obj, features = rownames(obj))
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)

cat("BCL6 in RNA variable features after re-normalization:", target_gene %in% VariableFeatures(obj[["RNA"]]), "\n")

expr <- get_assay_matrix(obj, assay = "RNA", layer = "scale.data")

features_to_export <- pick_features(rownames(expr), unique(c(target_candidates, chip_targets)))
if (length(features_to_export) < min_target_genes + 1) {
  stop(paste("Too few BCL6-related features available for export:", length(features_to_export)))
}
expr <- expr[features_to_export, , drop = FALSE]

write.csv(as.matrix(expr), config$paths$exported_expression_csv, quote = FALSE)

counts_path <- config$paths$exported_counts_csv
if (!is.null(counts_path) && nzchar(counts_path)) {
  counts_mat <- get_assay_matrix(obj, assay = "RNA", layer = "counts")
  counts_features <- pick_features(rownames(counts_mat), unique(c(target_candidates, chip_targets)))
  if (length(counts_features) >= min_target_genes + 1) {
    counts_sub <- counts_mat[counts_features, , drop = FALSE]
    counts_dense <- as.matrix(counts_sub)
    storage.mode(counts_dense) <- "integer"
    write.csv(counts_dense, counts_path, quote = FALSE)
    cat("Wrote RNA raw counts to:", counts_path, "(", nrow(counts_dense), "genes x", ncol(counts_dense), "cells )\n")
  } else {
    cat("Skipping raw counts export — too few BCL6-related features in counts slot.\n")
  }
}

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
