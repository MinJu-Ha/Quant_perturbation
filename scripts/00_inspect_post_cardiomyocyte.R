suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

config <- yaml::read_yaml("configs/config.yaml")
obj <- readRDS(config$paths$seurat_rds)

cat("\n===== Object summary =====\n")
print(obj)

cat("\n===== Assays =====\n")
print(Assays(obj))

cat("\n===== Default assay =====\n")
print(DefaultAssay(obj))

cat("\n===== Metadata columns =====\n")
print(colnames(obj@meta.data))

condition_col <- config$metadata$condition_column
if (condition_col %in% colnames(obj@meta.data)) {
  cat("\n===== Condition counts =====\n")
  print(table(obj@meta.data[[condition_col]], useNA = "ifany"))
} else {
  cat("\nCondition column not found:", condition_col, "\n")
}

subtype_col <- config$metadata$subtype_column
if (subtype_col %in% colnames(obj@meta.data)) {
  cat("\n===== CM subtype counts =====\n")
  print(table(obj@meta.data[[subtype_col]], useNA = "ifany"))
}

cat("\n===== Target gene name check =====\n")
target_candidates <- unlist(config$perturbation$target_gene_candidates)
for (g in target_candidates) {
  cat(g, "in rownames:", g %in% rownames(obj), "\n")
}

cat("\nGenes matching BCL6 ignore case:\n")
print(rownames(obj)[grepl("^bcl6$", rownames(obj), ignore.case = TRUE)])

cat("\n===== Assay slot check =====\n")
for (assay in Assays(obj)) {
  cat("\nAssay:", assay, "\n")
  counts_dim <- tryCatch(dim(GetAssayData(obj, assay = assay, slot = "counts")), error = function(e) NA)
  data_dim <- tryCatch(dim(GetAssayData(obj, assay = assay, slot = "data")), error = function(e) NA)
  scale_dim <- tryCatch(dim(GetAssayData(obj, assay = assay, slot = "scale.data")), error = function(e) NA)
  cat("counts:", paste(counts_dim, collapse = " x "), "\n")
  cat("data:", paste(data_dim, collapse = " x "), "\n")
  cat("scale.data:", paste(scale_dim, collapse = " x "), "\n")
}

cat("\nDone. Send this output if you want me to verify the setup.\n")
