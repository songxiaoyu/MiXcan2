suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  "scripts"
}
source(file.path(script_dir, "config.R"))

if (!file.exists(weights_file)) {
  stop("Weights file does not exist: ", weights_file, call. = FALSE)
}

weights <- fread(weights_file)
required <- c(
  "cell_type", "yName", "xNameMatrix", "type",
  "weight_cell_1", "weight_cell_2"
)
missing_columns <- setdiff(required, names(weights))
if (length(missing_columns)) {
  stop("Weights file is missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

if (weights[, anyNA(.SD), .SDcols = required]) {
  stop("Weights file contains missing values in required columns.", call. = FALSE)
}

cell_map <- c(
  epithelial = "epithelial",
  fibroblasts = "fibroblast",
  adipocytes = "adipose",
  predixcan = "predixcanlike"
)
unknown_cells <- setdiff(unique(weights$cell_type), names(cell_map))
if (length(unknown_cells)) {
  stop("Unknown cell_type values: ", paste(unknown_cells, collapse = ", "), call. = FALSE)
}

weights[, model_type := unname(cell_map[cell_type])]
weights[, gene_id := yName]
weights[, MarkerName_hg38 := sub("^chr", "", xNameMatrix, ignore.case = TRUE)]
weights[, MarkerName_hg38 := sub("_b38$", "", MarkerName_hg38)]
weights[, MarkerName_hg38 := gsub("_", ":", MarkerName_hg38, fixed = TRUE)]

valid_marker <- grepl("^[^:]+:[0-9]+:[^:]+:[^:]+$", weights$MarkerName_hg38)
if (!all(valid_marker)) {
  stop("Could not parse xNameMatrix for ", sum(!valid_marker), " rows.", call. = FALSE)
}

weights[, chr := tstrsplit(MarkerName_hg38, ":", fixed = TRUE)[[1]]]
weights[, position := as.integer(tstrsplit(MarkerName_hg38, ":", fixed = TRUE)[[2]])]

duplicate_n <- weights[, sum(duplicated(.SD)), .SDcols = c(
  "model_type", "gene_id", "MarkerName_hg38"
)]
if (duplicate_n) {
  stop("Duplicate model/gene/variant rows found: ", duplicate_n, call. = FALSE)
}

setcolorder(weights, c(
  "model_type", "gene_id", "MarkerName_hg38", "chr", "position",
  "type", "weight_cell_1", "weight_cell_2", "cell_type", "yName",
  "xNameMatrix"
))
fwrite(weights, normalized_weights_file)

manifest <- weights[, .(
  input_weight_rows = .N,
  input_weight_snps = uniqueN(MarkerName_hg38),
  model_type_selected = paste(sort(unique(type)), collapse = ";")
), by = .(model_type, gene_id, chr)]
fwrite(manifest, weights_manifest_file)

cat("Validated weights:", weights_file, "\n")
cat("Rows:", nrow(weights), "\n")
cat("Genes:", uniqueN(weights$gene_id), "\n")
cat("Variants:", uniqueN(weights$MarkerName_hg38), "\n")
print(weights[, .(rows = .N, genes = uniqueN(gene_id)), by = model_type])
