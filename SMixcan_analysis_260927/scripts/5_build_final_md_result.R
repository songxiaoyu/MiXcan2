suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  getwd()
}
workflow_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

new_result_file <- file.path(
  workflow_dir, "data", "derived",
  "md_results_corrected_reg0p001_all_models.csv"
)
new_result_file <- Sys.getenv(
  "SMIXCAN_CORRECTED_MD_FILE",
  unset = new_result_file
)
result3_file <- Sys.getenv(
  "SMIXCAN_RESULT3_MD_FILE",
  unset = file.path(
    workflow_dir, "data", "reference",
    "result3_md_wrong_direction_benchmark.csv"
  )
)
output_dir <- Sys.getenv(
  "SMIXCAN_FINAL_MD_DIR",
  unset = file.path(workflow_dir, "results", "final_md")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

missing_files <- c(new_result_file, result3_file)[
  !file.exists(c(new_result_file, result3_file))
]
if (length(missing_files)) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}

new_result <- fread(new_result_file)
result3 <- fread(result3_file)

setnames(
  result3,
  old = c(
    "phenotype", "celltype", "type", "filtered_snp_num",
    "Z_1_join", "p_1_join", "Z_2_join", "p_2_join"
  ),
  new = c(
    "outcome", "model_type", "training_model_type", "matched_snp_num",
    "Z_cell_1", "p_cell_1", "Z_cell_2", "p_cell_2"
  )
)

key <- c("outcome", "model_type", "gene_id", "training_model_type")
required_new <- c(
  "analysis", key, "chr", "input_snp_num", "matched_snp_num",
  "Z_cell_1", "p_cell_1", "Z_cell_2", "p_cell_2", "p_join"
)
required_old <- c(
  key, "matched_snp_num", "Z_cell_1", "p_cell_1",
  "Z_cell_2", "p_cell_2", "p_join"
)
if (length(setdiff(required_new, names(new_result)))) {
  stop("The corrected 0.001 result is missing required columns.")
}
if (length(setdiff(required_old, names(result3)))) {
  stop("The Result3 table is missing required columns.")
}

# PrediXcan-like rows remain a separate comparator and are not included in the
# collaborator-facing S-MiXcan table.
final <- new_result[model_type != "predixcanlike"]
result3_models <- unique(result3[, ..key])
result3_models[, in_result3 := TRUE]
result3_genes <- unique(result3$gene_id)

final <- merge(final, result3_models, by = key, all.x = TRUE)
final[, model_status := fcase(
  !is.na(in_result3), "result3_overlap",
  !gene_id %chin% result3_genes, "new_gene",
  default = "new_model_existing_gene"
)]
final[, in_result3 := NULL]

column_order <- c(
  "analysis", "outcome", "model_type", "gene_id", "chr",
  "training_model_type", "model_status", "input_snp_num",
  "matched_snp_num", "Z_cell_1", "p_cell_1", "Z_cell_2",
  "p_cell_2", "p_join"
)
setcolorder(final, column_order)
setorder(final, model_status, model_type, gene_id, outcome)

comparison <- merge(
  final[model_status == "result3_overlap"],
  result3[, ..required_old],
  by = key,
  suffixes = c("_corrected", "_result3")
)
comparison[, `:=`(
  matched_snp_diff = matched_snp_num_corrected - matched_snp_num_result3,
  Z_cell_1_direction_check = Z_cell_1_corrected + Z_cell_1_result3,
  Z_cell_2_direction_check = Z_cell_2_corrected + Z_cell_2_result3,
  p_cell_1_diff = p_cell_1_corrected - p_cell_1_result3,
  p_cell_2_diff = p_cell_2_corrected - p_cell_2_result3,
  p_join_diff = p_join_corrected - p_join_result3
)]

summary <- data.table(
  metric = c(
    "final_rows",
    "final_genes",
    "final_models",
    "result3_overlap_rows",
    "result3_overlap_genes",
    "result3_overlap_models",
    "new_gene_rows",
    "new_genes",
    "new_gene_models",
    "exact_matched_snp_rows",
    "p_join_equal_within_1e_12",
    "direction_reversed_cell_1_within_1e_12",
    "direction_reversed_cell_2_within_1e_12",
    "max_abs_p_join_diff",
    "max_abs_direction_check_cell_1",
    "max_abs_direction_check_cell_2"
  ),
  value = c(
    nrow(final),
    uniqueN(final$gene_id),
    uniqueN(final[, paste(model_type, gene_id, training_model_type)]),
    final[model_status == "result3_overlap", .N],
    uniqueN(final[model_status == "result3_overlap", gene_id]),
    uniqueN(final[
      model_status == "result3_overlap",
      paste(model_type, gene_id, training_model_type)
    ]),
    final[model_status == "new_gene", .N],
    uniqueN(final[model_status == "new_gene", gene_id]),
    uniqueN(final[
      model_status == "new_gene",
      paste(model_type, gene_id, training_model_type)
    ]),
    sum(comparison$matched_snp_diff == 0),
    sum(abs(comparison$p_join_diff) < 1e-12),
    sum(abs(comparison$Z_cell_1_direction_check) < 1e-12),
    sum(abs(comparison$Z_cell_2_direction_check) < 1e-12),
    max(abs(comparison$p_join_diff)),
    max(abs(comparison$Z_cell_1_direction_check)),
    max(abs(comparison$Z_cell_2_direction_check))
  )
)

stopifnot(
  nrow(final) == 177L,
  uniqueN(final$gene_id) == 21L,
  final[model_status == "result3_overlap", .N] == 123L,
  final[model_status == "new_gene", .N] == 54L,
  !anyNA(final$p_join),
  all(comparison$matched_snp_diff == 0),
  all(abs(comparison$p_join_diff) < 1e-12),
  all(abs(comparison$Z_cell_1_direction_check) < 1e-12),
  all(abs(comparison$Z_cell_2_direction_check) < 1e-12)
)

fwrite(
  final,
  file.path(output_dir, "smixcan_result3_method_direction_corrected.csv")
)
fwrite(
  comparison,
  file.path(output_dir, "result3_overlap_direction_validation.csv")
)
fwrite(summary, file.path(output_dir, "validation_summary.csv"))

cat("Result3-compatible direction-corrected result written to:", output_dir, "\n")
print(summary)
