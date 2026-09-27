suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  "scripts"
}
source(file.path(script_dir, "config.R"))

input_files <- file.path(results_dir, c(
  "md_results.csv", "bcac_overall_results.csv", "bcac_subtype_results.csv"
))
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing association results: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

output_columns <- c(
  "analysis", "outcome", "model_type", "gene_id", "chr",
  "training_model_type", "input_snp_num", "matched_snp_num",
  "Z_cell_1", "p_cell_1", "Z_cell_2", "p_cell_2", "p_join"
)
results <- rbindlist(lapply(input_files, fread), fill = TRUE)
missing_columns <- setdiff(output_columns, names(results))
if (length(missing_columns)) {
  stop("Missing result columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

results <- results[, ..output_columns]
setorder(results, analysis, outcome, model_type, p_join, na.last = TRUE)
output_file <- file.path(results_dir, "all_results.csv")
fwrite(results, output_file)

cat("Combined result table written to:", output_file, "\n")
print(results[, .N, by = .(analysis, outcome)])
