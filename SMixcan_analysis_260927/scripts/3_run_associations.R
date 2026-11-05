suppressPackageStartupMessages({
  library(data.table)
  library(SMiXcan)
})

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  "scripts"
}
source(file.path(script_dir, "config.R"))

if (!file.exists(normalized_weights_file)) stop("Run step 1 first.", call. = FALSE)

normalize_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

is_snp <- function(a) {
  a <- toupper(as.character(a))
  nchar(a) == 1L & a %chin% c("A", "C", "G", "T")
}

complement <- function(a) {
  map <- c(A = "T", T = "A", C = "G", G = "C")
  unname(map[toupper(as.character(a))])
}

harmonize_to_reference <- function(gwas, reference_markers, chr_col, pos_col,
                                   a1_col, a2_col, statistic_columns,
                                   p_col = NULL,
                                   flip_columns = statistic_columns) {
  g <- copy(gwas)
  g[, chr_h := normalize_chr(get(chr_col))]
  g[, pos_h := as.character(get(pos_col))]
  g[, a1_h := toupper(as.character(get(a1_col)))]
  g[, a2_h := toupper(as.character(get(a2_col)))]
  for (column in statistic_columns) {
    g[, (column) := suppressWarnings(as.numeric(get(column)))]
  }
  if (!is.null(p_col) && p_col %chin% names(g)) {
    g[, p_sort := suppressWarnings(as.numeric(get(p_col)))]
  } else {
    g[, p_sort := seq_len(.N)]
  }

  g[, id_same := paste(chr_h, pos_h, a1_h, a2_h, sep = ":")]
  g[, id_swap := paste(chr_h, pos_h, a2_h, a1_h, sep = ":")]
  complement_ok <- is_snp(g$a1_h) & is_snp(g$a2_h)
  g[, `:=`(id_comp_same = NA_character_, id_comp_swap = NA_character_)]
  if (any(complement_ok)) {
    ca1 <- complement(g$a1_h[complement_ok])
    ca2 <- complement(g$a2_h[complement_ok])
    g$id_comp_same[complement_ok] <- paste(
      g$chr_h[complement_ok], g$pos_h[complement_ok], ca1, ca2, sep = ":"
    )
    g$id_comp_swap[complement_ok] <- paste(
      g$chr_h[complement_ok], g$pos_h[complement_ok], ca2, ca1, sep = ":"
    )
  }

  reference <- data.table(MarkerName_hg38 = reference_markers)
  setkey(reference, MarkerName_hg38)

  map_one <- function(id_column, flip) {
    selected <- g[!is.na(get(id_column)), c(
      list(MarkerName_hg38 = get(id_column), p_sort = p_sort),
      mget(statistic_columns)
    )]
    selected <- reference[selected, on = "MarkerName_hg38", nomatch = 0L]
    if (flip && nrow(selected)) {
      for (column in flip_columns) selected[, (column) := -get(column)]
    }
    selected
  }

  aligned <- rbindlist(list(
    map_one("id_same", FALSE),
    map_one("id_swap", TRUE),
    map_one("id_comp_same", FALSE),
    map_one("id_comp_swap", TRUE)
  ), use.names = TRUE, fill = TRUE)
  if (!nrow(aligned)) return(aligned)
  setorder(aligned, MarkerName_hg38, p_sort)
  aligned[, .SD[1], by = MarkerName_hg38]
}

read_reference <- function(gene) {
  bim_file <- file.path(reference_dir, paste0(gene, "_eur_hg38.bim"))
  raw_file <- file.path(reference_dir, paste0(gene, "_eur_hg38_012.raw"))
  if (!file.exists(bim_file) || !file.exists(raw_file)) return(NULL)

  bim <- fread(bim_file, header = FALSE)
  if (ncol(bim) != 6L) stop("Unexpected BIM format for ", gene)
  setnames(bim, c("chr", "MarkerName_hg38", "cm", "position", "a1", "a2"))
  genotype <- fread(raw_file)
  if (ncol(genotype) < 7L) stop("Unexpected RAW format for ", gene)
  genotype <- genotype[, 7:ncol(genotype)]
  if (ncol(genotype) != nrow(bim)) stop("BIM/RAW variant count mismatch for ", gene)
  setnames(genotype, bim$MarkerName_hg38)
  list(bim = bim, genotype = genotype)
}

empty_result <- function(analysis, outcome, model_type, gene, status, reason) {
  data.table(
    analysis = analysis, outcome = outcome, model_type = model_type,
    gene_id = gene, chr = NA_character_, training_model_type = NA_character_,
    input_snp_num = NA_integer_, matched_snp_num = NA_integer_,
    Z_cell_1 = NA_real_, p_cell_1 = NA_real_, Z_cell_2 = NA_real_,
    p_cell_2 = NA_real_, p_join = NA_real_, Z_sep_cell_1 = NA_real_,
    p_sep_cell_1 = NA_real_, Z_sep_cell_2 = NA_real_,
    p_sep_cell_2 = NA_real_, reg_scale_selected = NA_real_,
    reg_condition = NA_real_, association_mode = NA_character_,
    status = status, reason = reason
  )
}

run_one <- function(analysis, outcome, model_type, gene, weights_gene,
                    reference, aligned, beta_col, se_col,
                    family, n_control = NULL, n_case = NULL) {
  if (is.null(reference)) {
    return(empty_result(analysis, outcome, model_type, gene, "skipped", "missing_reference"))
  }
  if (!nrow(aligned)) {
    return(empty_result(analysis, outcome, model_type, gene, "skipped", "no_gwas_reference_overlap"))
  }

  merged <- merge(weights_gene, aligned, by = "MarkerName_hg38")
  merged <- merged[MarkerName_hg38 %chin% colnames(reference$genotype)]
  merged <- merged[is.finite(get(beta_col)) & is.finite(get(se_col)) & get(se_col) > 0]
  if (!nrow(merged)) {
    return(empty_result(analysis, outcome, model_type, gene, "skipped", "no_usable_weight_gwas_overlap"))
  }

  setorder(merged, MarkerName_hg38)
  genotype <- as.matrix(reference$genotype[, merged$MarkerName_hg38, with = FALSE])
  W <- as.matrix(merged[, .(weight_cell_1, weight_cell_2)])
  result <- tryCatch(
    SMiXcan_assoc_test_K(
      W = W,
      gwas_results = list(
        Beta = merged[[beta_col]],
        se_Beta = merged[[se_col]]
      ),
      x_g = genotype,
      n0 = n_control,
      n1 = n_case,
      family = family,
      regularization = regularization_mode,
      reg_scale = regularization_scale,
      condition_max = condition_max
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    return(empty_result(analysis, outcome, model_type, gene, "failed", result$message))
  }

  data.table(
    analysis = analysis,
    outcome = outcome,
    model_type = model_type,
    gene_id = gene,
    chr = unique(weights_gene$chr)[1],
    training_model_type = paste(sort(unique(weights_gene$type)), collapse = ";"),
    input_snp_num = nrow(weights_gene),
    matched_snp_num = nrow(merged),
    Z_cell_1 = result$Z_join[[1]],
    p_cell_1 = result$p_join_vec[[1]],
    Z_cell_2 = result$Z_join[[2]],
    p_cell_2 = result$p_join_vec[[2]],
    p_join = result$p_join,
    Z_sep_cell_1 = result$Z_sep[[1]],
    p_sep_cell_1 = result$p_sep[[1]],
    Z_sep_cell_2 = result$Z_sep[[2]],
    p_sep_cell_2 = result$p_sep[[2]],
    reg_scale_selected = result$reg_scale_selected,
    reg_condition = result$reg_condition,
    association_mode = result$mode,
    status = "completed",
    reason = NA_character_
  )
}

weights <- fread(normalized_weights_file)

md_data <- lapply(names(md_files), function(phenotype) {
  file <- file.path(data_dir, md_files[[phenotype]])
  if (!file.exists(file)) stop("Missing MD input: ", file)
  data <- fread(file)
  data[, chr_norm := normalize_chr(chr)]
  list(name = phenotype, data = data)
})
names(md_data) <- names(md_files)

bcac <- fread(bcac_overall_file)
if ("V1" %chin% names(bcac)) bcac[, V1 := NULL]
bcac[, chr_norm := normalize_chr(CHR)]

subtype <- fread(bcac_subtype_file)
subtype_columns <- c(
  "Luminal_A", "Luminal_B", "Luminal_B_HER2Neg",
  "HER2_Enriched", "Triple_Neg"
)
missing_subtypes <- setdiff(subtype_columns, names(subtype))
if (length(missing_subtypes)) stop("Missing subtype columns: ", paste(missing_subtypes, collapse = ", "))
if (!all(c("Baseline.Meta", "Effect.Meta") %chin% names(subtype))) {
  parsed <- tstrsplit(subtype$var_name, "_", fixed = TRUE)
  subtype[, Baseline.Meta := toupper(parsed[[3]])]
  subtype[, Effect.Meta := toupper(parsed[[4]])]
}
if (!"CHR" %chin% names(subtype) && "chr" %chin% names(subtype)) {
  setnames(subtype, "chr", "CHR")
}
subtype[, chr_norm := normalize_chr(CHR)]

result_list <- list()
index <- 1L
groups <- unique(weights[, .(model_type, gene_id)])

for (row_index in seq_len(nrow(groups))) {
  model <- groups$model_type[[row_index]]
  gene <- groups$gene_id[[row_index]]
  cat(sprintf("Model %s | gene %s\n", model, gene))
  weight_gene <- weights[model_type == model & gene_id == gene]
  reference <- read_reference(gene)

  if (is.null(reference)) {
    for (phenotype in names(md_data)) {
      result_list[[index]] <- empty_result("MD", phenotype, model, gene, "skipped", "missing_reference")
      index <- index + 1L
    }
    result_list[[index]] <- empty_result("BCAC_overall", "Breast_cancer", model, gene, "skipped", "missing_reference")
    index <- index + 1L
    for (subtype_name in subtype_columns) {
      result_list[[index]] <- empty_result("BCAC_subtype", subtype_name, model, gene, "skipped", "missing_reference")
      index <- index + 1L
    }
    next
  }

  target_chr <- normalize_chr(reference$bim$chr[[1]])
  reference_markers <- reference$bim$MarkerName_hg38

  for (phenotype in names(md_data)) {
    md_chr <- md_data[[phenotype]]$data[chr_norm == target_chr]
    # MD Zscore is the effect of Allele1. The reference marker and PLINK
    # dosage are ordered non-effect/effect, so pass Allele2 before Allele1.
    aligned <- harmonize_to_reference(
      md_chr, reference_markers,
      chr_col = "chr", pos_col = "POS_hg38",
      a1_col = "Allele2", a2_col = "Allele1",
      statistic_columns = "Zscore", p_col = "P-value"
    )
    if (nrow(aligned)) aligned[, se_input := 1]
    result_list[[index]] <- run_one(
      "MD", phenotype, model, gene, weight_gene, reference, aligned,
      beta_col = "Zscore", se_col = "se_input", family = "gaussian"
    )
    index <- index + 1L
  }

  bcac_chr <- bcac[chr_norm == target_chr]
  aligned_bcac <- harmonize_to_reference(
    bcac_chr, reference_markers,
    chr_col = "CHR", pos_col = "POS_hg38",
    a1_col = "Baseline.Meta", a2_col = "Effect.Meta",
    statistic_columns = c("Beta.meta", "sdE.meta"),
    flip_columns = "Beta.meta"
  )
  result_list[[index]] <- run_one(
    "BCAC_overall", "Breast_cancer", model, gene, weight_gene,
    reference, aligned_bcac, beta_col = "Beta.meta", se_col = "sdE.meta",
    family = "binomial", n_control = bcac_n_control, n_case = bcac_n_case
  )
  index <- index + 1L

  subtype_chr <- subtype[chr_norm == target_chr]
  aligned_subtype <- harmonize_to_reference(
    subtype_chr, reference_markers,
    chr_col = "CHR", pos_col = "POS_hg38",
    a1_col = "Baseline.Meta", a2_col = "Effect.Meta",
    statistic_columns = subtype_columns
  )
  if (nrow(aligned_subtype)) aligned_subtype[, se_input := 1]
  for (subtype_name in subtype_columns) {
    subtype_input <- copy(aligned_subtype)
    if (nrow(subtype_input)) subtype_input[, z_input := get(subtype_name)]
    result_list[[index]] <- run_one(
      "BCAC_subtype", subtype_name, model, gene, weight_gene,
      reference, subtype_input, beta_col = "z_input", se_col = "se_input",
      family = "binomial", n_control = bcac_n_control, n_case = bcac_n_case
    )
    index <- index + 1L
  }
}

results <- rbindlist(result_list, fill = TRUE)
skipped_file <- file.path(results_dir, "skipped_models.csv")
skipped_results <- results[status != "completed"]
if (nrow(skipped_results)) {
  fwrite(skipped_results, skipped_file)
} else if (file.exists(skipped_file)) {
  unlink(skipped_file)
}

output_columns <- c(
  "analysis", "outcome", "model_type", "gene_id", "chr",
  "training_model_type", "input_snp_num", "matched_snp_num",
  "Z_cell_1", "p_cell_1", "Z_cell_2", "p_cell_2", "p_join"
)
completed_results <- results[status == "completed", ..output_columns]
fwrite(completed_results[analysis == "MD"], file.path(results_dir, "md_results.csv"))
fwrite(
  completed_results[analysis == "BCAC_overall"],
  file.path(results_dir, "bcac_overall_results.csv")
)
fwrite(
  completed_results[analysis == "BCAC_subtype"],
  file.path(results_dir, "bcac_subtype_results.csv")
)

cat("Association workflow finished.\n")
print(results[, .N, by = .(analysis, status)])
