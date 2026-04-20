# ============================================================
# FULL BCAC SUBTYPE S-MiXcan SCRIPT (HG38) — FIXED VERSION
# - Input GWAS file: tab-delimited, columns:
#   var_name, Effect.Meta, Baseline.Meta, Luminal_A, Luminal_B,
#   Luminal_B_HER2Neg, HER2_Enriched, Triple_Neg
# - Uses MarkerName_hg38 (chr:pos:ref:alt) as key
# - 4-way allele mapping: same / swapped / complement / comp+swapped
# - Keeps indels + palindromic (no exclusion)
# - Runs SMiXcan per subtype, per gene, per model_type:
#     adipose, fibroblast, epithelial, predixcanlike
# - Counts indel + palindromic among variants ACTUALLY USED (ref_vars)
# - Outputs:
#   (1) per-model results (long format: one row per gene × subtype)
#   (2) merged results across all models
#   (3) summary SNP composition table (all_models + each model_type)
# ============================================================

library(data.table)
library(SMiXcan)

# ----------------------------
# Paths
# ----------------------------
setwd('/Users/zhusinan/Downloads/adriana/')
dir_base   <- 'Result3/Ref'                 # per-gene .bim/.raw files
dir_output <- "Result3/Result_bcac_subtypes/"
dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

gwas_path <- "bcac_2020_subtype_zscores_for_s-mixcan_hg38.txt"  # <-- CHANGE to your actual filename/path

# ----------------------------
# Helpers
# ----------------------------
to_plink_id <- function(x) {
  x <- sub("^chr", "", x, ignore.case = TRUE)
  sub("^([^_]+)_([^_]+)_([^_]+)_([^_]+).*", "\\1:\\2:\\3:\\4", x)
}
normalize_chr <- function(chr) {
  chr <- as.character(chr)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr
}
is_snp_acgt <- function(a) {
  a <- toupper(as.character(a))
  nchar(a) == 1 & a %chin% c("A","C","G","T")
}
comp_allele <- function(a) {
  a <- toupper(as.character(a))
  map <- c(A="T", T="A", C="G", G="C")
  as.character(map[a])
}

# --- counting helpers (among USED predictors) ---
parse_ref_alt <- function(marker) {
  p <- tstrsplit(marker, ":", fixed = TRUE)
  data.table(ref = toupper(p[[3]]), alt = toupper(p[[4]]))
}
is_indel <- function(ref, alt) nchar(ref) != 1 | nchar(alt) != 1
is_palindromic <- function(ref, alt) {
  (ref=="A" & alt=="T") | (ref=="T" & alt=="A") |
    (ref=="C" & alt=="G") | (ref=="G" & alt=="C")
}

# ------------------------------------------------------------
# Harmonize subtype GWAS -> BIM MarkerName_hg38 (chr:pos:ref:alt)
# using 4-way matching and flip Z where needed.
#
# Returns:
#   data.table with columns: MarkerName_hg38 + subtype_cols (aligned Z per subtype)
# ------------------------------------------------------------
harmonize_subtypeZ_to_bim_by_id <- function(gwas_chr, ld_bim,
                                            subtype_cols,
                                            chr_col="CHR", pos_col="POS_hg38",
                                            base_col="Baseline.Meta", eff_col="Effect.Meta") {
  g <- copy(gwas_chr)

  g[, chr_norm := normalize_chr(get(chr_col))]
  g[, pos      := as.character(get(pos_col))]
  g[, A1       := toupper(as.character(get(base_col)))]
  g[, A2       := toupper(as.character(get(eff_col)))]

  # ensure subtype columns numeric
  for (sc in subtype_cols) {
    if (sc %chin% names(g)) g[, (sc) := suppressWarnings(as.numeric(get(sc)))]
  }

  # candidate IDs
  g[, id_same := paste0(chr_norm, ":", pos, ":", A1, ":", A2)]
  g[, id_swap := paste0(chr_norm, ":", pos, ":", A2, ":", A1)]

  okc <- is_snp_acgt(g$A1) & is_snp_acgt(g$A2)
  g[, id_comp_same := NA_character_]
  g[, id_comp_swap := NA_character_]
  if (any(okc)) {
    c1 <- comp_allele(g$A1[okc])
    c2 <- comp_allele(g$A2[okc])
    g$id_comp_same[okc] <- paste0(g$chr_norm[okc], ":", g$pos[okc], ":", c1, ":", c2)
    g$id_comp_swap[okc] <- paste0(g$chr_norm[okc], ":", g$pos[okc], ":", c2, ":", c1)
  }

  # BIM key
  b <- ld_bim[, .(MarkerName_hg38)]
  setkey(b, MarkerName_hg38)

  # IMPORTANT FIX: safe data.table selection with new column + subtype columns
  join_one <- function(dt, id_col, flip) {
    tmp <- dt[!is.na(get(id_col)), c(list(MarkerName_hg38 = get(id_col)), mget(subtype_cols))]
    if (nrow(tmp) == 0) return(tmp[0])
    tmp <- b[tmp, on="MarkerName_hg38", nomatch=0L]
    if (nrow(tmp) == 0) return(tmp)
    if (flip) {
      for (sc in subtype_cols) tmp[, (sc) := -get(sc)]
    }
    tmp
  }

  out <- rbindlist(list(
    join_one(g, "id_same",      flip = FALSE),
    join_one(g, "id_swap",      flip = TRUE),
    join_one(g, "id_comp_same", flip = FALSE),
    join_one(g, "id_comp_swap", flip = TRUE)
  ), use.names = TRUE, fill = TRUE)

  if (nrow(out) == 0) return(out)

  # Deduplicate per MarkerName_hg38:
  # Keep first non-NA across concatenated rows (same > swap > comp_same > comp_swap)
  setkey(out, MarkerName_hg38)
  out <- out[, lapply(.SD, function(x) {
    i <- which(!is.na(x))[1]
    if (length(i) == 0) NA_real_ else x[i]
  }), by = MarkerName_hg38, .SDcols = subtype_cols]

  out
}

# ----------------------------
# Load subtype GWAS
# ----------------------------
gwas_raw <- fread(gwas_path, sep = "\t")

subtype_list <- c("Luminal_A", "Luminal_B", "Luminal_B_HER2Neg", "HER2_Enriched", "Triple_Neg")
missing_sub <- setdiff(subtype_list, names(gwas_raw))
if (length(missing_sub) > 0) stop("Missing subtype columns in GWAS file: ", paste(missing_sub, collapse=", "))

# Parse var_name: chromosome_position_allele1_allele2
# gwas_raw already loaded with fread(...)
# var_name like: "9_46841825_A_G" or "6_32682510_I_D"

tmp <- tstrsplit(gwas_raw$var_name, "_", fixed = TRUE)
stopifnot(length(tmp) >= 4)

# gwas_raw[, CHR := paste0("chr", tmp[[1]])]

# Use allele1/allele2 from var_name
gwas_raw[, Baseline.Meta := toupper(tmp[[3]])]
gwas_raw[, Effect.Meta   := toupper(tmp[[4]])]

# Keep chr_norm for fast per-chr subset
gwas_raw[, chr_norm := normalize_chr(CHR)]

for (sc in subtype_list) gwas_raw[, (sc) := suppressWarnings(as.numeric(get(sc)))]

# ----------------------------
# Model configs: 3 celltypes + predixcanlike
# ----------------------------
model_cfg <- list(
  list(model_type="adipose",       weight_file="subset_weight_mixcan2_adipose.csv"),
  list(model_type="fibroblast",    weight_file="subset_weight_mixcan2_fibroblast.csv"),
  list(model_type="epithelial",    weight_file="subset_weight_mixcan2_epithelial.csv"),
  list(model_type="predixcanlike", weight_file="subset_weight_predixcanlike.csv")
)

all_results <- list()
k <- 1L

# ----------------------------
# Main loop
# ----------------------------
for (cfg in model_cfg) {

  model_type  <- cfg$model_type
  weight_file <- cfg$weight_file
  cat(sprintf("\n===== Model: %s | weights: %s =====\n", model_type, weight_file))

  W <- fread(weight_file)
  W[, MarkerName_hg38 := to_plink_id(xNameMatrix)]
  W[, chr_num := tstrsplit(MarkerName_hg38, ":", fixed = TRUE)[[1]]]

  gene_list <- unique(W$ID)

  # long-format results for this model
  res_model <- list()
  m <- 1L

  for (gene in gene_list) {

    W_s <- W[ID == gene]
    gene_chr  <- as.character(unique(W_s$chr_num))[1]
    gene_type <- as.character(unique(W_s$type))[1]
    mixcan_n  <- nrow(W_s)

    cat(sprintf("Processing gene: %s | model=%s\n", gene, model_type))

    bim_path <- file.path(dir_base, sprintf("%s_eur_hg38.bim", gene))
    raw_path <- file.path(dir_base, sprintf("%s_eur_hg38_012.raw", gene))
    if (!file.exists(bim_path) || !file.exists(raw_path)) {
      warning(sprintf("Missing BIM/RAW for gene %s. Skip.", gene))
      next
    }

    tryCatch({

      # BIM
      ld_bim <- fread(bim_path, header = FALSE)
      setnames(ld_bim, c('chr', 'MarkerName_hg38', 'mor', 'POS_hg38', 'A1_bim', 'A2_bim'))
      ld_bim[, chr_norm := normalize_chr(chr)]
      target_chr <- ld_bim$chr_norm[1]

      # RAW genotypes
      X_ref <- fread(raw_path, sep = " ")
      X_ref <- X_ref[, 7:ncol(X_ref)]
      setnames(X_ref, ld_bim$MarkerName_hg38)

      # GWAS restricted to chr
      gwas_chr <- gwas_raw[chr_norm == target_chr]
      if (nrow(gwas_chr) == 0) {
        warning(sprintf("No subtype GWAS records on chr %s for gene %s.", target_chr, gene))
        next
      }

      # Harmonize subtype Z columns to BIM MarkerName_hg38
      gwas_aligned <- harmonize_subtypeZ_to_bim_by_id(
        gwas_chr = gwas_chr,
        ld_bim   = ld_bim,
        subtype_cols = subtype_list,
        chr_col  = "chr",
        pos_col  = "POS_hg38",
        base_col = "Baseline.Meta",
        eff_col  = "Effect.Meta"
      )

      if (nrow(gwas_aligned) == 0) {
        warning(sprintf("No matched variants after harmonization for gene %s.", gene))
        next
      }

      # Merge BIM + aligned GWAS + weights by MarkerName_hg38
      total_input <- merge(
        merge(ld_bim, gwas_aligned, by = "MarkerName_hg38"),
        W_s,
        by = "MarkerName_hg38"
      )

      if (nrow(total_input) == 0) {
        warning(sprintf("No overlap among (BIM, GWAS, weights) for gene %s.", gene))
        next
      }

      # Variants actually USED (must exist in genotype columns)
      ref_vars <- total_input$MarkerName_hg38
      ref_vars <- ref_vars[ref_vars %chin% colnames(X_ref)]
      if (length(ref_vars) == 0) {
        warning(sprintf("No genotype columns matched after filtering for gene %s.", gene))
        next
      }

      # count used indel & palindromic
      info <- parse_ref_alt(ref_vars)
      used_indel <- sum(is_indel(info$ref, info$alt))
      used_pal   <- sum(is_palindromic(info$ref, info$alt) & !is_indel(info$ref, info$alt))

      # enforce consistent ordering
      total_input <- total_input[match(ref_vars, total_input$MarkerName_hg38)]
      X_ref_filtered <- as.matrix(X_ref[, ref_vars, with = FALSE])

      # weights
      if (!all(c("weight_cell_1","weight_cell_2") %chin% names(total_input))) {
        stop(sprintf("Missing weight_cell_1 / weight_cell_2 in weights file: %s", weight_file))
      }
      W1 <- as.matrix(total_input$weight_cell_1)
      W2 <- as.matrix(total_input$weight_cell_2)

      # Run S-MiXcan for each subtype (Z-score GWAS)
      for (subtype in subtype_list) {

        # IMPORTANT FIX: guard against missing/empty subtype column
        if (!subtype %chin% names(total_input)) next
        Z_sub <- suppressWarnings(as.numeric(total_input[[subtype]]))
        if (length(Z_sub) == 0 || all(is.na(Z_sub))) next

        gwas_results <- list(
          Beta    = Z_sub,
          se_Beta = rep(1, length(Z_sub))   # Z-scores -> se=1
        )

        SMiXcan_result <- SMiXcan_assoc_test(
          W1, W2, gwas_results, X_ref_filtered,
          n0= 91477, n1=106278, family='binomial'
        )

        res_model[[m]] <- data.table(
          model_type       = model_type,
          subtype          = subtype,
          gene_id          = gene,
          chr              = gene_chr,
          type             = gene_type,
          MiXcan_snp_num    = mixcan_n,
          filtered_snp_num  = length(ref_vars),
          used_indel        = used_indel,
          used_palindromic  = used_pal,
          Z_1_join          = SMiXcan_result[[1]],
          p_1_join          = SMiXcan_result[[2]],
          Z_2_join          = SMiXcan_result[[3]],
          p_2_join          = SMiXcan_result[[4]],
          p_join            = SMiXcan_result[[5]]
        )
        m <- m + 1L
      }

    }, error = function(e) {
      warning(sprintf("Gene %s failed: %s", gene, e$message))
    })
  }

  res_model_dt <- rbindlist(res_model, fill = TRUE)
  fwrite(res_model_dt, file.path(dir_output, paste0("bcac_subtypes_", model_type, "_SMiXcan_results_long.csv")))

  all_results[[k]] <- res_model_dt
  k <- k + 1L
}

# ----------------------------
# Merge all models
# ----------------------------
merge_result <- rbindlist(all_results, fill = TRUE)
fwrite(merge_result, file.path(dir_output, "bcac_subtypes_SMiXcan_results_long_ALLMODELS.csv"))

# ----------------------------
# SNP composition summary table (across all subtypes)
# To avoid double counting across subtypes, collapse by (model_type, gene_id).
# ----------------------------
uniq_gene <- unique(merge_result[, .(model_type, gene_id, filtered_snp_num, used_indel, used_palindromic)])

summary_by_model <- uniq_gene[, .(
  total       = sum(filtered_snp_num, na.rm = TRUE),
  indels      = sum(used_indel, na.rm = TRUE),
  palindromic = sum(used_palindromic, na.rm = TRUE)
), by = .(model_type)]

summary_by_model[, indels := sprintf("%d (%.1f%%)", indels, 100 * indels / total)]
summary_by_model[, palindromic := sprintf("%d (%.1f%%)", palindromic, 100 * palindromic / total)]

summary_all <- uniq_gene[, .(
  model_type  = "all_models",
  total       = sum(filtered_snp_num, na.rm = TRUE),
  indels      = sum(used_indel, na.rm = TRUE),
  palindromic = sum(used_palindromic, na.rm = TRUE)
)]
summary_all[, indels := sprintf("%d", indels)]
summary_all[, palindromic := sprintf("%d", palindromic)]

summary_table <- rbindlist(list(summary_all, summary_by_model), use.names = TRUE, fill = TRUE)
fwrite(summary_table, file.path(dir_output, "bcac_subtypes_SNP_composition_summary.csv"))
print(summary_table)

# ----------------------------
# Optional sanity check: genes tested per subtype
# ----------------------------
tested_gene_counts <- merge_result[, .(n_genes = uniqueN(gene_id)), by = .(model_type, subtype)]
fwrite(tested_gene_counts, file.path(dir_output, "bcac_subtypes_tested_gene_counts.csv"))
print(tested_gene_counts)
