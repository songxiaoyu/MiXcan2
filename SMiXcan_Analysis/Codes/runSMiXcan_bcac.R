# ============================================================
# BCAC S-MiXcan (include indels + palindromic) + add predixcanlike
# ============================================================

library(data.table)
library(lme4)
library(glmnet)
library(doRNG)
library(ACAT)
library(MiXcan)
library(tibble)
library(tidyr)
library(dplyr)
library(MASS)
library(SMiXcan)

setwd('/Users/zhusinan/Downloads/adriana/')
dir_base   <- 'Result3/Ref/'
dir_output <- "Result3/Result_bcac/"
dir.create(dir_output, showWarnings = FALSE, recursive = TRUE)

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

# Harmonize BCAC GWAS to BIM MarkerName_hg38 (chr:pos:ref:alt) using 4-way mapping.
# Flips Beta when swapped or comp+swapped.
harmonize_bcac_to_bim_by_id <- function(gwas_chr, ld_bim,
                                        chr_col="CHR", pos_col="POS_hg38",
                                        base_col="Baseline.Meta", eff_col="Effect.Meta",
                                        beta_col="Beta.meta", se_col="sdE.meta",
                                        p_col=NULL) {

  g <- copy(gwas_chr)

  g[, chr_norm := normalize_chr(get(chr_col))]
  g[, pos      := as.character(get(pos_col))]
  g[, A1       := toupper(as.character(get(base_col)))]
  g[, A2       := toupper(as.character(get(eff_col)))]
  g[, B        := suppressWarnings(as.numeric(get(beta_col)))]
  g[, SE       := suppressWarnings(as.numeric(get(se_col)))]

  if (!is.null(p_col) && p_col %chin% names(g)) {
    g[, P := suppressWarnings(as.numeric(get(p_col)))]
  } else {
    g[, P := NA_real_]
  }

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

  b <- ld_bim[, .(MarkerName_hg38)]
  setkey(b, MarkerName_hg38)

  join_one <- function(dt, id_col, flip_beta) {
    tmp <- dt[!is.na(get(id_col)), .(MarkerName_hg38 = get(id_col), B, SE, P)]
    if (nrow(tmp) == 0) return(tmp[0])
    tmp <- b[tmp, on="MarkerName_hg38", nomatch=0L]
    if (nrow(tmp) == 0) return(tmp)
    tmp[, B := if (flip_beta) -B else B]
    tmp
  }

  out <- rbindlist(list(
    join_one(g, "id_same",      flip_beta = FALSE),
    join_one(g, "id_swap",      flip_beta = TRUE),
    join_one(g, "id_comp_same", flip_beta = FALSE),
    join_one(g, "id_comp_swap", flip_beta = TRUE)
  ), use.names = TRUE, fill = TRUE)

  if (nrow(out) == 0) {
    return(out[, .(MarkerName_hg38, Beta_aligned = B, se_aligned = SE)])
  }

  setorder(out, MarkerName_hg38, P)
  out <- out[, .SD[1], by = MarkerName_hg38]

  out[, .(MarkerName_hg38, Beta_aligned = B, se_aligned = SE)]
}

# --- counting helpers ---
parse_ref_alt <- function(marker) {
  p <- tstrsplit(marker, ":", fixed = TRUE)
  data.table(ref = toupper(p[[3]]), alt = toupper(p[[4]]))
}
is_indel <- function(ref, alt) nchar(ref) != 1 | nchar(alt) != 1
is_palindromic <- function(ref, alt) {
  (ref=="A" & alt=="T") | (ref=="T" & alt=="A") |
    (ref=="C" & alt=="G") | (ref=="G" & alt=="C")
}

# ----------------------------
# Inputs
# ----------------------------
# Add predixcanlike as a "model"
model_cfg <- list(
  list(model_type = "adipose",      weight_file = "subset_weight_mixcan2_adipose.csv"),
  list(model_type = "fibroblast",   weight_file = "subset_weight_mixcan2_fibroblast.csv"),
  list(model_type = "epithelial",   weight_file = "subset_weight_mixcan2_epithelial.csv"),
  list(model_type = "predixcanlike",weight_file = "subset_weight_predixcanlike.csv")
)

result_total <- vector(mode = "list", length = length(model_cfg))

# Read BCAC GWAS (DO NOT filter out indels/palindromic)
gwas_raw <- fread("/Users/zhusinan/Downloads/adriana/plink_snplist_by_gene/bcac_2020_sumstats_for_s-mixcan_hg38.csv")
if ("V1" %chin% names(gwas_raw)) gwas_raw[, V1 := NULL]
gwas_raw[, chr_norm := normalize_chr(CHR)]
gwas_raw[, Baseline.Meta := toupper(as.character(Baseline.Meta))]
gwas_raw[, Effect.Meta   := toupper(as.character(Effect.Meta))]

t <- 1L
for (cfg in model_cfg) {

  model_type  <- cfg$model_type
  weight_file <- cfg$weight_file
  cat(sprintf("\n===== Running model_type: %s | weights: %s =====\n", model_type, weight_file))

  W <- fread(weight_file)
  W[, MarkerName_hg38 := to_plink_id(xNameMatrix)]
  W[, chr_num := tstrsplit(MarkerName_hg38, ":", fixed = TRUE)[[1]]]

  gene_list <- unique(W$ID)
  G <- length(gene_list)

  result <- data.table(
    gene_id          = gene_list,
    chr              = NA_character_,
    type             = NA_character_,
    model_type       = model_type,
    MiXcan_snp_num    = NA_integer_,
    filtered_snp_num  = NA_integer_,
    used_indel        = NA_integer_,
    used_palindromic  = NA_integer_,
    Z_1_join          = NA_real_,
    p_1_join          = NA_real_,
    Z_2_join          = NA_real_,
    p_2_join          = NA_real_,
    p_join            = NA_real_
  )

  for (g in seq_len(G)) {

    gene <- gene_list[g]
    cat(sprintf("Processing gene: %s | model=%s\n", gene, model_type))

    W_s <- W[ID == gene]
    result[g, chr := as.character(unique(W_s$chr_num))[1]]
    result[g, type := as.character(unique(W_s$type))[1]]
    result[g, MiXcan_snp_num := nrow(W_s)]

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
        warning(sprintf("No BCAC GWAS records on chr %s for gene %s.", target_chr, gene))
        next
      }

      # Harmonize GWAS -> BIM
      gwas_aligned <- harmonize_bcac_to_bim_by_id(
        gwas_chr = gwas_chr,
        ld_bim   = ld_bim,
        chr_col  = "CHR",
        pos_col  = "POS_hg38",
        base_col = "Baseline.Meta",
        eff_col  = "Effect.Meta",
        beta_col = "Beta.meta",
        se_col   = "sdE.meta",
        p_col    = NULL
      )

      if (nrow(gwas_aligned) == 0) {
        warning(sprintf("No matched variants after harmonization for gene %s.", gene))
        next
      }

      # Merge BIM + aligned GWAS + weights
      total_input <- merge(
        merge(ld_bim, gwas_aligned, by = "MarkerName_hg38"),
        W_s,
        by = "MarkerName_hg38"
      )

      if (nrow(total_input) == 0) {
        warning(sprintf("No overlap among (BIM, GWAS, weights) for gene %s.", gene))
        next
      }

      # variants actually used
      ref_vars <- total_input$MarkerName_hg38
      ref_vars <- ref_vars[ref_vars %chin% colnames(X_ref)]
      if (length(ref_vars) == 0) {
        warning(sprintf("No genotype columns matched after filtering for gene %s.", gene))
        next
      }

      # count used indel & palindromic
      info <- parse_ref_alt(ref_vars)
      n_indel <- sum(is_indel(info$ref, info$alt))
      n_pal   <- sum(is_palindromic(info$ref, info$alt) & !is_indel(info$ref, info$alt))

      result[g, `:=`(
        filtered_snp_num = length(ref_vars),
        used_indel       = n_indel,
        used_palindromic = n_pal
      )]

      # enforce ordering
      total_input <- total_input[match(ref_vars, total_input$MarkerName_hg38)]
      X_ref_filtered <- as.matrix(X_ref[, ref_vars, with = FALSE])

      # weights
      if (!all(c("weight_cell_1","weight_cell_2") %chin% names(total_input))) {
        stop(sprintf("Missing weight_cell_1 / weight_cell_2 in weights file: %s", weight_file))
      }
      W1 <- as.matrix(total_input$weight_cell_1)
      W2 <- as.matrix(total_input$weight_cell_2)

      gwas_results <- list(
        Beta    = as.numeric(total_input$Beta_aligned),
        se_Beta = as.numeric(total_input$se_aligned)
      )

      SMiXcan_result <- SMiXcan_assoc_test(
        W1, W2, gwas_results, X_ref_filtered,
        n0 = 91477, n1 = 106278, family = 'binomial'
      )

      result[g, c('Z_1_join','p_1_join','Z_2_join','p_2_join','p_join')] <- SMiXcan_result

    }, error = function(e) {
      warning(sprintf("Gene %s failed: %s", gene, e$message))
    })
  }

  result_total[[t]] <- result
  fwrite(result, paste0(dir_output, model_type, "_bcac_SMiXcan_result.csv"))
  t <- t + 1L
}

# ----------------------------
# Merge + summary
# ----------------------------
merge_result <- rbindlist(result_total, fill = TRUE)
fwrite(merge_result, paste0(dir_output, "bcac_SMiXcan_result_merged_allmodels.csv"))

