# ============================================================
# S-MiXcan pipeline (include indels + palindromic) + counting
# Counts indel & palindromic predictors actually USED (ref_vars)
# ============================================================

library(data.table)
library(data.table)
library(lme4)
library(glmnet)
library(doRNG)
library(ACAT)
library(tibble)
library(tidyr)
library(dplyr)
library(MASS)
library(SMiXcan)

setwd('/Users/zhusinan/Downloads/adriana/')

# Input reference (per-gene plink files)
dir_base   <- 'Result3/Ref/'
# Output path
dir_output <- "Result3/Result_PrediXcan/"
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

# Harmonize GWAS to BIM MarkerName_hg38 (chr:pos:ref:alt) using 4-way matching.
# Builds candidate IDs from GWAS (same/swapped/complement/comp+swapped), matches to BIM IDs,
# flips Zscore when swapped or comp+swapped, keeps one record per BIM ID (min P).
harmonize_gwas_to_bim_by_id <- function(gwas_chr, ld_bim,
                                        chr_col="chr", pos_col="POS_hg38",
                                        a1_col="Allele1", a2_col="Allele2",
                                        z_col="Zscore", p_col="P-value") {
  g <- copy(gwas_chr)

  g[, chr_norm := normalize_chr(get(chr_col))]
  g[, pos      := as.character(get(pos_col))]
  g[, A1       := toupper(as.character(get(a1_col)))]
  g[, A2       := toupper(as.character(get(a2_col)))]
  g[, Z        := suppressWarnings(as.numeric(get(z_col)))]
  g[, P        := suppressWarnings(as.numeric(get(p_col)))]

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

  join_one <- function(dt, id_col, flip) {
    tmp <- dt[!is.na(get(id_col)), .(MarkerName_hg38 = get(id_col), Z, P)]
    if (nrow(tmp) == 0) return(tmp[0])
    tmp <- b[tmp, on="MarkerName_hg38", nomatch=0L]
    if (nrow(tmp) == 0) return(tmp)
    tmp[, Z := if (flip) -Z else Z]
    tmp
  }

  out <- rbindlist(list(
    join_one(g, "id_same",      flip = FALSE),
    join_one(g, "id_swap",      flip = TRUE),
    join_one(g, "id_comp_same", flip = FALSE),
    join_one(g, "id_comp_swap", flip = TRUE)
  ), use.names=TRUE, fill=TRUE)

  if (nrow(out) == 0) return(out[, .(MarkerName_hg38, Zscore_aligned=Z, Pvalue_aligned=P)])

  setorder(out, MarkerName_hg38, P)
  out <- out[, .SD[1], by=MarkerName_hg38]

  out[, .(MarkerName_hg38, Zscore_aligned=Z, Pvalue_aligned=P)]
}

# --- counting helpers (for variants actually used) ---
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
# Main
# ----------------------------
pheno_list <- c('DA', 'NDA', 'PMD')

result_total <- vector(mode = "list", length = length(pheno_list))
t <- 1L

# Read PrediXcan-like weights once
W <- fread("subset_weight_predixcanlike.csv")
W[, MarkerName_hg38 := to_plink_id(xNameMatrix)]
W[, chr_num := tstrsplit(MarkerName_hg38, ":", fixed = TRUE)[[1]]]
gene_list <- unique(W$ID)
G <- length(gene_list)

for (pheno in pheno_list) {

  # GWAS: DO NOT exclude indels / palindromic
  gwas <- fread(paste0("MD_results_2021/", pheno, "_MetaResultswInfo_20210609_hg38.txt"))
  gwas[, chr_norm := normalize_chr(chr)]
  gwas[, Allele1 := toupper(as.character(Allele1))]
  gwas[, Allele2 := toupper(as.character(Allele2))]

  result <- data.table(
    gene_id          = gene_list,
    chr              = NA_character_,
    type             = NA_character_,
    phenotype        = pheno,
    PrediXcan_snp_num = NA_integer_,
    filtered_snp_num  = NA_integer_,   # actually used
    used_indel        = NA_integer_,   # among used: indels
    used_palindromic  = NA_integer_,   # among used: palindromic SNPs
    Z_1_join          = NA_real_,
    p_1_join          = NA_real_,
    Z_2_join          = NA_real_,
    p_2_join          = NA_real_,
    p_join            = NA_real_
  )

  for (g in seq_len(G)) {

    gene <- gene_list[g]
    cat(sprintf("Processing gene: %s | pheno=%s\n", gene, pheno))

    W_s <- W[ID == gene]
    result[g, chr := as.character(unique(W_s$chr_num))[1]]
    result[g, type := as.character(unique(W_s$type))[1]]
    result[g, PrediXcan_snp_num := nrow(W_s)]

    # per-gene ref files
    bim_path <- file.path(dir_base, sprintf("%s_eur_hg38.bim", gene))
    raw_path <- file.path(dir_base, sprintf("%s_eur_hg38_012.raw", gene))
    if (!file.exists(bim_path) || !file.exists(raw_path)) {
      warning(sprintf("Missing BIM/RAW for gene %s. Skip.", gene))
      next
    }

    tryCatch({

      # BIM (MarkerName_hg38 is chr:pos:ref:alt)
      ld_bim <- fread(bim_path, header = FALSE)
      setnames(ld_bim, c('CHR', 'MarkerName_hg38', 'mor', 'POS_hg38', 'A1_bim', 'A2_bim'))
      ld_bim[, CHR_norm := normalize_chr(CHR)]
      target_chr <- ld_bim$CHR_norm[1]

      # RAW genotypes
      X_ref <- fread(raw_path, sep = " ")
      X_ref <- X_ref[, 7:ncol(X_ref)]
      setnames(X_ref, ld_bim$MarkerName_hg38)

      # GWAS restricted to gene chr
      gwas_chr <- gwas[chr_norm == target_chr]
      if (nrow(gwas_chr) == 0) {
        warning(sprintf("No GWAS records on chr %s for gene %s.", target_chr, gene))
        next
      }

      # Harmonize GWAS to BIM by MarkerName_hg38
      gwas_aligned <- harmonize_gwas_to_bim_by_id(
        gwas_chr = gwas_chr,
        ld_bim   = ld_bim,
        chr_col  = "chr",
        pos_col  = "POS_hg38",
        a1_col   = "Allele1",
        a2_col   = "Allele2",
        z_col    = "Zscore",
        p_col    = "P-value"
      )

      if (nrow(gwas_aligned) == 0) {
        warning(sprintf("No matched variants after harmonization for gene %s.", gene))
        next
      }

      # Merge (BIM + aligned GWAS) and then with weights by MarkerName_hg38
      total_input <- merge(
        merge(ld_bim, gwas_aligned, by = "MarkerName_hg38"),
        W_s,
        by = "MarkerName_hg38"
      )

      if (nrow(total_input) == 0) {
        warning(sprintf("No overlap among (BIM, GWAS, weights) for gene %s.", gene))
        next
      }

      # final set of variants actually USED
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

      # consistent ordering
      total_input <- total_input[match(ref_vars, total_input$MarkerName_hg38)]
      X_ref_filtered <- as.matrix(X_ref[, ref_vars, with = FALSE])

      # weights: for predixcanlike it might still have weight_cell_1/2 columns in your file
      # If your predixcanlike file uses a single weight column, tell me its column name.
      if (!all(c("weight_cell_1","weight_cell_2") %chin% names(total_input))) {
        stop("Missing weight_cell_1 / weight_cell_2 in predixcanlike weights file.")
      }
      W1 <- as.matrix(total_input$weight_cell_1)
      W2 <- as.matrix(total_input$weight_cell_2)

      gwas_results <- list(
        Beta    = as.numeric(total_input$Zscore_aligned),
        se_Beta = rep(1, length(ref_vars))
      )

      SMiXcan_result <- SMiXcan_assoc_test(
        W1, W2, gwas_results, X_ref_filtered,
        n0=NULL, n1=NULL, family='gaussian'
      )

      result[g, c('Z_1_join','p_1_join','Z_2_join','p_2_join','p_join')] <- SMiXcan_result

    }, error = function(e) {
      warning(sprintf("Gene %s failed: %s", gene, e$message))
    })
  }

  result_total[[t]] <- result
  t <- t + 1L
  fwrite(result, paste0(dir_output, pheno, "_PrediXcan_result.csv"))
}

merge_result <- rbindlist(result_total, fill = TRUE)
fwrite(merge_result, paste0(dir_output, "PrediXcan_result_merged.csv"))

# ----------------------------
# One-line summary table for predixcanlike
# ----------------------------
summary_predixcan <- merge_result[, .(
  total       = sum(filtered_snp_num, na.rm = TRUE)/3,
  indels      = sum(used_indel, na.rm = TRUE)/3,
  palindromic = sum(used_palindromic, na.rm = TRUE)/3
)]

summary_predixcan[, indels := sprintf("%d (%.1f%%)", indels, 100 * indels / total)]
summary_predixcan[, palindromic := sprintf("%d (%.1f%%)", palindromic, 100 * palindromic / total)]
summary_predixcan[, model_type := "predixcanlike"]
summary_predixcan <- summary_predixcan[, .(model_type, total, indels, palindromic)]

fwrite(summary_predixcan, paste0(dir_output, "SNP_composition_predixcanlike_summary.csv"))
print(summary_predixcan)

