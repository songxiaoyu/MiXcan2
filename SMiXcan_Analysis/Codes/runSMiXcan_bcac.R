# ============================================================
# BCAC S-MiXcan (include indels + palindromic) + add predixcanlike
# ============================================================

library(data.table)
library(SMiXcan)

analysis_dir <- normalizePath(Sys.getenv("MIXCAN_ANALYSIS_DIR", unset = "/Users/zhusinan/Downloads/adriana"), mustWork = TRUE)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  path <- sub(file_arg, "", args[startsWith(args, file_arg)])
  if (length(path) > 0) return(dirname(normalizePath(path[1], mustWork = FALSE)))
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)))
  }
  getwd()
}

repo_data_dir <- normalizePath(file.path(get_script_dir(), "..", "..", "Data"), mustWork = TRUE)
weight_dir <- repo_data_dir
dir_base <- file.path(repo_data_dir, "1000Genome_Ref")

setwd(analysis_dir)
dir_output <- "Result/Result_bcac/"
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

build_bcac_candidate_table <- function(gwas_chr,
                                       chr_col = "CHR", pos_col = "POS_hg38",
                                       base_col = "Baseline.Meta", eff_col = "Effect.Meta",
                                       beta_col = "Beta.meta", se_col = "sdE.meta",
                                       p_col = NULL,
                                       marker_filter = NULL) {
  g <- copy(gwas_chr)

  g[, chr_norm := normalize_chr(get(chr_col))]
  g[, pos      := as.character(get(pos_col))]
  g[, A1       := toupper(as.character(get(base_col)))]
  g[, A2       := toupper(as.character(get(eff_col)))]
  g[, B        := suppressWarnings(as.numeric(get(beta_col)))]
  g[, SE       := suppressWarnings(as.numeric(get(se_col)))]
  g[, source_order := .I]

  if (!is.null(p_col) && p_col %chin% names(g)) {
    g[, P := suppressWarnings(as.numeric(get(p_col)))]
  } else {
    g[, P := NA_real_]
  }

  make_candidate <- function(dt, marker, flip_beta, priority) {
    out <- dt[!is.na(marker), .(
      MarkerName_hg38 = marker[!is.na(marker)],
      B = B[!is.na(marker)],
      SE = SE[!is.na(marker)],
      P = P[!is.na(marker)],
      priority = priority,
      source_order = source_order[!is.na(marker)]
    )]
    if (flip_beta && nrow(out) > 0) out[, B := -B]
    out
  }

  id_same <- paste0(g$chr_norm, ":", g$pos, ":", g$A1, ":", g$A2)
  id_swap <- paste0(g$chr_norm, ":", g$pos, ":", g$A2, ":", g$A1)

  id_comp_same <- rep(NA_character_, nrow(g))
  id_comp_swap <- rep(NA_character_, nrow(g))
  okc <- is_snp_acgt(g$A1) & is_snp_acgt(g$A2)
  if (any(okc)) {
    c1 <- comp_allele(g$A1[okc])
    c2 <- comp_allele(g$A2[okc])
    id_comp_same[okc] <- paste0(g$chr_norm[okc], ":", g$pos[okc], ":", c1, ":", c2)
    id_comp_swap[okc] <- paste0(g$chr_norm[okc], ":", g$pos[okc], ":", c2, ":", c1)
  }

  if (!is.null(marker_filter)) {
    keep <- id_same %chin% marker_filter |
      id_swap %chin% marker_filter |
      id_comp_same %chin% marker_filter |
      id_comp_swap %chin% marker_filter
    if (!any(keep)) return(data.table())
    g <- g[keep]
    id_same <- id_same[keep]
    id_swap <- id_swap[keep]
    id_comp_same <- id_comp_same[keep]
    id_comp_swap <- id_comp_swap[keep]
  }

  out <- rbindlist(list(
    make_candidate(g, id_same,      flip_beta = FALSE, priority = 1L),
    make_candidate(g, id_swap,      flip_beta = TRUE,  priority = 2L),
    make_candidate(g, id_comp_same, flip_beta = FALSE, priority = 3L),
    make_candidate(g, id_comp_swap, flip_beta = TRUE,  priority = 4L)
  ), use.names = TRUE, fill = TRUE)

  if (nrow(out) == 0) return(out)

  setorder(out, MarkerName_hg38, P, priority, source_order)
  out <- out[, .SD[1], by = MarkerName_hg38]
  setkey(out, MarkerName_hg38)
  out[, .(MarkerName_hg38, Beta_aligned = B, se_aligned = SE)]
}

harmonize_bcac_to_bim_cached <- function(gwas_candidates, ld_bim) {
  b <- ld_bim[, .(MarkerName_hg38)]
  setkey(b, MarkerName_hg38)
  out <- gwas_candidates[b, on = "MarkerName_hg38", nomatch = 0L]
  if (nrow(out) == 0) return(out[, .(MarkerName_hg38, Beta_aligned, se_aligned)])
  out[, .(MarkerName_hg38, Beta_aligned, se_aligned)]
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

run_smixcan_assoc <- function(W1, W2, gwas_results, X_ref, n0 = NULL, n1 = NULL,
                              family = c("binomial", "gaussian")) {
  family <- match.arg(family)
  W <- cbind(as.numeric(W1), as.numeric(W2))
  res <- SMiXcan_assoc_test_K(
    W = W,
    gwas_results = gwas_results,
    x_g = X_ref,
    n0 = n0,
    n1 = n1,
    family = family,
    regularization = "fixed",
    reg_scale = 0.001,
    weight_eps = 0
  )
  c(res$Z_join[1], res$p_join_vec[1], res$Z_join[2], res$p_join_vec[2], res$p_join)
}

# ----------------------------
# Inputs
# ----------------------------
# Add predixcanlike as a "model"
model_cfg <- list(
  list(model_type = "adipose",      weight_file = file.path(weight_dir, "subset_weight_mixcan2_adipose.csv")),
  list(model_type = "fibroblast",   weight_file = file.path(weight_dir, "subset_weight_mixcan2_fibroblast.csv")),
  list(model_type = "epithelial",   weight_file = file.path(weight_dir, "subset_weight_mixcan2_epithelial.csv")),
  list(model_type = "predixcanlike",weight_file = file.path(weight_dir, "subset_weight_predixcanlike.csv"))
)

result_total <- vector(mode = "list", length = length(model_cfg))

build_target_markers_by_chr <- function(ref_dir) {
  bim_files <- list.files(ref_dir, pattern = "[.]bim$", full.names = TRUE)
  targets <- rbindlist(lapply(bim_files, function(path) {
    x <- fread(path, header = FALSE, select = c(1, 2))
    setnames(x, c("chr", "MarkerName_hg38"))
    x[, chr_norm := normalize_chr(chr)]
    x[, .(chr_norm, MarkerName_hg38)]
  }), use.names = TRUE)
  targets <- unique(targets)
  by_chr <- targets[, .(markers = list(MarkerName_hg38)), by = chr_norm]
  setNames(by_chr$markers, by_chr$chr_norm)
}

target_markers_by_chr <- build_target_markers_by_chr(dir_base)

# Read BCAC GWAS (DO NOT filter out indels/palindromic)
gwas_path <- file.path(repo_data_dir, "bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv")
gwas_raw <- fread(gwas_path)
if ("V1" %chin% names(gwas_raw)) gwas_raw[, V1 := NULL]
gwas_raw[, chr_norm := normalize_chr(CHR)]
gwas_raw[, Baseline.Meta := toupper(as.character(Baseline.Meta))]
gwas_raw[, Effect.Meta   := toupper(as.character(Effect.Meta))]

gwas_candidate_cache <- new.env(parent = emptyenv())
get_gwas_candidates_for_chr <- function(chr_norm) {
  cache_key <- as.character(chr_norm)
  if (!exists(cache_key, envir = gwas_candidate_cache, inherits = FALSE)) {
    gwas_chr <- gwas_raw[chr_norm == cache_key]
    marker_filter <- target_markers_by_chr[[cache_key]]
    if (nrow(gwas_chr) == 0 || length(marker_filter) == 0) {
      assign(cache_key, data.table(), envir = gwas_candidate_cache)
    } else {
      assign(
        cache_key,
        build_bcac_candidate_table(
          gwas_chr = gwas_chr,
          chr_col  = "CHR",
          pos_col  = "POS_hg38",
          base_col = "Baseline.Meta",
          eff_col  = "Effect.Meta",
          beta_col = "Beta.meta",
          se_col   = "sdE.meta",
          p_col    = NULL,
          marker_filter = marker_filter
        ),
        envir = gwas_candidate_cache
      )
    }
  }
  get(cache_key, envir = gwas_candidate_cache, inherits = FALSE)
}

ref_data_cache <- new.env(parent = emptyenv())
get_ref_data_for_gene <- function(gene) {
  if (!exists(gene, envir = ref_data_cache, inherits = FALSE)) {
    bim_path <- file.path(dir_base, sprintf("%s_eur_hg38.bim", gene))
    raw_path <- file.path(dir_base, sprintf("%s_eur_hg38_012.raw", gene))
    if (!file.exists(bim_path) || !file.exists(raw_path)) {
      assign(gene, NULL, envir = ref_data_cache)
    } else {
      ld_bim <- fread(bim_path, header = FALSE)
      setnames(ld_bim, c('chr', 'MarkerName_hg38', 'mor', 'POS_hg38', 'A1_bim', 'A2_bim'))
      ld_bim[, chr_norm := normalize_chr(chr)]

      X_ref <- fread(raw_path, sep = " ")
      X_ref <- X_ref[, 7:ncol(X_ref)]
      setnames(X_ref, ld_bim$MarkerName_hg38)

      assign(gene, list(ld_bim = ld_bim, X_ref = X_ref), envir = ref_data_cache)
    }
  }
  get(gene, envir = ref_data_cache, inherits = FALSE)
}

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

    ref_data <- get_ref_data_for_gene(gene)
    if (is.null(ref_data)) {
      warning(sprintf("Missing BIM/RAW for gene %s. Skip.", gene))
      next
    }

    tryCatch({

      ld_bim <- ref_data$ld_bim
      X_ref <- ref_data$X_ref
      target_chr <- ld_bim$chr_norm[1]

      # GWAS candidates are built once per chromosome and reused across genes/models.
      gwas_candidates <- get_gwas_candidates_for_chr(target_chr)
      if (nrow(gwas_candidates) == 0) {
        warning(sprintf("No BCAC GWAS records on chr %s for gene %s.", target_chr, gene))
        next
      }

      # Harmonize GWAS -> BIM
      gwas_aligned <- harmonize_bcac_to_bim_cached(
        gwas_candidates = gwas_candidates,
        ld_bim = ld_bim
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

      SMiXcan_result <- run_smixcan_assoc(
        W1, W2, gwas_results, X_ref_filtered,
        n0 = 91477, n1 = 106278, family = 'binomial'
      )

      result[g, `:=`(
        Z_1_join = SMiXcan_result[[1]],
        p_1_join = SMiXcan_result[[2]],
        Z_2_join = SMiXcan_result[[3]],
        p_2_join = SMiXcan_result[[4]],
        p_join   = SMiXcan_result[[5]]
      )]

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
