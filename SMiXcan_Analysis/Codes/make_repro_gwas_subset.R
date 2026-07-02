library(data.table)

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
ref_dir <- file.path(repo_data_dir, "1000Genome_Ref")

normalize_chr <- function(chr) {
  chr <- as.character(chr)
  sub("^chr", "", chr, ignore.case = TRUE)
}

is_snp_acgt <- function(a) {
  a <- toupper(as.character(a))
  nchar(a) == 1 & a %chin% c("A", "C", "G", "T")
}

comp_allele <- function(a) {
  a <- toupper(as.character(a))
  map <- c(A = "T", T = "A", C = "G", G = "C")
  as.character(map[a])
}

target_markers <- unique(unlist(lapply(
  list.files(ref_dir, pattern = "[.]bim$", full.names = TRUE),
  function(path) fread(path, header = FALSE, select = 2)[[1]]
), use.names = FALSE))

marker_keep <- function(chr, pos, allele1, allele2) {
  chr <- normalize_chr(chr)
  pos <- as.character(pos)
  allele1 <- toupper(as.character(allele1))
  allele2 <- toupper(as.character(allele2))

  id_same <- paste0(chr, ":", pos, ":", allele1, ":", allele2)
  id_swap <- paste0(chr, ":", pos, ":", allele2, ":", allele1)

  id_comp_same <- rep(NA_character_, length(chr))
  id_comp_swap <- rep(NA_character_, length(chr))
  ok_comp <- is_snp_acgt(allele1) & is_snp_acgt(allele2)
  if (any(ok_comp)) {
    comp1 <- comp_allele(allele1[ok_comp])
    comp2 <- comp_allele(allele2[ok_comp])
    id_comp_same[ok_comp] <- paste0(chr[ok_comp], ":", pos[ok_comp], ":", comp1, ":", comp2)
    id_comp_swap[ok_comp] <- paste0(chr[ok_comp], ":", pos[ok_comp], ":", comp2, ":", comp1)
  }

  id_same %chin% target_markers |
    id_swap %chin% target_markers |
    id_comp_same %chin% target_markers |
    id_comp_swap %chin% target_markers
}

filter_bcac_overall <- function() {
  infile <- file.path(repo_data_dir, "bcac_2020_sumstats_for_s-mixcan_hg38.csv")
  outfile <- file.path(repo_data_dir, "bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv")

  gwas <- fread(infile)
  keep <- marker_keep(
    chr = gwas[["CHR"]],
    pos = gwas[["POS_hg38"]],
    allele1 = gwas[["Baseline.Meta"]],
    allele2 = gwas[["Effect.Meta"]]
  )

  fwrite(gwas[keep], outfile)
  cat(sprintf(
    "BCAC overall: kept %d of %d rows -> %s\n",
    sum(keep), nrow(gwas), outfile
  ))
}

filter_bcac_subtype <- function() {
  infile <- file.path(repo_data_dir, "bcac_2020_subtype_zscores_for_s-mixcan_hg38.txt")
  outfile <- file.path(repo_data_dir, "bcac_2020_subtype_zscores_for_s-mixcan_hg38_repro_subset.txt")

  gwas <- fread(infile)
  var_parts <- tstrsplit(gwas[["var_name"]], "_", fixed = TRUE)
  stopifnot(length(var_parts) >= 4)

  chr <- if ("CHR" %chin% names(gwas)) {
    gwas[["CHR"]]
  } else if ("chr" %chin% names(gwas)) {
    gwas[["chr"]]
  } else {
    var_parts[[1]]
  }

  keep <- marker_keep(
    chr = chr,
    pos = gwas[["POS_hg38"]],
    allele1 = var_parts[[3]],
    allele2 = var_parts[[4]]
  )

  fwrite(gwas[keep], outfile, sep = "\t")
  cat(sprintf(
    "BCAC subtype: kept %d of %d rows -> %s\n",
    sum(keep), nrow(gwas), outfile
  ))
}

cat(sprintf("Target reference markers: %d\n", length(target_markers)))
filter_bcac_overall()
filter_bcac_subtype()
