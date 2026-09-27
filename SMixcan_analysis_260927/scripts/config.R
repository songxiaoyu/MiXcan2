default_workflow_dir <- if (exists("script_dir", inherits = FALSE)) {
  file.path(script_dir, "..")
} else {
  "."
}

workflow_dir <- normalizePath(
  Sys.getenv(
    "MD_BCAC_WORKFLOW_DIR",
    unset = default_workflow_dir
  ),
  mustWork = FALSE
)

weights_file <- Sys.getenv(
  "MD_BCAC_WEIGHTS_FILE",
  unset = file.path(workflow_dir, "inputs", "mixcan_predixcan_weights_v3_2026-09-21.csv")
)

data_dir <- file.path(workflow_dir, "data")
reference_dir <- file.path(data_dir, "1000Genome_Ref")
legacy_reference_dir <- Sys.getenv(
  "MD_BCAC_LEGACY_REFERENCE_DIR",
  unset = ""
)
work_dir <- file.path(workflow_dir, "work")
results_dir <- file.path(workflow_dir, "results_sinan_20260927")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

normalized_weights_file <- file.path(work_dir, "weights_normalized.csv")
weights_manifest_file <- file.path(results_dir, "weight_coverage_manifest.csv")

github_raw_base <- "https://raw.githubusercontent.com/songxiaoyu/MiXcan2/main/Data"

md_repro_files <- c(
  DA = "DA_MetaResultswInfo_20210609_hg38_repro_subset.txt",
  NDA = "NDA_MetaResultswInfo_20210609_hg38_repro_subset.txt",
  PMD = "PMD_MetaResultswInfo_20210609_hg38_repro_subset.txt"
)
md_full_weight_files <- c(
  DA = file.path("full_md", "DA_MetaResultswInfo_20210609_hg38_weight_subset.csv"),
  NDA = file.path("full_md", "NDA_MetaResultswInfo_20210609_hg38_weight_subset.csv"),
  PMD = file.path("full_md", "PMD_MetaResultswInfo_20210609_hg38_weight_subset.csv")
)
md_files <- if (all(file.exists(file.path(data_dir, md_full_weight_files)))) {
  md_full_weight_files
} else {
  md_repro_files
}

full_bcac_dir <- file.path(data_dir, "full_gwas")
bcac_overall_full_subset <- file.path(
  full_bcac_dir,
  "bcac_2020_sumstats_for_s-mixcan_hg38_weight_subset.csv"
)
bcac_subtype_full_subset <- file.path(
  full_bcac_dir,
  "bcac_2020_subtype_zscores_for_s-mixcan_hg38_weight_subset.txt"
)

bcac_overall_file <- Sys.getenv(
  "MD_BCAC_OVERALL_FILE",
  unset = if (file.exists(bcac_overall_full_subset)) {
    bcac_overall_full_subset
  } else {
    file.path(data_dir, "bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv")
  }
)
bcac_subtype_file <- Sys.getenv(
  "MD_BCAC_SUBTYPE_FILE",
  unset = if (file.exists(bcac_subtype_full_subset)) {
    bcac_subtype_full_subset
  } else {
    file.path(data_dir, "bcac_2020_subtype_zscores_for_s-mixcan_hg38_repro_subset.txt")
  }
)

bcac_n_case <- 106278
bcac_n_control <- 91477

regularization_mode <- Sys.getenv("SMIXCAN_REGULARIZATION", unset = "fixed")
regularization_scale <- as.numeric(Sys.getenv("SMIXCAN_REG_SCALE", unset = "0.1"))
condition_max <- as.numeric(Sys.getenv("SMIXCAN_CONDITION_MAX", unset = "100"))
try_download_missing_reference <- tolower(Sys.getenv(
  "MD_BCAC_TRY_DOWNLOAD_MISSING_REF",
  unset = "true"
)) %in% c("true", "1", "yes")

if (!regularization_mode %in% c("fixed", "estimate")) {
  stop("SMIXCAN_REGULARIZATION must be fixed or estimate.", call. = FALSE)
}
if (!is.finite(regularization_scale) || regularization_scale < 0) {
  stop("SMIXCAN_REG_SCALE must be a non-negative number.", call. = FALSE)
}
