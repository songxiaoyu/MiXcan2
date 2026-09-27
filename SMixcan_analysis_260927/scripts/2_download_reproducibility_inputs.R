suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  "scripts"
}
source(file.path(script_dir, "config.R"))

if (!file.exists(normalized_weights_file)) {
  stop("Run 1_validate_weights.R first.", call. = FALSE)
}

download_one <- function(url, destination, required = TRUE) {
  if (file.exists(destination) && file.info(destination)$size > 0) {
    return(TRUE)
  }
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temp_file <- tempfile(tmpdir = dirname(destination))
  status <- tryCatch({
    suppressWarnings(download.file(url, temp_file, mode = "wb", quiet = TRUE))
    file.exists(temp_file) && file.info(temp_file)$size > 0
  }, error = function(e) FALSE)
  if (isTRUE(status)) {
    invisible(file.rename(temp_file, destination))
    return(TRUE)
  }
  unlink(temp_file)
  if (required) stop("Could not download required input: ", url, call. = FALSE)
  FALSE
}

for (phenotype in names(md_files)) {
  filename <- md_files[[phenotype]]
  download_one(
    paste(github_raw_base, "MD_results_2021", filename, sep = "/"),
    file.path(data_dir, filename)
  )
}

download_one(
  paste(github_raw_base, basename(bcac_overall_file), sep = "/"),
  bcac_overall_file
)
download_one(
  paste(github_raw_base, basename(bcac_subtype_file), sep = "/"),
  bcac_subtype_file
)

weights <- fread(normalized_weights_file)
genes <- sort(unique(weights$gene_id))
reference_status <- rbindlist(lapply(genes, function(gene) {
  bim_name <- paste0(gene, "_eur_hg38.bim")
  raw_name <- paste0(gene, "_eur_hg38_012.raw")
  bim_destination <- file.path(reference_dir, bim_name)
  raw_destination <- file.path(reference_dir, raw_name)

  if (nzchar(legacy_reference_dir)) {
    legacy_bim <- file.path(legacy_reference_dir, bim_name)
    legacy_raw <- file.path(legacy_reference_dir, raw_name)
    if (!file.exists(bim_destination) && file.exists(legacy_bim)) {
      file.copy(legacy_bim, bim_destination)
    }
    if (!file.exists(raw_destination) && file.exists(legacy_raw)) {
      file.copy(legacy_raw, raw_destination)
    }
  }

  bim_ok <- file.exists(bim_destination)
  raw_ok <- file.exists(raw_destination)
  if (try_download_missing_reference && !bim_ok) {
    bim_ok <- download_one(
      paste(github_raw_base, "1000Genome_Ref", bim_name, sep = "/"),
      bim_destination,
      required = FALSE
    )
  }
  if (try_download_missing_reference && !raw_ok) {
    raw_ok <- download_one(
      paste(github_raw_base, "1000Genome_Ref", raw_name, sep = "/"),
      raw_destination,
      required = FALSE
    )
  }
  data.table(gene_id = gene, reference_available = bim_ok && raw_ok)
}))

cat("Downloaded public reproducibility inputs.\n")
print(reference_status[, .N, by = reference_available])
