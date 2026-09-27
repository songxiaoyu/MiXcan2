suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  "scripts"
}
source(file.path(script_dir, "config.R"))

if (!file.exists(normalized_weights_file)) stop("Run step 1 first.", call. = FALSE)

weights <- fread(normalized_weights_file)
list_dir <- file.path(work_dir, "ld_reference_snplists")
dir.create(list_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- unique(weights[, .(gene_id, chr)])
manifest[, reference_exists :=
  file.exists(file.path(reference_dir, paste0(gene_id, "_eur_hg38.bim"))) &
  file.exists(file.path(reference_dir, paste0(gene_id, "_eur_hg38_012.raw")))]

for (gene in manifest[reference_exists == FALSE, gene_id]) {
  markers <- sort(unique(weights[gene_id == gene, MarkerName_hg38]))
  fwrite(
    data.table(marker = markers),
    file.path(list_dir, paste0(gene, "_snplist.txt")),
    col.names = FALSE
  )
}

fwrite(
  manifest[reference_exists == FALSE, .(gene_id, chr)],
  file.path(list_dir, "missing_reference_manifest.tsv"),
  sep = "\t"
)

cat("Prepared SNP lists for", manifest[reference_exists == FALSE, .N], "genes.\n")
print(manifest[reference_exists == FALSE])
