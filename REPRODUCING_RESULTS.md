# Reproducing S-MiXcan Results

This document describes how to reproduce:

1. S-MiXcan based MD TWAS using Chen 2022 GWAS summary statistics .
2. S-MiXcan based breast cancer risk TWAS using BCAC GWAS summary statistics.

## Software Requirements

### R Packages

Install R and the R packages used by the scripts:

```r
install.packages(c(
  "data.table", "lme4", "glmnet", "doRNG", "tibble",
  "tidyr", "dplyr", "MASS", "devtools"
))

devtools::install_github("petraf01/BayesDeBulk/BayesDeBulk")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("limma")

install.packages("matrixStats")
install.packages("nnls")
install.packages("R.methodsS3")
install.packages("lcmix", repos = "http://r-forge.r-project.org")
install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)

devtools::install_github("yaowuliu/ACAT")
devtools::install_github("kjgleason/Primo")
devtools::install_github("songxiaoyu/MiXcan2/Package")
devtools::install_github("songxiaoyu/SMiXcan")
```

On macOS, if package compilation fails, install the command line tools and gfortran:

```bash
sudo xcode-select --install
open https://mac.r-project.org/tools/gfortran-14.2-universal.pkg
```

### PLINK

The optional reference-regeneration script requires PLINK2. One convenient option is to install it with conda:

```bash
conda create -n plink2 -c bioconda -c conda-forge plink2
conda activate plink2
plink2 --version
```

Alternatively, download the PLINK2 binary directly and add it to your `PATH`:

- PLINK 2: <https://www.cog-genomics.org/plink/2.0/>

<!--The helper script `SMiXcan_Analysis/Codes/run_plink_ld.sh` assumes `plink2` is available after activating the `plink2` conda environment. If you use a different environment name or install location, update the `conda activate plink2` lines in that script or make sure the executable is already on your `PATH`.-->

### Python Liftover Dependencies

As we need to liftover SNPs between hg19 to hg38, install Python with `pandas` and `pyliftover`:

```bash
python3 -m pip install pandas pyliftover
python3 -c "import pandas, pyliftover; print('liftover dependencies OK')"
```

## Input Data

### 1. Model Weight Files

The model weight files from MiXcan2 used by these scripts are included in `Data/`:

- `subset_weight_mixcan2_adipose.csv`
- `subset_weight_mixcan2_epithelial.csv`
- `subset_weight_mixcan2_fibroblast.csv`
- `subset_weight_predixcanlike.csv`

### 2. 1000 Genomes Reference Files

The per-gene 1000 Genomes reference files used by these scripts are included in `Data/`:

- `1000Genome_Ref/{gene_id}_eur_hg38.{bed,bim,fam}`
- `1000Genome_Ref/{gene_id}_eur_hg38_012.raw`
- `1000Genome_kgp_eur_unrelated_samples_ID.txt`

### 3. GWAS Summary Statistics (Chen 2022) for MD

The hg38 mammographic density GWAS subsets required to reproduce these TWAS results are included under `Data/`:

```text
Data/MD_results_2021/DA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/NDA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/PMD_MetaResultswInfo_20210609_hg38_repro_subset.txt
```

These subset files were created from the dense area (DA), non-dense area (NDA), and percent mammographic density (PMD) GWAS summary statistics from Chen et al., 2022, Breast Cancer Research 24:27, for the European ancestry replication population:

- Link: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results-1>
- Article DOI: <https://doi.org/10.1186/s13058-022-01524-0>

### 4. BCAC Overall Breast Cancer GWAS Summary Statistics

Download the 2020 overall breast cancer summary statistics from the BCAC summary-results page:

- BCAC summary results: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results>

Use the European ancestry overall breast cancer meta-analysis summary statistics for 133,384 cases and 113,789 controls. The hg38 BCAC overall breast cancer GWAS subset required to reproduce these TWAS results is included under `Data/`:

```text
Data/bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv
```

### 5. BCAC Breast Cancer Subtype GWAS Z-Scores

Download the breast cancer subtype summary statistics from the same BCAC summary-results page:

- BCAC summary results: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results>

Use the European ancestry subtype results for luminal A-like, luminal B-like, luminal B/HER2-negative-like, HER2-enriched-like, and triple-negative breast cancer. The hg38 subtype z-score subset required to reproduce these TWAS results is included under `Data/`:

```text
Data/bcac_2020_subtype_zscores_for_s-mixcan_hg38_repro_subset.txt
```

The subtype analysis expects z-score columns for `Luminal_A`, `Luminal_B`, `Luminal_B_HER2Neg`, `HER2_Enriched`, and `Triple_Neg`.

## Running the Analyses

Set `MIXCAN_ANALYSIS_DIR` from the repository root and run:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir

Rscript SMiXcan_Analysis/Codes/runSMiXcan_MD.R
Rscript SMiXcan_Analysis/Codes/runPrediXcan_MD.R
Rscript SMiXcan_Analysis/Codes/runSMiXcan_bcac.R
Rscript SMiXcan_Analysis/Codes/runSMiXcan_subtype.R
```

## Expected Outputs

The scripts write output files under `analysis_dir/Result/`:

- `Result_SMiXcan/`: MD S-MiXcan results for DA, NDA, and PMD.
- `Result_PrediXcan/`: MD PrediXcan-like comparison results.
- `Result_bcac/`: BCAC overall S-MiXcan and PrediXcan-like results.
- `Result_bcac_subtypes/`: BCAC subtype-specific S-MiXcan results.

The corresponding result tables used in the manuscript are included in this repository under `SMiXcan_Analysis/Result_*`.

The analysis scripts are in `SMiXcan_Analysis/Codes/`.

## Optional: Recreate the MD GWAS Subsets

The repository includes the MD rows that match the bundled model weights and 1000 Genomes reference files. To recreate these subset files from the full Chen 2022 mammographic density GWAS files, download the original hg19 files and place them in `analysis_dir/MD_results_2021/`:

```text
/path/to/analysis_dir/MD_results_2021/DA_MetaResultswInfo_20210609.txt
/path/to/analysis_dir/MD_results_2021/NDA_MetaResultswInfo_20210609.txt
/path/to/analysis_dir/MD_results_2021/PMD_MetaResultswInfo_20210609.txt
```

Then lift them from hg19 to hg38:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir
python3 SMiXcan_Analysis/Codes/liftover_md_gwas.py
```

The liftover script keeps the original columns and adds `POS_hg38`. It writes:

```text
/path/to/analysis_dir/MD_results_2021/DA_MetaResultswInfo_20210609_hg38.txt
/path/to/analysis_dir/MD_results_2021/NDA_MetaResultswInfo_20210609_hg38.txt
/path/to/analysis_dir/MD_results_2021/PMD_MetaResultswInfo_20210609_hg38.txt
```

Create the reproducibility subsets:

```bash
Rscript SMiXcan_Analysis/Codes/make_repro_md_gwas_subset.R /path/to/analysis_dir/MD_results_2021
```

This writes:

```text
Data/MD_results_2021/DA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/NDA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/PMD_MetaResultswInfo_20210609_hg38_repro_subset.txt
```

## Optional Reference Regeneration

The TWAS scripts use the bundled reference files in `Data/1000Genome_Ref/` by default. If you need to regenerate these reference files from genome-wide 1000 Genomes PLINK2 files, add the following external inputs to `analysis_dir/`:

```text
analysis_dir/
├── gene_snplists/
│   └── {gene_id}_snplist_new.txt
└── reference_pgen/
    └── chr{CHR}_hg38.{pgen,pvar,psam}
```

Then run:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir
bash SMiXcan_Analysis/Codes/run_plink_ld.sh
```

By default, `run_plink_ld.sh` reads SNP-list files from `analysis_dir/gene_snplists/`, reads genome-wide PLINK2 files from `analysis_dir/reference_pgen/`, writes regenerated files to `Data/1000Genome_Ref/`, and uses `Data/1000Genome_kgp_eur_unrelated_samples_ID.txt` as the EUR sample list. To use different input folders, set `MIXCAN_SNPLIST_DIR=/path/to/gene_snplists` or `MIXCAN_PGEN_DIR=/path/to/reference_pgen`. To write regenerated references elsewhere, set `MIXCAN_REF_DIR=/path/to/ref_dir`.

## Optional: Recreate the Prepared BCAC GWAS Subsets

The repository includes the BCAC rows that match the bundled model weights and 1000 Genomes reference files. To recreate these subset files from the full prepared BCAC GWAS files, download the 2020 overall breast cancer and subtype summary statistics from the BCAC summary-results page:

- BCAC summary results: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results>

Prepare the full hg38 files using the columns expected by the analysis scripts:

- Overall BCAC: `CHR`, `POS_hg38`, `Baseline.Meta`, `Effect.Meta`, `Beta.meta`, and `sdE.meta`.
- BCAC subtype: `var_name`, `POS_hg38`, and z-score columns for `Luminal_A`, `Luminal_B`, `Luminal_B_HER2Neg`, `HER2_Enriched`, and `Triple_Neg`.

Save the full prepared files as:

```text
Data/bcac_2020_sumstats_for_s-mixcan_hg38.csv
Data/bcac_2020_subtype_zscores_for_s-mixcan_hg38.txt
```

Then create the reproducibility subsets:

```bash
Rscript SMiXcan_Analysis/Codes/make_repro_gwas_subset.R
```

This writes:

```text
Data/bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv
Data/bcac_2020_subtype_zscores_for_s-mixcan_hg38_repro_subset.txt
```
