# Reproducing the MD and BCAC TWAS Analyses

This document describes how to reproduce the MD and BCAC TWAS analyses reported in the manuscript.

The analysis scripts are in `SMiXcan_Analysis/Codes/`. The model weight files and per-gene 1000 Genomes reference files used by these scripts are included in `Data/`:

- `subset_weight_mixcan2_adipose.csv`
- `subset_weight_mixcan2_epithelial.csv`
- `subset_weight_mixcan2_fibroblast.csv`
- `subset_weight_predixcanlike.csv`
- `1000Genome_Ref/{gene_id}_eur_hg38.{bed,bim,fam}`
- `1000Genome_Ref/{gene_id}_eur_hg38_012.raw`
- `1000Genome_kgp_eur_unrelated_samples_ID.txt`

The scripts read these files using paths relative to the repository, so the weights and per-gene reference files do not need to be copied into the working analysis directory.

## Required External Input Files

The GWAS summary statistics are not bundled with the repository because they are external data resources. To reproduce the MD and BCAC analyses, create a local analysis directory with the following layout:

```text
analysis_dir/
├── MD_results_2021/
│   ├── DA_MetaResultswInfo_20210609_hg38.txt
│   ├── NDA_MetaResultswInfo_20210609_hg38.txt
│   └── PMD_MetaResultswInfo_20210609_hg38.txt
├── plink_snplist_by_gene/
│   └── bcac_2020_sumstats_for_s-mixcan_hg38.csv
└── bcac_2020_subtype_zscores_for_s-mixcan_hg38.txt
```

### Downloading the external GWAS files

The external files are not redistributed in this repository because they are large GWAS summary statistics from external consortium resources. Download them from the source pages below, then rename or preprocess them to match the filenames expected by the scripts.

Create the local folders first:

```bash
mkdir -p /path/to/analysis_dir/MD_results_2021
mkdir -p /path/to/analysis_dir/plink_snplist_by_gene
```

Then download and prepare the following files:

1. Mammographic density GWAS summary statistics

   Download the dense area (DA), non-dense area (NDA), and percent mammographic density (PMD) GWAS summary statistics from Chen et al., 2022, Breast Cancer Research 24:27, for the European ancestry replication population:

   - Chen 2022 mammographic density GWAS summary results: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results-1>
   - Article DOI: <https://doi.org/10.1186/s13058-022-01524-0>

   Download the original hg19 files from that page and place them in `analysis_dir/MD_results_2021/`:

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

   The liftover script keeps the original columns and adds `POS_hg38`, which is required by the MD TWAS scripts. It writes the final files as:

   ```text
   /path/to/analysis_dir/MD_results_2021/DA_MetaResultswInfo_20210609_hg38.txt
   /path/to/analysis_dir/MD_results_2021/NDA_MetaResultswInfo_20210609_hg38.txt
   /path/to/analysis_dir/MD_results_2021/PMD_MetaResultswInfo_20210609_hg38.txt
   ```

2. BCAC overall breast cancer GWAS summary statistics

   Download the 2020 overall breast cancer summary statistics from the BCAC summary-results page:

   - BCAC summary results: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results>

   Use the European ancestry overall breast cancer meta-analysis summary statistics for 133,384 cases and 113,789 controls. Prepare the file for S-MiXcan using hg38 coordinates and the columns expected by `SMiXcan_Analysis/Codes/runSMiXcan_bcac.R`, including `CHR`, `POS_hg38`, `Baseline.Meta`, `Effect.Meta`, `Beta.meta`, and `sdE.meta`. Save the final file as:

   ```text
   /path/to/analysis_dir/plink_snplist_by_gene/bcac_2020_sumstats_for_s-mixcan_hg38.csv
   ```

3. BCAC breast cancer subtype GWAS z-scores

   Download the breast cancer subtype summary statistics from the same BCAC summary-results page:

   - BCAC summary results: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results>

   Use the European ancestry subtype results for luminal A-like, luminal B-like, luminal B/HER2-negative-like, HER2-enriched-like, and triple-negative breast cancer. If the file is still in hg19, save the pre-liftover file as:

   ```text
   /path/to/analysis_dir/bcac_2020_subtype_zscores_for_s-mixcan.txt
   ```

   Then run:

   ```bash
   export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir
   python3 SMiXcan_Analysis/Codes/liftover_bcac_subtype.py
   ```

   This creates the hg38 file used by `runSMiXcan_subtype.R`:

   ```text
   /path/to/analysis_dir/bcac_2020_subtype_zscores_for_s-mixcan_hg38.txt
   ```

   The subtype analysis expects z-score columns for `Luminal_A`, `Luminal_B`, `Luminal_B_HER2Neg`, `HER2_Enriched`, and `Triple_Neg`.

The TWAS scripts use the bundled reference files in `Data/1000Genome_Ref/` by default. If you need to regenerate these reference files from genome-wide 1000 Genomes PLINK2 files, add the following external inputs to `analysis_dir/`:

```text
analysis_dir/
├── Filtered_Ref/
│   └── {gene_id}_snplist_new.txt
└── plink_snplist_by_gene/
    └── chr{CHR}_hg38.{pgen,pvar,psam}
```

Then run:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir
bash SMiXcan_Analysis/Codes/run_plink_ld.sh
```

By default, `run_plink_ld.sh` writes regenerated files to `Data/1000Genome_Ref/` and uses `Data/1000Genome_kgp_eur_unrelated_samples_ID.txt` as the EUR sample list. To write regenerated references elsewhere, set `MIXCAN_REF_DIR=/path/to/ref_dir`.

## Software Requirements

Install R and the R packages used by the scripts:

```r
install.packages(c(
  "data.table", "lme4", "glmnet", "doRNG", "tibble",
  "tidyr", "dplyr", "MASS", "devtools"
))
devtools::install_github("yaowuliu/ACAT")
devtools::install_github("songxiaoyu/MiXcan2/Package")
devtools::install_github("songxiaoyu/SMiXcan")
```

The analysis scripts also require PLINK2 and PLINK 1.9. One convenient option is to install both with conda:

```bash
conda create -n plink2 -c bioconda -c conda-forge plink plink2
conda activate plink2
plink2 --version
plink --version
```

Alternatively, download the binaries directly from the PLINK websites and add them to your `PATH`:

- PLINK 2: <https://www.cog-genomics.org/plink/2.0/>
- PLINK 1.9: <https://www.cog-genomics.org/plink/1.9/>

The helper script `SMiXcan_Analysis/Codes/run_plink_ld.sh` assumes `plink` and `plink2` are available after activating the `plink2` conda environment. If you use a different environment name or install location, update the `conda activate plink2` lines in that script or make sure both executables are already on your `PATH`.

For the MD GWAS liftover step and the optional BCAC subtype liftover step, install Python with `pandas` and `pyliftover`:

```bash
python3 -m pip install pandas pyliftover
python3 -c "import pandas, pyliftover; print('liftover dependencies OK')"
```

## Running the Analyses

After the external files have been downloaded and the required hg38 files have been created, set `MIXCAN_ANALYSIS_DIR` from the repository root and run:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir

Rscript SMiXcan_Analysis/Codes/runSMiXcan_MD.R
Rscript SMiXcan_Analysis/Codes/runPrediXcan_MD.R
Rscript SMiXcan_Analysis/Codes/runSMiXcan_bcac.R
Rscript SMiXcan_Analysis/Codes/runSMiXcan_subtype.R
```

For BCAC subtype data that still need to be lifted from hg19 to hg38, run the liftover script before `runSMiXcan_subtype.R`:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir
python3 SMiXcan_Analysis/Codes/liftover_bcac_subtype.py
```

## Expected Outputs

The scripts write output files under `analysis_dir/Result3/`:

- `Result_SMiXcan/`: MD S-MiXcan results for DA, NDA, and PMD.
- `Result_PrediXcan/`: MD PrediXcan-like comparison results.
- `Result_bcac/`: BCAC overall S-MiXcan and PrediXcan-like results.
- `Result_bcac_subtypes/`: BCAC subtype-specific S-MiXcan results.

The corresponding result tables used in the manuscript are included in this repository under `SMiXcan_Analysis/Result_*`.
