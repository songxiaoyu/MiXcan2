# Workflow for Reproducing the Published Results

## 1. Reproducing MiXcan2 (Discovery) Results 

The full pipeline used in discovery analyses has two stages: (1) training the
MiXcan2 ensemble prediction models, and (2) running the TWAS (prediction
+ association) using those trained models. Scripts for both stages are
included in the `scripts/` folder.

**1. Train the ensemble models.** Please note that the genomics data are protected due to participant privacy and data governance requirements. To access these data, please refer to the data access request procedures. 
Once you have access to the GTEx and GWAS data (see the preprint for a full
description of the training data and study population), run
[`scripts/run_ensemble_mixcan2.R`](scripts/run_ensemble_mixcan2.R) to build MiXcan2
models trained using GTEx v8 mammary tissue samples from female subjects
with European ancestry.

**2. Run the TWAS.** With trained models in hand, the TWAS itself is run
in two steps, also in `scripts/`:

- Step 1 (prediction): [`scripts/predict_expression.R`](scripts/predict_expression.R)
  applies the trained MiXcan2 weights to genotype data to predict
  cell-type-level gene expression (GReX) for each sample.
- Step 2 (association): [`scripts/association_analysis.R`](scripts/association_analysis.R)
  tests the predicted expression from Step 1 for association with the
  phenotype(s) of interest.


# 2. Reproducing S-MiXcan (Validation & Enrichment) Results

We applied S-MiXcan to  (1) validate significant genes in  MD TWAS using Chen 2022 GWAS summary statistics 
and (2) evaluate their association with breast cancer risk using BCAC GWAS summary statistics.

As this analysis only requires prediction weights from discovery analysis and a subset of GWAS summary
statistics, which can be made publicly available. We uploaded data together with step-by-step codes for 
reproducible analysis. The details are as follows: 

### Software Installation

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

### Input Data

#### 1. Model Weight Files

The model weight files from MiXcan2 used by these scripts are included in `Data/`:

- `subset_weight_mixcan2_adipose.csv`
- `subset_weight_mixcan2_epithelial.csv`
- `subset_weight_mixcan2_fibroblast.csv`
- `subset_weight_predixcanlike.csv`

#### 2. 1000 Genomes Reference Files

The per-gene 1000 Genomes reference files used by these scripts are included in `Data/`:

- `1000Genome_Ref/{gene_id}_eur_hg38.{bed,bim,fam}`
- `1000Genome_Ref/{gene_id}_eur_hg38_012.raw`
- `1000Genome_kgp_eur_unrelated_samples_ID.txt`

#### 3. GWAS Summary Statistics (Chen 2022) for MD

The hg38 mammographic density GWAS subsets required to reproduce these TWAS results are included under `Data/`:

```text
Data/MD_results_2021/DA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/NDA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/PMD_MetaResultswInfo_20210609_hg38_repro_subset.txt
```
More details about the generation of this dataset is in the later section. 

#### 4. BCAC Overall Breast Cancer GWAS Summary Statistics
The hg38 BCAC overall breast cancer GWAS subset required to reproduce these TWAS results is included under `Data/`:

```text
Data/bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv
```
More details about the generation of this dataset is in the later section. 

#### 5. BCAC Breast Cancer Subtype GWAS Z-Scores

The hg38 subtype z-score subset required to reproduce these TWAS results is included under `Data/`:

```text
Data/bcac_2020_subtype_zscores_for_s-mixcan_hg38_repro_subset.txt
```
Download the breast cancer subtype summary statistics from the same BCAC summary-results page:

- BCAC summary results: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results>
- The analysis focuses on significant genes and uses hg38 coordinates. See below for details.

Use the European ancestry subtype results for luminal A-like, luminal B-like, luminal B/HER2-negative-like, HER2-enriched-like, and triple-negative breast cancer. 


The subtype analysis expects z-score columns for `Luminal_A`, `Luminal_B`, `Luminal_B_HER2Neg`, `HER2_Enriched`, and `Triple_Neg`.

### Running the Analyses

Set `MIXCAN_ANALYSIS_DIR` from the repository root and run:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir

Rscript SMiXcan_Analysis/Codes/runSMiXcan_MD.R
Rscript SMiXcan_Analysis/Codes/runPrediXcan_MD.R
Rscript SMiXcan_Analysis/Codes/runSMiXcan_bcac.R
Rscript SMiXcan_Analysis/Codes/runSMiXcan_subtype.R
```

### Expected Outputs

The scripts write output files under `analysis_dir/Result/`:

- `Result_SMiXcan/`: MD S-MiXcan results for DA, NDA, and PMD.
- `Result_PrediXcan/`: MD PrediXcan-like comparison results.
- `Result_bcac/`: BCAC overall S-MiXcan and PrediXcan-like results.
- `Result_bcac_subtypes/`: BCAC subtype-specific S-MiXcan results.

The corresponding result tables used in the manuscript are included in this repository under `SMiXcan_Analysis/Result_*`.


### Additional Information of the Input Data (Not Needed for Reproducing Results)

#### Creation of the MD GWAS Subsets
<!--The repository includes the MD rows that match the bundled model weights and 1000 Genomes reference files. To recreate these subset files from the full Chen 2022 mammographic density GWAS files, download the original hg19 files and place them in :
-->


- Link: <https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results-1>

- The analysis focuses on significant genes and uses hg38 coordinates. See below for details.


 To create the input data, download the full Chen 2022 mammographic density GWAS files from 
 [here](https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results/gwas-summary-results-1) 
 and save as following:
```text
/path/to/analysis_dir/MD_results_2021/DA_MetaResultswInfo_20210609.txt
/path/to/analysis_dir/MD_results_2021/NDA_MetaResultswInfo_20210609.txt
/path/to/analysis_dir/MD_results_2021/PMD_MetaResultswInfo_20210609.txt
```
Article DOI is <https://doi.org/10.1186/s13058-022-01524-0>. These files include  dense area (DA), 
non-dense area (NDA), and percent mammographic density (PMD) GWAS summary statistics from Chen et al., 2022, 
Breast Cancer Research 24:27, for the European ancestry replication population:

The summary statistics is based on hg19. First, lift them from hg19 to hg38:

Install Python with `pandas` and `pyliftover` before running the liftover script:

```bash
python3 -m pip install pandas pyliftover
python3 -c "import pandas, pyliftover; print('liftover dependencies OK')"
```

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir
python3 SMiXcan_Analysis/Codes/liftover_md_gwas.py
```

The liftover script saves the results as:

```text
/path/to/analysis_dir/MD_results_2021/DA_MetaResultswInfo_20210609_hg38.txt
/path/to/analysis_dir/MD_results_2021/NDA_MetaResultswInfo_20210609_hg38.txt
/path/to/analysis_dir/MD_results_2021/PMD_MetaResultswInfo_20210609_hg38.txt
```

Then, subset to significant genes:

```bash
Rscript SMiXcan_Analysis/Codes/make_repro_md_gwas_subset.R /path/to/analysis_dir/MD_results_2021
```

Then we get the input data:

```text
Data/MD_results_2021/DA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/NDA_MetaResultswInfo_20210609_hg38_repro_subset.txt
Data/MD_results_2021/PMD_MetaResultswInfo_20210609_hg38_repro_subset.txt
```

#### Creation of the Input Reference LD

Download reference genome from https://www.cog-genomics.org/plink/2.0/resources. 
Use 2022-08-04 Byrska-Bishop et al. (build 38, 3202 samples, contigs unphased). 
Please pick‘Split by chromosome’, ‘info annotations’, ‘King-based pedigree corrections. Save files in `analysis_dir/reference_pgen/`. 

The reference-regeneration script requires PLINK2. One convenient option is to install it with conda:

```bash
conda create -n plink2 -c bioconda -c conda-forge plink2
conda activate plink2
plink2 --version
```

Alternatively, download the PLINK2 binary directly and add it to your `PATH`:

- PLINK 2: <https://www.cog-genomics.org/plink/2.0/>

Get the snps used in analysis here in `analysis_dir/gene_snplists/`

Get the EUR sample ID here: `Data/1000Genome_kgp_eur_unrelated_samples_ID.txt`

Then run:

```bash
export MIXCAN_ANALYSIS_DIR=/path/to/analysis_dir
bash SMiXcan_Analysis/Codes/run_plink_ld.sh
```

The output will be LD matrix for each gene, and saved in `Data/1000Genome_Ref/`.

### Creation of the Input BCAC GWAS Subsets




To create the input data, download the orginal input data BCAC 2020 GWAS files from 
[here](https://www.ccge.medschl.cam.ac.uk/breast-cancer-association-consortium-bcac/data-data-access/summary-results). 
Use the European ancestry overall breast cancer meta-analysis summary statistics for 133,384 cases and 113,789 controls. 
Keep these columns for BC and its subtype analyses:

- Overall BC: `CHR`, `POS_hg38`, `Baseline.Meta`, `Effect.Meta`, `Beta.meta`, and `sdE.meta`.
- BC subtype: `var_name`, `POS_hg38`, and z-score columns for `Luminal_A`, `Luminal_B`, `Luminal_B_HER2Neg`, `HER2_Enriched`, and `Triple_Neg`.

Then, prepare the  hg38 files from hg19 as above and save the results as:

```text
Data/bcac_2020_sumstats_for_s-mixcan_hg38.csv
Data/bcac_2020_subtype_zscores_for_s-mixcan_hg38.txt
```

Then subset to the significant gene related subset as input using:

```bash
Rscript SMiXcan_Analysis/Codes/make_repro_gwas_subset.R
```

The output will be in:

```text
Data/bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv
Data/bcac_2020_subtype_zscores_for_s-mixcan_hg38_repro_subset.txt
```
