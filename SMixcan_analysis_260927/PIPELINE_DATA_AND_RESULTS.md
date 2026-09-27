# Pipeline, data, and results

This document describes the independent S-MiXcan analysis using the updated
MiXcan and PrediXcan-like weights. It records the complete workflow, input data
locations, generated files, and the status of the latest run.

## Workflow directory

```text
Dropbox/Paper_MiXcan2/GitHub/SMixcan_analysis_260927
```

This workflow is separate from the original publication analysis under
`SMiXcan_Analysis/` and does not overwrite the original code or results.

All executable workflow code and `config.R` are stored under:

```text
SMixcan_analysis_260927/scripts/
```

## Analyses

The workflow tests the updated prediction weights against:

1. Mammographic density: DA, NDA, and PMD.
2. BCAC overall breast cancer risk.
3. BCAC subtypes: Luminal A, Luminal B, Luminal B HER2-negative,
   HER2 enriched, and triple negative.

The weights include three MiXcan cell-type models and one PrediXcan-like model:

| Model | Gene models |
|---|---:|
| Epithelial MiXcan | 19 |
| Fibroblast MiXcan | 21 |
| Adipose MiXcan | 19 |
| PrediXcan-like | 20 |
| Total | 79 |

## Input data

### Prediction weights

```text
inputs/mixcan_predixcan_weights_v3_2026-09-21.csv
```

The file contains 21 genes, 79 gene-model combinations, and 931 unique hg38
variants. Step 1 writes the standardized version to:

```text
work/weights_normalized.csv
```

### Mammographic-density GWAS

The complete hg38 DA, NDA, and PMD summary statistics are stored outside the
workflow because each file is approximately 1 GB:

```text
/Users/zhusinan/PhD_Projects/adriana/MD_results_2021/
  DA_MetaResultswInfo_20210609_hg38.txt
  NDA_MetaResultswInfo_20210609_hg38.txt
  PMD_MetaResultswInfo_20210609_hg38.txt
```

Step 2e extracts only the current weight loci and stores the compact inputs at:

```text
data/full_md/
  DA_MetaResultswInfo_20210609_hg38_weight_subset.csv
  NDA_MetaResultswInfo_20210609_hg38_weight_subset.csv
  PMD_MetaResultswInfo_20210609_hg38_weight_subset.csv
```

On another computer, set the full-data directory before running:

```bash
export MD_BCAC_MD_FULL_DIR=/path/to/MD_results_2021
```

### BCAC GWAS

The complete downloaded hg19 files are stored at:

```text
data/full_gwas/
  bcac_2020_sumstats_for_s-mixcan.txt.gz
  bcac_2020_subtype_zscores_for_s-mixcan.txt.gz
```

Step 2d uses `data/liftover/hg38ToHg19.over.chain.gz` to map the current hg38
weight loci to the hg19 BCAC coordinates. It writes compact analysis inputs:

```text
data/full_gwas/
  bcac_2020_sumstats_for_s-mixcan_hg38_weight_subset.csv
  bcac_2020_subtype_zscores_for_s-mixcan_hg38_weight_subset.txt
```

### LD reference

The hg38 1000 Genomes EUR per-gene reference files are stored under:

```text
data/1000Genome_Ref/
```

Each gene uses a PLINK BIM file for variant annotation and a corresponding
`_012.raw` dosage file. References are available for all 21 current genes.

## Pipeline steps

### Step 1: validate and normalize weights

```bash
Rscript scripts/1_validate_weights.R
```

Input:

```text
inputs/mixcan_predixcan_weights_v3_2026-09-21.csv
```

Outputs:

```text
work/weights_normalized.csv
results_sinan_20260927/weight_coverage_manifest.csv
```

### Step 2: obtain public reproducibility inputs

```bash
Rscript scripts/2_download_reproducibility_inputs.R
```

This provides fallback MD subsets and the public LD-reference files. The full
compact inputs are preferred whenever steps 2d and 2e have been completed.

### Steps 2b and 2c: generate missing LD references

```bash
Rscript scripts/2b_prepare_ld_reference_lists.R
bash scripts/2c_generate_missing_ld_reference.sh
```

These steps are required only after adding genes without a matching reference
under `data/1000Genome_Ref/`.

### Step 2d: prepare complete BCAC inputs

```bash
python3 scripts/2d_prepare_full_bcac.py
```

This step requires `pyliftover`. It extracts the current weight SNPs from the
complete BCAC files and validates available coordinates against the publication
subset.

### Step 2e: prepare complete MD inputs

```bash
python3 scripts/2e_prepare_full_md.py
```

This streams the three complete MD files without loading them fully into
memory. The original files are not modified or copied.

### Step 3: association analysis

```bash
Rscript scripts/3_run_associations.R
```

Inputs:

- Normalized prediction weights.
- Compact full MD and BCAC summary statistics.
- Per-gene 1000 Genomes EUR LD references.

Outputs:

```text
results_sinan_20260927/md_results.csv
results_sinan_20260927/bcac_overall_results.csv
results_sinan_20260927/bcac_subtype_results.csv
```

`results_sinan_20260927/skipped_models.csv` is created only when a model is skipped or fails.
The current complete run has no skipped models, so this file is absent.

The three result tables contain columns only through `p_join`. Diagnostic and
multiple-testing columns are not included in the collaborator-facing outputs.

### Step 4: combine final tables

```bash
Rscript scripts/4_combine_results.R
```

Outputs:

```text
results_sinan_20260927/all_results.csv
```

This step combines the three association tables without applying FDR,
Bonferroni, or any other multiple-testing correction. The collaborator will
perform the required correction downstream.

## Run the full workflow

From the workflow directory:

```bash
export MD_BCAC_PYTHON=python3
export MD_BCAC_MD_FULL_DIR=/path/to/MD_results_2021
bash scripts/run_all.sh
```

The runner automatically prepares full BCAC and MD compact inputs when their
source files are available, then runs association testing and combines the
results.

## Latest completed run

Run date: 2026-09-27

| Analysis | Completed | Skipped or failed |
|---|---:|---:|
| MD DA | 79/79 | 0 |
| MD NDA | 79/79 | 0 |
| MD PMD | 79/79 | 0 |
| BCAC overall | 79/79 | 0 |
| BCAC Luminal A | 79/79 | 0 |
| BCAC Luminal B | 79/79 | 0 |
| BCAC Triple negative | 79/79 | 0 |
| BCAC HER2 enriched | 79/79 | 0 |
| BCAC Luminal B HER2-negative | 79/79 | 0 |

All 711 requested results completed and all 711 have finite `p_join` values.

## Final result locations

The primary result for downstream review is:

```text
results_sinan_20260927/all_results.csv
```

Raw per-analysis output remains in `md_results.csv`,
`bcac_overall_results.csv`, and `bcac_subtype_results.csv`.

## Git and large files

The workflow code, weights, compact analysis-ready data, LD references, and
final results are located inside the Paper_MiXcan2 Git repository. The local
`.gitignore` excludes `work/` and the two complete raw BCAC GWAS archives.
Those raw archives exceed GitHub's normal per-file limit and can be restored
from the source listed above. Generated intermediates and the incomplete
`.partial` download remain local under `work/`.
