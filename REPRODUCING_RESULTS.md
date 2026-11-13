# Reproducing the Current Results

This document describes the two maintained workflows in this repository.
Legacy publication snapshots and superseded result directories are not part of
the active GitHub tree.

## 1. Software

Install R and the MiXcan2 package:

```r
install.packages("remotes")
remotes::install_github("songxiaoyu/MiXcan2", subdir = "Package")
library(MiXcan2)
```

The S-MiXcan workflow also requires the `SMiXcan` package and the R packages
listed in [`SMixcan_analysis_260927/README.md`](SMixcan_analysis_260927/README.md).

## 2. MiXcan2 Individual-Level Workflow

The general scripts are retained under [`scripts/`](scripts/):

1. [`run_ensemble_mixcan2.R`](scripts/run_ensemble_mixcan2.R): train MiXcan2 ensemble prediction models.
2. [`predict_expression.R`](scripts/predict_expression.R): predict cell-type-level GReX in application data.
3. [`association_analysis.R`](scripts/association_analysis.R): test predicted GReX for phenotype association.

These scripts require controlled individual-level genotype, expression,
cell-fraction, covariate, and phenotype data. Those protected datasets are not
distributed through GitHub. Set the input and output paths inside each script
before running it.

## 3. Current S-MiXcan Summary-Statistics Workflow

The complete maintained workflow is under
[`SMixcan_analysis_260927/`](SMixcan_analysis_260927/). Start with its
[`README.md`](SMixcan_analysis_260927/README.md) for input definitions,
configuration, commands, and output descriptions.

Run the full pipeline from the repository root with:

```bash
bash SMixcan_analysis_260927/scripts/run_all.sh
```

The pipeline performs:

1. weight validation;
2. reproducibility-input preparation;
3. LD-reference list preparation and optional reference generation;
4. GWAS preparation and association testing;
5. result combination; and
6. construction and validation of the corrected Result3-compatible MD table.

## 4. Final MD Result

The final reviewable delivery is stored directly in:

```text
SMixcan_analysis_260927/MD_SMiXcan_delivery_Sinan_20261113/
```

It contains:

```text
smixcan_result3_method_direction_corrected_20261113.csv
RESULT3_DIRECTION_CORRECTED_REPRODUCTION_GUIDE_20261113.md
```

Supporting validation files are in:

```text
SMixcan_analysis_260927/results/final_md/
```

The final method preserves the newer Result3 SNP matching, 503-sample European
LD reference, and indel support, while correcting the MD Z-score direction and
adding results for genes present only in the September 2026 weight table.

## 5. Archived Material

Large raw GWAS inputs, superseded outputs, generated PLINK files, and the older
publication snapshot are maintained outside the active GitHub repository. They
are not required to review the final delivery, and adding them back would make
the repository ambiguous and unnecessarily large.
