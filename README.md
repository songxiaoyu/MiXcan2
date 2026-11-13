# MiXcan2 and S-MiXcan Analysis

This repository contains the MiXcan2 R package, the general MiXcan2
individual-level workflow, and the latest S-MiXcan mammographic-density
analysis.

## Repository Structure

Only the current code and compact reproducibility files are kept in GitHub:

| Path | Purpose |
| --- | --- |
| [`Package/`](Package/) | Source code for the `MiXcan2` R package. |
| [`scripts/`](scripts/) | General MiXcan2 model-training, prediction, and association scripts. |
| [`SMixcan_analysis_260927/`](SMixcan_analysis_260927/) | Latest S-MiXcan workflow using the September 2026 weights, including compact inputs, corrected MD results, validation files, and documentation. |
| [`REPRODUCING_RESULTS.md`](REPRODUCING_RESULTS.md) | Entry point for reproducing the current workflows. |

Older publication snapshots, superseded outputs, raw GWAS files, generated
PLINK artifacts, and IDE files are intentionally excluded from the active
repository.

The project-owner archive is stored outside GitHub at:

```text
Dropbox/Paper_MiXcan2/Archive/Paper_MiXcan2_local_archive_20261113/
```

## Install MiXcan2

Install the package directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("songxiaoyu/MiXcan2", subdir = "Package")
library(MiXcan2)
```

## Current S-MiXcan Workflow

The maintained summary-statistics workflow is documented in:

- [`SMixcan_analysis_260927/README.md`](SMixcan_analysis_260927/README.md)
- [`SMixcan_analysis_260927/docs/WORKFLOW.md`](SMixcan_analysis_260927/docs/WORKFLOW.md)

The dated collaborator delivery is:

- [`SMixcan_analysis_260927/MD_SMiXcan_delivery_Sinan_20261113/`](SMixcan_analysis_260927/MD_SMiXcan_delivery_Sinan_20261113/)

The final MD table uses the Result3 marker-matching strategy with 503-sample
European LD references and indel support, while correcting the Result3
Z-score direction.

## Publications

- MiXcan: <https://www.nature.com/articles/s41467-023-35888-4>
- Cell-type-aware mammographic-density TWAS preprint:
  <https://www.medrxiv.org/content/10.1101/2025.11.25.25341027v1>
- S-MiXcan preprint:
  <https://www.medrxiv.org/content/10.64898/2026.03.22.26349035v1>
