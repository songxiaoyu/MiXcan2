# Current workflow

This document describes the active S-MiXcan workflow using the September 2026
prediction weights. Historical Result2 analyses and earlier result directories
are intentionally excluded from the GitHub project.

## Method used for the final MD result

The final mammographic-density result uses:

- Result3 marker-ID matching;
- the 503-sample 1000 Genomes EUR reference;
- indel and palindromic-variant support;
- fixed regularization scale `0.001`;
- corrected MD effect/non-effect allele orientation;
- the September 2026 weights, including six genes absent from Result3.

Result3 itself had the wrong MD Z-score direction. Its matched SNPs and
two-sided p-values are used for validation, while the corrected workflow
reverses the Result3 Z direction.

## Directory layout

```text
SMixcan_analysis_260927/
  README.md
  docs/                         Current method and audit documentation
  scripts/                      Executable workflow, Steps 1-5
  inputs/                       September 2026 prediction weights
  data/                         Compact inputs and EUR LD references
  results/final_md/             Validated final MD result and diagnostics
  MD_SMiXcan_delivery_Sinan_20261113/  Files prepared for collaborators
```

Generated intermediate files are written to `work/` and
`results/current_run/`; both are ignored by Git.

## Inputs

Prediction weights:

```text
inputs/mixcan_predixcan_weights_v3_2026-09-21.csv
```

The repository includes compact DA, NDA, PMD, BCAC overall, and BCAC subtype
inputs as well as per-gene `.bim` and `_012.raw` EUR reference files for all 21
genes.

The Result3 validation benchmark is:

```text
data/reference/result3_md_wrong_direction_benchmark.csv
```

The benchmark is retained only to check matched SNP counts, p-values, and the
expected sign reversal. It must not be used as the corrected direction output.

## Run the complete workflow

From the repository root:

```bash
bash SMixcan_analysis_260927/scripts/run_all.sh
```

The runner performs:

1. weight validation and normalization;
2. reproducibility-input checks;
3. MD, BCAC overall, and BCAC subtype association analysis;
4. combination of association outputs;
5. construction and validation of the final Result3-compatible MD table.

The default association settings are fixed regularization with scale `0.001`.

## Run individual steps

```bash
Rscript SMixcan_analysis_260927/scripts/1_validate_weights.R
Rscript SMixcan_analysis_260927/scripts/2_download_reproducibility_inputs.R
Rscript SMixcan_analysis_260927/scripts/2b_prepare_ld_reference_lists.R
bash SMixcan_analysis_260927/scripts/2c_generate_missing_ld_reference.sh
python3 SMixcan_analysis_260927/scripts/2d_prepare_full_bcac.py
python3 SMixcan_analysis_260927/scripts/2e_prepare_full_md.py
Rscript SMixcan_analysis_260927/scripts/3_run_associations.R
Rscript SMixcan_analysis_260927/scripts/4_combine_results.R
Rscript SMixcan_analysis_260927/scripts/5_build_final_md_result.R
```

Steps 2b and 2c are only needed after adding genes without an LD reference.
Steps 2d and 2e are only needed when rebuilding compact inputs from complete
GWAS files.

## Outputs

Association outputs from a fresh run:

```text
results/current_run/
  md_results.csv
  bcac_overall_results.csv
  bcac_subtype_results.csv
  all_results.csv
  weight_coverage_manifest.csv
```

Validated final MD outputs:

```text
results/final_md/
  smixcan_result3_method_direction_corrected.csv
  result3_overlap_direction_validation.csv
  validation_summary.csv
```

The final MD table contains 177 rows: 123 Result3-overlap rows and 54 rows for
six newly added genes. It contains no PrediXcan-like rows and no
multiple-testing columns.

## Acceptance criteria

For all 123 Result3-overlap rows:

- matched SNP counts equal Result3;
- `p_join` equals Result3 within `1e-12`;
- corrected `Z_cell_1` equals negative Result3 `Z_1_join` within `1e-12`;
- corrected `Z_cell_2` equals negative Result3 `Z_2_join` within `1e-12`.

The Step 5 script enforces these checks and stops if any criterion fails.
