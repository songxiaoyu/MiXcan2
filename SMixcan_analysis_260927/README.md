# S-MiXcan analysis with updated prediction weights

This directory contains the current, reviewable S-MiXcan workflow for the
September 2026 prediction weights. Historical Result2 compatibility code,
earlier result directories, caches, and raw large downloads are not part of
the active GitHub project.

## Final collaborator delivery

Submit the two files under:

```text
MD_SMiXcan_delivery_Sinan_20261113/
```

The final MD result uses:

- Result3 allele-ID matching and the 503-sample EUR LD reference;
- indel and palindromic-variant support;
- fixed regularization scale `0.001`;
- corrected MD effect/non-effect allele orientation;
- September 2026 weights, including six genes absent from Result3.

Result3 had the wrong MD Z-score direction. For 123 overlapping rows, the
final table reproduces its matched SNP counts and two-sided p-values while
reversing both cell-component Z scores to the corrected direction.

Final composition:

| Status | Genes | Models | Rows |
|---|---:|---:|---:|
| Result3 overlap | 15 | 41 | 123 |
| New genes | 6 | 18 | 54 |
| Total | 21 | 59 | 177 |

The final table excludes PrediXcan-like rows and does not apply
multiple-testing correction.

## Project layout

```text
SMixcan_analysis_260927/
  README.md
  docs/
    WORKFLOW.md
    MD_HARMONIZATION_AUDIT.md
    RESULT3_DIRECTION_CORRECTED_WORKFLOW.md
  scripts/
    1_validate_weights.R
    2_download_reproducibility_inputs.R
    2b_prepare_ld_reference_lists.R
    2c_generate_missing_ld_reference.sh
    2d_prepare_full_bcac.py
    2e_prepare_full_md.py
    3_run_associations.R
    4_combine_results.R
    5_build_final_md_result.R
    config.R
    run_all.sh
  inputs/
  data/
  results/final_md/
  MD_SMiXcan_delivery_Sinan_20261113/
```

## Run

From the repository root:

```bash
bash SMixcan_analysis_260927/scripts/run_all.sh
```

The default regularization is fixed scale `0.001`. Intermediate outputs are
written to ignored directories under `work/` and `results/current_run/`.
The validated MD table is written to `results/final_md/`.

To rebuild only the final MD table from the included corrected input and
Result3 benchmark:

```bash
Rscript SMixcan_analysis_260927/scripts/5_build_final_md_result.R
```

Step 5 performs all row-count, SNP-count, p-value, and direction checks and
stops if validation fails.

## Documentation

- [Current workflow](docs/WORKFLOW.md)
- [MD harmonization audit](docs/MD_HARMONIZATION_AUDIT.md)
- [Result3 direction-corrected result](docs/RESULT3_DIRECTION_CORRECTED_WORKFLOW.md)

## Inputs retained in Git

- September 2026 weights;
- compact MD and BCAC analysis inputs;
- per-gene 1000 Genomes EUR `.bim` and `_012.raw` files;
- Result3 validation benchmark;
- corrected `0.001` MD input used to build the final table.

Complete raw GWAS downloads and generated caches are intentionally excluded.
They are not required to reproduce the included final MD table.
