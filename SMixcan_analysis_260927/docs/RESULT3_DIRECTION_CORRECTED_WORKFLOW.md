# Result3-method results with corrected MD direction

## Purpose

This workflow is the collaborator-facing continuation of
`Dropbox_from_Sinan_20251218_Result3 2` using the September 2026 weights.

The December Result3 table has the wrong MD association direction. Its SNP
matching and two-sided p-values are retained as the historical benchmark, but
all cell-component Z scores are corrected to the effect-allele direction.

The workflow uses:

- the September 2026 weights table
  `inputs/mixcan_predixcan_weights_v3_2026-09-21.csv`;
- the Result3 allele-ID matching behavior and 503-sample EUR reference;
- fixed regularization scale `0.001`;
- corrected MD effect/non-effect allele order;
- six newly available genes in addition to Result3-overlap genes.

## Build the result

From the repository root:

```bash
Rscript \
  SMixcan_analysis_260927/scripts/5_build_final_md_result.R
```

The validation benchmark is stored in the repository at:

```text
SMixcan_analysis_260927/data/reference/
result3_md_wrong_direction_benchmark.csv
```

Set `SMIXCAN_RESULT3_MD_FILE` to use another copy of that file.

The corrected association input is:

```text
SMixcan_analysis_260927/data/derived/
md_results_corrected_reg0p001_all_models.csv
```

## Primary output

```text
results/final_md/smixcan_result3_method_direction_corrected.csv
```

It contains only S-MiXcan rows and ends at `p_join`; no multiple-testing
columns are added.

Expected composition:

| Status | Genes | Models | Rows across DA/NDA/PMD |
|---|---:|---:|---:|
| `result3_overlap` | 15 | 41 | 123 |
| `new_gene` | 6 | 18 | 54 |
| Total | 21 | 59 | 177 |

Five Result3 genes are absent from the September weights and are not copied
into the new result.

The six genes newly added from the September weights are:

```text
ENSG00000068137.14
ENSG00000100393.12
ENSG00000105325.13
ENSG00000140400.16
ENSG00000163634.11
ENSG00000224238.2
```

## Validation meaning

For all 123 Result3-overlap rows:

- matched SNP counts must equal Result3;
- `p_join` must equal Result3 within `1e-12`;
- corrected `Z_cell_1` must equal `-Result3 Z_1_join` within `1e-12`;
- corrected `Z_cell_2` must equal `-Result3 Z_2_join` within `1e-12`.

The sign reversal is intentional. Result3 used the wrong MD allele direction;
two-sided p-values are unchanged when the Z-score sign is reversed.

Diagnostic outputs:

- `results/final_md/result3_overlap_direction_validation.csv`
- `results/final_md/validation_summary.csv`
