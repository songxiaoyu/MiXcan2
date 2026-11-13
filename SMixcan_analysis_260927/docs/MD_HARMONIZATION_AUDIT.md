# Mammographic density allele-harmonization audit

## Finding

The September 2026 mammographic density (MD) analysis passed the MD allele
columns to the reference harmonizer in the wrong order. In the MD summary
statistics, `Allele1` is the effect allele for `Zscore` and `Allele2` is the
non-effect allele. The PLINK marker/reference dosage is oriented to the second
allele in marker IDs such as `1:11548365:G:A_A`. Therefore MD must be supplied
to the harmonizer as non-effect/effect: `Allele2`, then `Allele1`.

The previous call used `Allele1`, then `Allele2`. This preserved the p-values
but reversed the reported directions of association. BCAC is not affected:
its `Baseline.Meta`/`Effect.Meta` allele order was already passed correctly.

## Correction

The MD call in `scripts/3_run_associations.R` now uses:

```r
a1_col = "Allele2"
a2_col = "Allele1"
```

The superseded publication-era workflow was moved out of the active GitHub
tree during repository cleanup. All maintained reruns should use
`scripts/3_run_associations.R` in this project.

## Rerun verification

The September workflow was rerun after the correction.

- MD: 237 of 237 rows completed.
- Every `Z_cell_1` and `Z_cell_2` is exactly the negative of its pre-correction
  value.
- Every `p_cell_1`, `p_cell_2`, and `p_join` is unchanged.
- The BCAC overall and subtype tables are unchanged.
- On the 123 rows shared with the archived manuscript-correction table,
  Pearson correlations are 0.9814 for `Z_cell_1` and 0.9809 for `Z_cell_2`;
  directions agree for 114/123 and 117/123 rows, respectively.

The remaining differences from the archived manuscript results are expected
because the September analysis uses a different prediction-weight/model set,
different SNP overlap for many models, and a different association
regularization setting. A direction-only allele correction cannot change a
two-sided p-value.

## Corrected files retained in the current project

- `data/derived/md_results_corrected_reg0p001_all_models.csv`
- `results/final_md/smixcan_result3_method_direction_corrected.csv`

The former dated result directories were removed from the active GitHub
project during the latest-only reorganization. They are not required to
rebuild or validate the final MD result.

The compact benchmark retained for direction validation is:

```text
data/reference/result3_md_wrong_direction_benchmark.csv
```

The original publication snapshot is stored only in the dated local archive;
it is not required for the maintained workflow.

## Result3-compatible final result

The December Result3 table was subsequently confirmed to have the wrong MD
direction. For continuity with its newer marker-ID harmonization, indel
support, and 503-sample EUR reference, the final collaborator result uses the
Result3 method with fixed regularization `0.001`, but reverses the Result3 Z
directions through the corrected effect/non-effect allele order.

Across the 123 rows shared with the September weights, matched SNP counts and
two-sided p-values reproduce Result3, while both cell-component Z scores are
the negatives of Result3. Six genes newly available in the September weights
are added with the same corrected workflow.

Current delivery:

```text
MD_SMiXcan_delivery_Sinan_20261113/
```
