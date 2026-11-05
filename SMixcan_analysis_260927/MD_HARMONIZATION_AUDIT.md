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

The publication-era script `SMiXcan_Analysis/Codes/runSMiXcan_MD.R` has the
same correction to prevent future reruns from reintroducing the error.

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

## Corrected files

- `results_sinan_20260927/md_results.csv`
- `results_sinan_20260927/all_results.csv`

The archived manuscript-correction table used to verify direction is:

```text
SMiXcan_Analysis/Archive/Result2/Result_SMiXcan/SMiXcan_result_merged.csv
```

Local file history identifies this as the corrected result used by the
manuscript. It does not establish whether the file was delivered by email,
Google Drive, or Dropbox in January/February 2026.
