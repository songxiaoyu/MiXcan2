# MD and BCAC analysis with updated prediction weights

This folder is an independent, reviewable rerun of the published mammographic
density (MD) and BCAC S-MiXcan analyses. It does not modify the original
`MiXcan2/SMiXcan_Analysis` code or results.

For the complete step-by-step workflow, input locations, output locations, and
latest run status, see [`PIPELINE_DATA_AND_RESULTS.md`](PIPELINE_DATA_AND_RESULTS.md).
All executable workflow files are organized under `scripts/`.

## Scope

The workflow analyzes the updated weights against:

1. Mammographic density phenotypes: `DA`, `NDA`, and `PMD`.
2. BCAC overall breast cancer risk.
3. BCAC subtypes: `Luminal_A`, `Luminal_B`, `Luminal_B_HER2Neg`,
   `HER2_Enriched`, and `Triple_Neg`.

The 1000 Genomes EUR reference files and public reproducibility subsets are
downloaded from `songxiaoyu/MiXcan2`. Complete MD, BCAC overall, and BCAC
subtype files are then reduced to the SNPs required by the current weights.
Genes without a complete LD reference are reported in
`results_sinan_20260927/skipped_models.csv` and are not silently treated as null results.

## Input weights

Default input, stored inside this independent workflow:

```text
SMixcan_analysis_260927/inputs/mixcan_predixcan_weights_v3_2026-09-21.csv
```

Required columns:

```text
cell_type,yName,xNameMatrix,type,weight_cell_1,weight_cell_2
```

Override the path on another computer with:

```bash
export MD_BCAC_WEIGHTS_FILE=/path/to/new_weights.csv
```

## Data sources

The previous BCAC files were retrieved from the shared Google Drive folder:

```text
https://drive.google.com/drive/folders/1qbXbnuO11oU_0HBZQIVaXR5aMf-92NeQ
```

Complete downloads are stored locally under `data/full_gwas/`:

- `bcac_2020_sumstats_for_s-mixcan.txt.gz`
- `bcac_2020_subtype_zscores_for_s-mixcan.txt.gz`

These files use hg19 coordinates. Step 2d maps only the 931 variants used by
the current weights to hg19, extracts those records, and writes compact inputs
with verified hg38 positions. It validates the chain mappings against the
existing publication subset before writing any analysis input. The two complete
raw archives exceed GitHub's normal per-file limit and are therefore excluded
from Git; their compact analysis-ready subsets are included in the repository.
The interrupted
download of the larger raw subtype file is retained under
`work/download_cache/` with a `.partial` suffix and must not be used as input.

The complete hg38 MODE/MD summary statistics are stored locally under:

```text
/Users/zhusinan/PhD_Projects/adriana/MD_results_2021
```

Step 2e streams the complete DA, NDA, and PMD files and writes compact inputs
containing the loci used by the current weights. On another computer, specify
the source directory with:

```bash
export MD_BCAC_MD_FULL_DIR=/path/to/MD_results_2021
```

Step 2d requires Python 3 and `pyliftover`:

```bash
python3 -m pip install pyliftover
```

## Run one step at a time

```bash
Rscript SMixcan_analysis_260927/scripts/1_validate_weights.R
Rscript SMixcan_analysis_260927/scripts/2_download_reproducibility_inputs.R
# Optional when updated weights contain genes absent from the public reference:
Rscript SMixcan_analysis_260927/scripts/2b_prepare_ld_reference_lists.R
bash SMixcan_analysis_260927/scripts/2c_generate_missing_ld_reference.sh
python3 SMixcan_analysis_260927/scripts/2d_prepare_full_bcac.py
python3 SMixcan_analysis_260927/scripts/2e_prepare_full_md.py
Rscript SMixcan_analysis_260927/scripts/3_run_associations.R
Rscript SMixcan_analysis_260927/scripts/4_combine_results.R
```

Or run the complete workflow:

```bash
export MD_BCAC_PYTHON=python3
bash SMixcan_analysis_260927/scripts/run_all.sh
```

`scripts/run_all.sh` runs step 2d automatically when both full BCAC files and the
liftOver chain are present. It runs step 2e when the complete MD directory is
available. The association step prefers the compact full inputs and falls back
to the public reproducibility subsets only when they are absent. Set
`MD_BCAC_OVERALL_FILE` or `MD_BCAC_SUBTYPE_FILE` to override the BCAC inputs.

## Outputs

Intermediate files are written to `work/`. Final files are written directly to
`results_sinan_20260927/`:

- `md_results.csv`
- `bcac_overall_results.csv`
- `bcac_subtype_results.csv`
- `all_results.csv`
- `weight_coverage_manifest.csv`

`skipped_models.csv` is created only when at least one model is skipped or
fails. It is absent from the current complete run.

The association tables contain only the model identifiers, SNP counts,
cell-component joint statistics, and `p_join`. They end at the `p_join` column.
Multiple-testing correction is intentionally not performed in this workflow.

## Statistical settings

- MD uses `family = "gaussian"`.
- BCAC overall and subtype analyses use `family = "binomial"`.
- BCAC counts are preserved from the publication workflow:
  `n_case = 106278`, `n_control = 91477`.
- Regularization defaults to fixed scale `0.1`, matching the current package
  default. Override with `SMIXCAN_REG_SCALE`.
- Alleles are harmonized by same, swapped, strand-complement, and
  complement-plus-swapped matching. Effect estimates are sign-flipped for
  swapped matches.

## Reproducibility limitation

The two complete raw BCAC GWAS archives are large and are not tracked by Git.
All compact analysis-ready inputs under `data/full_md/` and `data/full_gwas/`,
the liftOver chain, the public reproducibility subsets, and the per-gene LD
references are tracked. The compact inputs are specific to the current weights.
After replacing the weights, rerun steps 1, 2d, and 2e before association
testing.

The required per-gene references for the current 21 genes are already stored in
`data/1000Genome_Ref/`. Optional step 2c is only needed after adding new genes.
It can generate additional references from chromosome PLINK2 files using these
default paths inside the workflow:

```text
data/plink_input/chr{1..22}_hg38.{pgen,pvar}
data/plink_input/hg38_corrected.psam
data/plink_input/eur_ids.txt
```

For PLINK files stored outside the workflow, set:

```bash
export MD_BCAC_PGEN_DIR=/local/path/to/plink_snplist_by_gene
export MD_BCAC_EUR_SAMPLES=/local/path/to/eur_ids.txt
export MD_BCAC_PSAM_FILE=/local/path/to/hg38_corrected.psam
```

On macOS, step 2c rejects Dropbox placeholders whose contents have not been
downloaded locally.
