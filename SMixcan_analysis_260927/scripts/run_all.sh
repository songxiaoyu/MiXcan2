#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

Rscript "${SCRIPT_DIR}/1_validate_weights.R"
Rscript "${SCRIPT_DIR}/2_download_reproducibility_inputs.R"

PYTHON_BIN="${MD_BCAC_PYTHON:-python3}"
FULL_BCAC_DIR="${WORKFLOW_DIR}/data/full_gwas"
if [[ -f "${FULL_BCAC_DIR}/bcac_2020_sumstats_for_s-mixcan.txt.gz" && \
      -f "${FULL_BCAC_DIR}/bcac_2020_subtype_zscores_for_s-mixcan.txt.gz" && \
      -f "${WORKFLOW_DIR}/data/liftover/hg38ToHg19.over.chain.gz" ]]; then
  "${PYTHON_BIN}" "${SCRIPT_DIR}/2d_prepare_full_bcac.py"
fi

MD_FULL_DIR="${MD_BCAC_MD_FULL_DIR:-/Users/zhusinan/PhD_Projects/adriana/MD_results_2021}"
if [[ -f "${MD_FULL_DIR}/DA_MetaResultswInfo_20210609_hg38.txt" && \
      -f "${MD_FULL_DIR}/NDA_MetaResultswInfo_20210609_hg38.txt" && \
      -f "${MD_FULL_DIR}/PMD_MetaResultswInfo_20210609_hg38.txt" ]]; then
  MD_BCAC_MD_FULL_DIR="${MD_FULL_DIR}" \
    "${PYTHON_BIN}" "${SCRIPT_DIR}/2e_prepare_full_md.py"
fi

Rscript "${SCRIPT_DIR}/3_run_associations.R"
Rscript "${SCRIPT_DIR}/4_combine_results.R"
SMIXCAN_CORRECTED_MD_FILE="${WORKFLOW_DIR}/results/current_run/md_results.csv" \
  Rscript "${SCRIPT_DIR}/5_build_final_md_result.R"

echo "Association results: ${WORKFLOW_DIR}/results/current_run"
echo "Final validated MD result: ${WORKFLOW_DIR}/results/final_md"
