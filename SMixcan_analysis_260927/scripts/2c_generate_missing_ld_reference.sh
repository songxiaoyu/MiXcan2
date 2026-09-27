#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PGEN_DIR="${MD_BCAC_PGEN_DIR:-${WORKFLOW_DIR}/data/plink_input}"
EUR_SAMPLES="${MD_BCAC_EUR_SAMPLES:-${PGEN_DIR}/eur_ids.txt}"
PSAM_FILE="${MD_BCAC_PSAM_FILE:-${PGEN_DIR}/hg38_corrected.psam}"
PLINK2_BIN="${PLINK2_BIN:-/opt/anaconda3/envs/plink2/bin/plink2}"
OUTPUT_DIR="${WORKFLOW_DIR}/data/1000Genome_Ref"
LIST_DIR="${WORKFLOW_DIR}/work/ld_reference_snplists"
PFILE_LINK_DIR="${WORKFLOW_DIR}/work/plink_pfile_links"
MANIFEST="${LIST_DIR}/missing_reference_manifest.tsv"

Rscript "${SCRIPT_DIR}/2b_prepare_ld_reference_lists.R"

check_local_file() {
  local path="$1"
  if [[ ! -s "${path}" ]]; then
    echo "ERROR: missing or empty input: ${path}" >&2
    exit 1
  fi
  if [[ "$(uname -s)" == "Darwin" ]]; then
    local blocks
    blocks="$(stat -f '%b' "${path}")"
    if [[ "${blocks}" == "0" ]]; then
      echo "ERROR: ${path} is a cloud-only placeholder (dataless)." >&2
      echo "Download it locally first or set MD_BCAC_PGEN_DIR to a local copy." >&2
      exit 1
    fi
  fi
}

check_local_file "${EUR_SAMPLES}"
check_local_file "${PSAM_FILE}"
if [[ ! -x "${PLINK2_BIN}" ]]; then
  echo "ERROR: plink2 is not executable: ${PLINK2_BIN}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${PFILE_LINK_DIR}"

while IFS=$'\t' read -r gene chr; do
  [[ -z "${gene}" ]] && continue

  source_prefix="${PGEN_DIR}/chr${chr}_hg38"
  pgen_prefix="${PFILE_LINK_DIR}/chr${chr}_hg38"
  snplist="${LIST_DIR}/${gene}_snplist.txt"
  output_prefix="${OUTPUT_DIR}/${gene}_eur_hg38"

  check_local_file "${source_prefix}.pgen"
  check_local_file "${source_prefix}.pvar"
  check_local_file "${snplist}"

  # The chromosome-specific .psam files in this dataset are empty stubs.
  # Build a temporary PLINK prefix that reuses the shared corrected sample file.
  ln -sfn "${source_prefix}.pgen" "${pgen_prefix}.pgen"
  ln -sfn "${source_prefix}.pvar" "${pgen_prefix}.pvar"
  ln -sfn "${PSAM_FILE}" "${pgen_prefix}.psam"

  echo "Generating EUR reference for ${gene} on chr${chr}"
  "${PLINK2_BIN}" \
    --pfile "${pgen_prefix}" \
    --extract "${snplist}" \
    --keep "${EUR_SAMPLES}" \
    --make-bed \
    --out "${output_prefix}"

  "${PLINK2_BIN}" \
    --bfile "${output_prefix}" \
    --export A \
    --out "${output_prefix}_012"

  if [[ ! -s "${output_prefix}.bim" || ! -s "${output_prefix}_012.raw" ]]; then
    echo "ERROR: reference generation failed for ${gene}" >&2
    exit 1
  fi
done < <(tail -n +2 "${MANIFEST}")

echo "All available missing references were generated in ${OUTPUT_DIR}."
