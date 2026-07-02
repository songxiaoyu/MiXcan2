#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
# INPUT folders:
#   gene_snplists/: {gene_id}_snplist_new.txt
#   reference_pgen/: chr{CHR}_hg38.{pgen,psam,pvar}
# Generate psam file in PGEN folder
#for i in {1..22} X Y; do
#    ln -s hg38_corrected.psam chr${i}_hg38.psam
#done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DATA_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)/Data"

ANALYSIS_DIR="${MIXCAN_ANALYSIS_DIR:-/Users/zhusinan/Downloads/adriana}"
SNPLIST_DIR="${MIXCAN_SNPLIST_DIR:-${ANALYSIS_DIR}/gene_snplists}"
PGEN_DIR="${MIXCAN_PGEN_DIR:-${ANALYSIS_DIR}/reference_pgen}"

# Output directory where reference files will be written.
# Results include Ref Genome: geneID_snplist_hg38_012.raw
#                 SNP info: geneID_snplist_hg38.bim
# OUT_DIR="/Users/zhusinan/Downloads/adriana/1000Genome_Ref"
OUT_DIR="${MIXCAN_REF_DIR:-${REPO_DATA_DIR}/1000Genome_Ref}"
# Threads
THREADS="4"
# EUR sample ID path
EUR_samples="${MIXCAN_EUR_SAMPLES:-${REPO_DATA_DIR}/1000Genome_kgp_eur_unrelated_samples_ID.txt}"

#optional (if conda is needed)
source /opt/anaconda3/etc/profile.d/conda.sh

# --- INPUT LISTS (parallel arrays; same order and length) ---
GENES=(
"ENSG00000134917.9"
"ENSG00000165269.12"
"ENSG00000110092.3"
"ENSG00000257800.1"
"ENSG00000176920.11"
"ENSG00000166333.13"
"ENSG00000236036.1"
"ENSG00000277147.4"
"ENSG00000084676.15"
"ENSG00000142233.11"
"ENSG00000224397.5"
"ENSG00000198914.3"
"ENSG00000232871.8"
"ENSG00000159164.9"
"ENSG00000092607.13"
"ENSG00000226445.1"
"ENSG00000048707.13"
"ENSG00000149054.15"
"ENSG00000103199.13"
"ENSG00000183779.6"
)

CHRS=(
"11"
"9"
"11"
"2"
"19"
"11"
"13"
"1"
"2"
"19"
"20"
"2"
"19"
"1"
"1"
"6"
"1"
"11"
"16"
"8"
)


if [[ ${#GENES[@]} -ne ${#CHRS[@]} ]]; then
  echo "[ERROR] GENES and CHRS must have the same length."
  exit 1
fi

mkdir -p "${OUT_DIR}"

for i in "${!GENES[@]}"; do
  GENE="${GENES[$i]}"
  CHR="${CHRS[$i]}"
  echo "==> ${GENE} (chr${CHR})"

  SNPLIST="${SNPLIST_DIR}/${GENE}_snplist_new.txt"
  OUT_PREFIX="${OUT_DIR}/${GENE}_eur_hg38"
  PFILE="${PGEN_DIR}/chr${CHR}_hg38"

  # Basic existence checks
  if [[ ! -s "${SNPLIST}" ]]; then
    echo "   [WARN] Missing snplist: ${SNPLIST} — skipping."
    continue
  fi
  if [[ ! -s "${PFILE}.pgen" || ! -s "${PFILE}.pvar" || ! -s "${PFILE}.psam" ]]; then
    echo "   [WARN] Missing pfile set for chr${CHR}: ${PFILE}.{pgen,pvar.zst,psam} — skipping."
    continue
  fi

  # ---- PLINK2: extract variants & make BED/BIM/FAM ----
  conda activate plink2
  plink2 --pfile "${PFILE}" \
         --extract "${SNPLIST}" \
         --make-bed \
         --keep "${EUR_samples}" \
         --threads "${THREADS}" \
         --out "${OUT_PREFIX}"

  # ---- PLINK2: export additive genotype counts (.raw) ----
  plink2 --bfile "${OUT_PREFIX}" \
         --export A \
         --threads "${THREADS}" \
         --out "${OUT_PREFIX}_012"

  # ---- Optional: LD r^2 matrix (square) ----
 # plink2 --bfile "${OUT_PREFIX}" \
 #        --r2 square \
 #        --threads "${THREADS}" \
 #        --out "${LD_PREFIX}"

  # Optional: compress LD matrix
  # gzip -f "${LD_PREFIX}.ld"

  echo "   [OK] Done ${GENE}"
done

echo "All done."
