#!/usr/bin/env bash
set -euo pipefail

# --- CONFIG ---
# INPUT PGEN Directory containing reference genome: chr{CHR}_hg38.{pgen,psam,pvar} , SNP lists: gene_ID_snplist_new.txt
# Generate psam file in PGEN folder
#for i in {1..22} X Y; do
#    ln -s hg38_corrected.psam chr${i}_hg38.psam
#done

ID_DIR="/Users/zhusinan/Downloads/adriana/Filtered_Ref"
PGEN_DIR="/Users/zhusinan/Downloads/adriana/plink_snplist_by_gene"

# Output directory where  live and results will be written 
# results including Ref Genome: geneID_snplist_hg38_012.raw; 
#                   SNP info: geneID_snplist_hg38.bim
# OUT_DIR="/Users/zhusinan/Downloads/adriana/plink_snplist_by_gene"
OUT_DIR="/Users/zhusinan/Downloads/adriana/Result3/Ref"
# Threads
THREADS="4"
# EUR sample ID path
EUR_samples="/Users/zhusinan/Downloads/adriana/Result3/kgp_eur_unrelated_samples.txt"

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

# mkdir -p "${OUT_DIR}"

for i in "${!GENES[@]}"; do
  GENE="${GENES[$i]}"
  CHR="${CHRS[$i]}"
  echo "==> ${GENE} (chr${CHR})"

  SNPLIST="${ID_DIR}/${GENE}_snplist_new.txt"
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
  # (if conda is needed)
  conda activate plink2 
  plink2 --pfile "${PFILE}" \
         --extract "${SNPLIST}" \
         --make-bed \
         --keep "${EUR_samples}" \
         --threads "${THREADS}" \
         --out "${OUT_PREFIX}"

  # ---- PLINK1.9: recode A (012) ----
  # (if conda is needed)
  conda activate plink
  plink --bfile "${OUT_PREFIX}" \
        --recodeA \
        --threads "${THREADS}" \
        --out "${OUT_PREFIX}_012"

  # ---- PLINK1.9: LD r^2 matrix (square) ----
 # plink --bfile "${OUT_PREFIX}" \
 #       --r2 square \
 #       --threads "${THREADS}" \
 #       --out "${LD_PREFIX}"

  # Optional: compress LD matrix
  # gzip -f "${LD_PREFIX}.ld"

  echo "   [OK] Done ${GENE}"
done

echo "All done."
