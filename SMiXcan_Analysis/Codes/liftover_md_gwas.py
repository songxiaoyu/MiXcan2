#!/usr/bin/env python3
"""Lift Chen 2022 mammographic density GWAS files from hg19 to hg38.

Expected input files under <analysis_dir>/MD_results_2021/:
  DA_MetaResultswInfo_20210609.txt
  NDA_MetaResultswInfo_20210609.txt
  PMD_MetaResultswInfo_20210609.txt

Outputs:
  DA_MetaResultswInfo_20210609_hg38.txt
  NDA_MetaResultswInfo_20210609_hg38.txt
  PMD_MetaResultswInfo_20210609_hg38.txt
"""

import argparse
import os
from pathlib import Path


PHENOTYPES = ("DA", "NDA", "PMD")
REQUIRED_COLUMNS = {
    "MarkerName",
    "rsid",
    "chr",
    "position",
    "Allele1",
    "Allele2",
    "freq",
    "info",
    "Zscore",
    "Weight",
    "P-value",
    "StudyNo",
    "HetPVal",
}

pd = None


def normalize_chr(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "":
        return None
    if text.lower().startswith("chr"):
        return text
    return f"chr{text}"


def lift_position(lo, chrom, position):
    chrom = normalize_chr(chrom)
    if chrom is None or pd.isna(position):
        return pd.NA

    try:
        pos = int(position)
    except (TypeError, ValueError):
        return pd.NA

    lifted = lo.convert_coordinate(chrom, pos)
    if not lifted:
        return pd.NA
    return int(lifted[0][1])


def validate_columns(columns, input_path):
    missing = sorted(REQUIRED_COLUMNS.difference(columns))
    if missing:
        missing_text = ", ".join(missing)
        raise ValueError(f"{input_path} is missing required columns: {missing_text}")


def liftover_file(lo, input_path, output_path, chunksize):
    print(f"Processing {input_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    total_rows = 0
    lifted_rows = 0
    first_chunk = True

    reader = pd.read_csv(
        input_path,
        sep=r"\s+",
        compression="infer",
        chunksize=chunksize,
        low_memory=False,
    )

    for chunk in reader:
        if first_chunk:
            validate_columns(set(chunk.columns), input_path)

        chunk["POS_hg38"] = [
            lift_position(lo, chrom, pos)
            for chrom, pos in zip(chunk["chr"], chunk["position"])
        ]

        total_rows += len(chunk)
        chunk = chunk.dropna(subset=["POS_hg38"]).copy()
        chunk["POS_hg38"] = chunk["POS_hg38"].astype("int64")
        lifted_rows += len(chunk)

        chunk.to_csv(
            output_path,
            sep="\t",
            index=False,
            mode="w" if first_chunk else "a",
            header=first_chunk,
            na_rep="",
        )
        first_chunk = False

    print(f"Saved {lifted_rows:,} of {total_rows:,} rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Lift Chen 2022 MD GWAS summary statistics from hg19 to hg38."
    )
    parser.add_argument(
        "--analysis-dir",
        default=os.environ.get("MIXCAN_ANALYSIS_DIR", "/Users/zhusinan/Downloads/adriana"),
        help="Analysis directory containing MD_results_2021/.",
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help="Directory containing original MD files. Defaults to <analysis-dir>/MD_results_2021.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory for hg38 output files. Defaults to the input directory.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=500000,
        help="Number of rows to process per chunk.",
    )
    args = parser.parse_args()

    global pd
    try:
        import pandas as pandas_module
        from pyliftover import LiftOver
    except ImportError as exc:
        raise SystemExit(
            "Missing Python dependency. Install dependencies with: "
            "python3 -m pip install pandas pyliftover"
        ) from exc

    pd = pandas_module

    analysis_dir = Path(args.analysis_dir).expanduser().resolve()
    input_dir = Path(args.input_dir).expanduser().resolve() if args.input_dir else analysis_dir / "MD_results_2021"
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else input_dir

    lo = LiftOver("hg19", "hg38")

    for pheno in PHENOTYPES:
        input_path = input_dir / f"{pheno}_MetaResultswInfo_20210609.txt"
        output_path = output_dir / f"{pheno}_MetaResultswInfo_20210609_hg38.txt"
        if not input_path.exists():
            raise FileNotFoundError(f"Missing input file: {input_path}")
        liftover_file(lo, input_path, output_path, args.chunksize)


if __name__ == "__main__":
    main()
