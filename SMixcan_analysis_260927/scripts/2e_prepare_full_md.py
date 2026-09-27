#!/usr/bin/env python3
"""Extract the current prediction-weight loci from complete hg38 MD GWAS files."""

from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
WEIGHTS_FILE = WORKFLOW_DIR / "work" / "weights_normalized.csv"
DEFAULT_SOURCE_DIR = Path("/Users/zhusinan/PhD_Projects/adriana/MD_results_2021")
SOURCE_DIR = Path(os.environ.get("MD_BCAC_MD_FULL_DIR", DEFAULT_SOURCE_DIR)).expanduser()
OUTPUT_DIR = WORKFLOW_DIR / "data" / "full_md"
PHENOTYPES = ("DA", "NDA", "PMD")


def normalize_chr(value: str) -> str:
    value = value.strip()
    return value[3:] if value.lower().startswith("chr") else value


def parse_position(value: str) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def complement(allele: str) -> str | None:
    allele = allele.upper()
    if len(allele) != 1 or allele not in "ACGT":
        return None
    return allele.translate(str.maketrans("ACGT", "TGCA"))


def alleles_compatible(ref: str, alt: str, allele1: str, allele2: str) -> bool:
    ref = ref.upper()
    alt = alt.upper()
    allele1 = allele1.upper()
    allele2 = allele2.upper()
    if (allele1, allele2) in ((ref, alt), (alt, ref)):
        return True
    comp1 = complement(allele1)
    comp2 = complement(allele2)
    return (
        comp1 is not None
        and comp2 is not None
        and (comp1, comp2) in ((ref, alt), (alt, ref))
    )


def read_markers() -> dict[str, dict[str, object]]:
    if not WEIGHTS_FILE.is_file():
        raise SystemExit(f"Missing normalized weights: {WEIGHTS_FILE}. Run step 1 first.")
    markers: dict[str, dict[str, object]] = {}
    with WEIGHTS_FILE.open(newline="") as handle:
        for row in csv.DictReader(handle):
            marker = row["MarkerName_hg38"]
            chrom, pos, ref, alt = marker.split(":", 3)
            markers.setdefault(
                marker,
                {
                    "marker": marker,
                    "chr": normalize_chr(chrom),
                    "pos": int(pos),
                    "ref": ref.upper(),
                    "alt": alt.upper(),
                },
            )
    return markers


def extract_phenotype(
    phenotype: str, markers: dict[str, dict[str, object]]
) -> tuple[int, int]:
    source = SOURCE_DIR / f"{phenotype}_MetaResultswInfo_20210609_hg38.txt"
    if not source.is_file():
        raise SystemExit(f"Missing full MD input: {source}")

    output = OUTPUT_DIR / (
        f"{phenotype}_MetaResultswInfo_20210609_hg38_weight_subset.csv"
    )
    targets: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for item in markers.values():
        targets[(str(item["chr"]), int(item["pos"]))].append(item)

    extracted = 0
    compatible = 0
    with source.open(newline="", errors="replace") as source_handle, output.open(
        "w", newline=""
    ) as output_handle:
        reader = csv.DictReader(source_handle)
        if reader.fieldnames is None:
            raise SystemExit(f"No header found in {source}")
        required = {"chr", "POS_hg38", "Allele1", "Allele2", "Zscore", "P-value"}
        missing = sorted(required.difference(reader.fieldnames))
        if missing:
            raise SystemExit(f"Missing columns in {source}: {', '.join(missing)}")

        output_fields = list(reader.fieldnames) + ["MarkerName_hg38_target"]
        writer = csv.DictWriter(output_handle, fieldnames=output_fields)
        writer.writeheader()

        for row in reader:
            pos = parse_position(row["POS_hg38"])
            if pos is None:
                continue
            candidates = targets.get((normalize_chr(row["chr"]), pos))
            if not candidates:
                continue
            for item in candidates:
                current = dict(row)
                current["MarkerName_hg38_target"] = item["marker"]
                writer.writerow(current)
                extracted += 1
                count_name = f"{phenotype}_rows"
                item[count_name] = int(item.get(count_name, 0)) + 1
                if alleles_compatible(
                    str(item["ref"]),
                    str(item["alt"]),
                    row["Allele1"],
                    row["Allele2"],
                ):
                    compatible += 1
                    compatible_name = f"{phenotype}_allele_compatible"
                    item[compatible_name] = int(item.get(compatible_name, 0)) + 1

    print(
        f"{phenotype}: extracted {extracted} rows; "
        f"allele-compatible rows {compatible}; output {output}"
    )
    return extracted, compatible


def main() -> None:
    markers = read_markers()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for phenotype in PHENOTYPES:
        extract_phenotype(phenotype, markers)


if __name__ == "__main__":
    main()
