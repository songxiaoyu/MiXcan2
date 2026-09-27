#!/usr/bin/env python3
"""Extract current model SNPs from full hg19 BCAC files and assign hg38 loci.

The prediction weights and LD references use hg38. The complete BCAC files use
hg19 positions in var_name. Only the loci used by the current weights are
lifted, which avoids creating another genome-wide copy of the GWAS files.
"""

from __future__ import annotations

import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

try:
    from pyliftover import LiftOver
except ImportError as exc:
    raise SystemExit(
        "pyliftover is required. Install it with: python3 -m pip install pyliftover"
    ) from exc


SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
WEIGHTS_FILE = WORKFLOW_DIR / "work" / "weights_normalized.csv"
CHAIN_FILE = WORKFLOW_DIR / "data" / "liftover" / "hg38ToHg19.over.chain.gz"
FULL_DIR = WORKFLOW_DIR / "data" / "full_gwas"
OVERALL_SOURCE = FULL_DIR / "bcac_2020_sumstats_for_s-mixcan.txt.gz"
SUBTYPE_SOURCE = FULL_DIR / "bcac_2020_subtype_zscores_for_s-mixcan.txt.gz"
OVERALL_OUTPUT = FULL_DIR / "bcac_2020_sumstats_for_s-mixcan_hg38_weight_subset.csv"
SUBTYPE_OUTPUT = FULL_DIR / "bcac_2020_subtype_zscores_for_s-mixcan_hg38_weight_subset.txt"

KNOWN_OVERALL = (
    WORKFLOW_DIR / "data" / "bcac_2020_sumstats_for_s-mixcan_hg38_repro_subset.csv"
)
KNOWN_SUBTYPE = (
    WORKFLOW_DIR / "data" / "bcac_2020_subtype_zscores_for_s-mixcan_hg38_repro_subset.txt"
)


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise SystemExit("Missing required input(s):\n  " + "\n  ".join(missing))


def normalize_chr(value: str) -> str:
    value = value.strip()
    return value[3:] if value.lower().startswith("chr") else value


def read_weight_markers() -> dict[str, dict[str, object]]:
    markers: dict[str, dict[str, object]] = {}
    with WEIGHTS_FILE.open(newline="") as handle:
        for row in csv.DictReader(handle):
            marker = row["MarkerName_hg38"]
            chrom, pos, ref, alt = marker.split(":", 3)
            markers.setdefault(
                marker,
                {
                    "marker": marker,
                    "chr38": normalize_chr(chrom),
                    "pos38": int(pos),
                    "ref38": ref.upper(),
                    "alt38": alt.upper(),
                },
            )
    return markers


def read_known_mapping() -> dict[tuple[str, int], int]:
    mapping: dict[tuple[str, int], int] = {}
    if KNOWN_OVERALL.is_file():
        with KNOWN_OVERALL.open(newline="") as handle:
            for row in csv.DictReader(handle):
                mapping[(normalize_chr(row["CHR"]), int(row["POS_hg38"]))] = int(
                    row["POS_hg19"]
                )
    if KNOWN_SUBTYPE.is_file():
        with KNOWN_SUBTYPE.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                mapping[(normalize_chr(row["chr"]), int(row["POS_hg38"]))] = int(
                    row["position"]
                )
    return mapping


def liftover_markers(
    markers: dict[str, dict[str, object]], known: dict[tuple[str, int], int]
) -> None:
    converter = LiftOver(str(CHAIN_FILE))
    validation_total = 0
    validation_exact = 0

    for item in markers.values():
        chrom = str(item["chr38"])
        pos38 = int(item["pos38"])
        hits = converter.convert_coordinate(f"chr{chrom}", pos38 - 1) or []
        same_chr_hits = [
            hit for hit in hits if normalize_chr(str(hit[0])) == chrom
        ]
        selected_hits = same_chr_hits or hits
        selected = max(selected_hits, key=lambda hit: hit[3]) if selected_hits else None

        known_pos = known.get((chrom, pos38))
        if selected is not None:
            pos19 = int(selected[1]) + 1
            item["chr19"] = normalize_chr(str(selected[0]))
            item["pos19"] = pos19
            item["strand"] = str(selected[2])
            item["mapping_source"] = "UCSC_chain"
            if known_pos is not None:
                validation_total += 1
                if pos19 != known_pos:
                    raise SystemExit(
                        f"Liftover validation failed for {item['marker']}: "
                        f"chain={pos19}, existing={known_pos}"
                    )
                validation_exact += 1
        elif known_pos is not None:
            item["chr19"] = chrom
            item["pos19"] = known_pos
            item["strand"] = "unknown"
            item["mapping_source"] = "existing_verified_subset"
        else:
            item["chr19"] = ""
            item["pos19"] = ""
            item["strand"] = ""
            item["mapping_source"] = "unmapped"

    print(
        f"Validated {validation_exact}/{validation_total} chain mappings "
        "against existing hg38 subset coordinates."
    )


def complement(allele: str) -> str | None:
    table = str.maketrans("ACGT", "TGCA")
    allele = allele.upper()
    if len(allele) != 1 or allele not in "ACGT":
        return None
    return allele.translate(table)


def alleles_compatible(item: dict[str, object], baseline: str, effect: str) -> bool:
    ref = str(item["ref38"])
    alt = str(item["alt38"])
    baseline = baseline.upper()
    effect = effect.upper()
    if (baseline, effect) in ((ref, alt), (alt, ref)):
        return True
    comp_baseline = complement(baseline)
    comp_effect = complement(effect)
    return (
        comp_baseline is not None
        and comp_effect is not None
        and (comp_baseline, comp_effect) in ((ref, alt), (alt, ref))
    )


def extract_gwas(
    source: Path,
    output: Path,
    delimiter: str,
    markers: dict[str, dict[str, object]],
    count_field: str,
) -> tuple[int, int]:
    targets: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for item in markers.values():
        if item["mapping_source"] != "unmapped":
            targets[(str(item["chr19"]), int(item["pos19"]))].append(item)

    row_count = 0
    compatible_count = 0
    with gzip.open(source, "rt", newline="") as source_handle, output.open(
        "w", newline=""
    ) as output_handle:
        reader = csv.DictReader(source_handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"No header found in {source}")
        output_fields = list(reader.fieldnames) + [
            "CHR",
            "POS_hg19",
            "POS_hg38",
            "MarkerName_hg38",
            "liftover_strand",
            "liftover_source",
        ]
        writer = csv.DictWriter(output_handle, fieldnames=output_fields, delimiter=delimiter)
        writer.writeheader()

        for row in reader:
            parts = row["var_name"].split("_")
            if len(parts) < 4:
                continue
            try:
                key = (normalize_chr(parts[0]), int(parts[1]))
            except ValueError:
                continue
            candidates = targets.get(key)
            if not candidates:
                continue

            for item in candidates:
                current = dict(row)
                current.update(
                    {
                        "CHR": f"chr{item['chr38']}",
                        "POS_hg19": item["pos19"],
                        "POS_hg38": item["pos38"],
                        "MarkerName_hg38": item["marker"],
                        "liftover_strand": item["strand"],
                        "liftover_source": item["mapping_source"],
                    }
                )
                writer.writerow(current)
                row_count += 1
                item[count_field] = int(item.get(count_field, 0)) + 1
                if alleles_compatible(
                    item, row["Baseline.Meta"], row["Effect.Meta"]
                ):
                    compatible_count += 1
                    compatible_field = f"{count_field}_allele_compatible"
                    item[compatible_field] = int(item.get(compatible_field, 0)) + 1

    return row_count, compatible_count


def main() -> None:
    require_files([WEIGHTS_FILE, CHAIN_FILE, OVERALL_SOURCE, SUBTYPE_SOURCE])
    markers = read_weight_markers()
    known = read_known_mapping()
    liftover_markers(markers, known)

    overall_rows, overall_compatible = extract_gwas(
        OVERALL_SOURCE, OVERALL_OUTPUT, ",", markers, "overall_rows"
    )
    subtype_rows, subtype_compatible = extract_gwas(
        SUBTYPE_SOURCE, SUBTYPE_OUTPUT, "\t", markers, "subtype_rows"
    )

    mapped = sum(item["mapping_source"] != "unmapped" for item in markers.values())
    print(f"Mapped weight variants: {mapped}/{len(markers)}")
    print(
        f"Overall BCAC extracted rows: {overall_rows}; "
        f"allele-compatible rows: {overall_compatible}"
    )
    print(
        f"Subtype BCAC extracted rows: {subtype_rows}; "
        f"allele-compatible rows: {subtype_compatible}"
    )
    print(f"Wrote {OVERALL_OUTPUT}")
    print(f"Wrote {SUBTYPE_OUTPUT}")


if __name__ == "__main__":
    main()
