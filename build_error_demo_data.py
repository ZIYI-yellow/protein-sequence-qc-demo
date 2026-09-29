#!/usr/bin/env python3
"""Build a controlled, disposable error dataset for the live QC demonstration.

The original UniProtKB files are never changed. This script copies them into a
separate folder and introduces four documented errors so the QC workflow can
show that it detects problems rather than silently accepting them.
"""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parent
SOURCE_DIR = PROJECT_DIR / "data"
ERROR_DIR = PROJECT_DIR / "error_demo_data"


def read_fasta_blocks(path: Path) -> list[tuple[str, str]]:
    blocks: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8-sig") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    blocks.append((header, "".join(chunks)))
                header = line
                chunks = []
            else:
                if header is None:
                    raise ValueError("Sequence appeared before the first FASTA header")
                chunks.append(line)
    if header is not None:
        blocks.append((header, "".join(chunks)))
    return blocks


def write_fasta_blocks(path: Path, blocks: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for header, sequence in blocks:
        lines.append(header)
        lines.extend(sequence[index : index + 60] for index in range(0, len(sequence), 60))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_corrupted_fasta(source: Path, destination: Path) -> tuple[str, str]:
    blocks = read_fasta_blocks(source)
    if len(blocks) < 5:
        raise ValueError("At least five FASTA records are required for the error demo")

    # Controlled error 1: an invalid amino-acid character in record 3.
    invalid_header, invalid_sequence = blocks[2]
    blocks[2] = (invalid_header, "Z" + invalid_sequence[1:])

    # Controlled error 2: record 4 is replaced by record 5's sequence, creating
    # an exact duplicate while preserving two different accessions.
    duplicate_header, _ = blocks[3]
    duplicate_source_header, duplicate_sequence = blocks[4]
    blocks[3] = (duplicate_header, duplicate_sequence)

    write_fasta_blocks(destination, blocks)
    return invalid_header[1:].split()[0], (
        f"{duplicate_header[1:].split()[0]} duplicates "
        f"{duplicate_source_header[1:].split()[0]}"
    )


def make_corrupted_metadata(source: Path, destination: Path) -> tuple[str, str]:
    with source.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)
    if not fieldnames or len(rows) < 2:
        raise ValueError("At least two TSV records are required for the error demo")

    # Controlled error 3: the first metadata accession no longer matches FASTA.
    original_entry = rows[0]["Entry"]
    rows[0]["Entry"] = "BROKEN_ID_001"

    # Controlled error 4: the second metadata length is deliberately wrong.
    length_entry = rows[1]["Entry"]
    rows[1]["Length"] = str(int(rows[1]["Length"]) + 7)

    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return f"{original_entry} changed to BROKEN_ID_001", length_entry


def main() -> int:
    ERROR_DIR.mkdir(parents=True, exist_ok=True)

    positive_fasta = SOURCE_DIR / "known_conotoxins.fasta"
    positive_metadata = SOURCE_DIR / "known_conotoxins_metadata.tsv"
    negative_fasta = SOURCE_DIR / "non_toxin_neuropeptides.fasta"
    negative_metadata = SOURCE_DIR / "non_toxin_neuropeptides_metadata.tsv"

    invalid_record, duplicate_detail = make_corrupted_fasta(
        positive_fasta, ERROR_DIR / positive_fasta.name
    )
    id_detail, length_record = make_corrupted_metadata(
        positive_metadata, ERROR_DIR / positive_metadata.name
    )

    shutil.copy2(negative_fasta, ERROR_DIR / negative_fasta.name)
    shutil.copy2(negative_metadata, ERROR_DIR / negative_metadata.name)

    manifest = [
        "CONTROLLED ERROR INJECTION MANIFEST",
        "The original data files were not changed.",
        "",
        f"1. ID mismatch: {id_detail}",
        f"2. Metadata length increased by 7 for: {length_record}",
        f"3. Invalid amino-acid symbol Z inserted in: {invalid_record}",
        f"4. Exact duplicate sequence created: {duplicate_detail}",
        "",
        "Expected result: the QC workflow should flag these errors.",
    ]
    (ERROR_DIR / "ERRORS_INTRODUCED.txt").write_text(
        "\n".join(manifest) + "\n", encoding="utf-8"
    )
    print("Controlled error data created in:", ERROR_DIR)
    for line in manifest[3:7]:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
