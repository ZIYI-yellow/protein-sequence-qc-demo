#!/usr/bin/env python3
"""Small, reproducible protein-sequence QC benchmark.

This script is deliberately a quality-control and descriptive comparison demo.
It does not classify unknown sequences, discover toxins, or make drug claims.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import html
import json
import math
import platform
import statistics
import sys
from pathlib import Path


VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
FEATURE_GROUPS = {
    "basic_pct": set("KRH"),
    "acidic_pct": set("DE"),
    "hydrophobic_pct": set("AVILMFWY"),
    "aromatic_pct": set("FWY"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a small, auditable protein-sequence QC benchmark."
    )
    parser.add_argument("--config", default="config.json", help="Path to JSON config")
    parser.add_argument("--output", default="outputs", help="Output directory")
    return parser.parse_args()


def read_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    required = {"project_title", "purpose", "boundary", "datasets"}
    missing = sorted(required - set(config))
    if missing:
        raise ValueError(f"Config is missing: {', '.join(missing)}")
    if not config["datasets"]:
        raise ValueError("Config must contain at least one dataset")
    return config


def clean_sequence(text: str) -> str:
    # Remove layout whitespace only. Other unexpected characters are retained
    # so the QC step can report them instead of silently deleting them.
    return "".join(text.upper().split())


def fasta_accession(header: str) -> str:
    """Return accession from UniProt-style header, with a safe fallback."""
    first = header.split()[0]
    parts = first.split("|")
    return parts[1] if len(parts) >= 3 else first


def read_fasta(path: Path) -> dict[str, dict]:
    records: dict[str, dict] = {}
    header: str | None = None
    chunks: list[str] = []

    def save_record() -> None:
        if header is None:
            return
        accession = fasta_accession(header)
        if accession in records:
            raise ValueError(f"Duplicate FASTA accession in {path.name}: {accession}")
        records[accession] = {
            "fasta_header": header,
            "fasta_sequence": clean_sequence("".join(chunks)),
        }

    with path.open("r", encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save_record()
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence before first FASTA header in {path.name}, line {line_number}"
                    )
                chunks.append(line)
    save_record()
    return records


def read_tsv(path: Path) -> dict[str, dict]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Entry",
            "Entry Name",
            "Protein names",
            "Organism",
            "Length",
            "Sequence",
            "Subcellular location [CC]",
            "Keywords",
        }
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path.name} is missing columns: {', '.join(missing)}")
        rows: dict[str, dict] = {}
        for row in reader:
            accession = (row.get("Entry") or "").strip()
            if not accession:
                raise ValueError(f"Empty Entry value in {path.name}")
            if accession in rows:
                raise ValueError(f"Duplicate TSV accession in {path.name}: {accession}")
            rows[accession] = row
    return rows


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def amino_acid_percentage(sequence: str, selected: set[str]) -> float:
    if not sequence:
        return math.nan
    return 100.0 * sum(char in selected for char in sequence) / len(sequence)


def kmers(sequence: str, k: int) -> set[str]:
    if len(sequence) < k:
        return set()
    return {sequence[index : index + k] for index in range(len(sequence) - k + 1)}


def jaccard_similarity(first: str, second: str, k: int) -> float:
    first_set = kmers(first, k)
    second_set = kmers(second, k)
    union = first_set | second_set
    return len(first_set & second_set) / len(union) if union else math.nan


def analyse_dataset(dataset: dict, project_dir: Path) -> tuple[list[dict], dict, list[dict]]:
    label = dataset["label"]
    fasta_path = (project_dir / dataset["fasta"]).resolve()
    metadata_path = (project_dir / dataset["metadata"]).resolve()
    fasta = read_fasta(fasta_path)
    metadata = read_tsv(metadata_path)

    fasta_ids = set(fasta)
    metadata_ids = set(metadata)
    missing_in_metadata = sorted(fasta_ids - metadata_ids)
    missing_in_fasta = sorted(metadata_ids - fasta_ids)
    shared_ids = sorted(fasta_ids & metadata_ids)

    rows: list[dict] = []
    for accession in shared_ids:
        fasta_record = fasta[accession]
        meta = metadata[accession]
        fasta_sequence = fasta_record["fasta_sequence"]
        metadata_sequence = clean_sequence(meta["Sequence"])
        invalid = sorted(set(fasta_sequence) - VALID_AMINO_ACIDS)
        try:
            metadata_length = int(meta["Length"])
        except ValueError as exc:
            raise ValueError(f"Invalid metadata Length for {accession}: {meta['Length']}") from exc

        required_nonempty = ["Entry Name", "Protein names", "Organism", "Length", "Sequence"]
        missing_required_fields = [
            field for field in required_nonempty if not (meta.get(field) or "").strip()
        ]
        row = {
            "entry": accession,
            "entry_name": meta["Entry Name"],
            "protein_name": meta["Protein names"],
            "organism": meta["Organism"],
            "group": label,
            "benchmark_role": dataset["benchmark_role"],
            "fasta_length": len(fasta_sequence),
            "metadata_length": metadata_length,
            "length_match": len(fasta_sequence) == metadata_length,
            "sequence_match": fasta_sequence == metadata_sequence,
            "invalid_residues": "".join(invalid),
            "invalid_residue_count": sum(char not in VALID_AMINO_ACIDS for char in fasta_sequence),
            "missing_required_fields": ";".join(missing_required_fields),
            "missing_required_field_count": len(missing_required_fields),
            "cysteine_count": fasta_sequence.count("C"),
            "cysteine_pct": amino_acid_percentage(fasta_sequence, {"C"}),
            "annotated_secreted": "secreted" in meta["Subcellular location [CC]"].lower(),
            "toxin_keyword_present": "toxin" in {
                item.strip().lower() for item in meta["Keywords"].split(";")
            },
            "fasta_header": fasta_record["fasta_header"],
            "sequence": fasta_sequence,
        }
        for feature_name, selected in FEATURE_GROUPS.items():
            row[feature_name] = amino_acid_percentage(fasta_sequence, selected)
        rows.append(row)

    id_report = []
    for accession in sorted(fasta_ids | metadata_ids):
        id_report.append(
            {
                "group": label,
                "entry": accession,
                "in_fasta": accession in fasta_ids,
                "in_metadata": accession in metadata_ids,
                "matched": accession in shared_ids,
            }
        )

    manifest_rows = [
        {
            "group": label,
            "role": dataset["benchmark_role"],
            "file_type": "FASTA",
            "file_name": fasta_path.name,
            "relative_path": str(fasta_path.relative_to(project_dir)),
            "sha256": file_sha256(fasta_path),
            "file_size_bytes": fasta_path.stat().st_size,
            "record_count": len(fasta),
            "source": dataset["source"],
            "download_date": dataset["download_date"],
        },
        {
            "group": label,
            "role": dataset["benchmark_role"],
            "file_type": "TSV metadata",
            "file_name": metadata_path.name,
            "relative_path": str(metadata_path.relative_to(project_dir)),
            "sha256": file_sha256(metadata_path),
            "file_size_bytes": metadata_path.stat().st_size,
            "record_count": len(metadata),
            "source": dataset["source"],
            "download_date": dataset["download_date"],
        },
    ]

    dataset_qc = {
        "group": label,
        "fasta_records": len(fasta),
        "metadata_records": len(metadata),
        "matched_records": len(shared_ids),
        "missing_in_metadata": missing_in_metadata,
        "missing_in_fasta": missing_in_fasta,
    }
    return rows, dataset_qc, id_report + manifest_rows


def mark_global_duplicates(rows: list[dict]) -> None:
    """Flag identical sequences across all groups, including cross-group leakage."""
    sequence_counts: dict[str, int] = {}
    for row in rows:
        sequence_counts[row["sequence"]] = sequence_counts.get(row["sequence"], 0) + 1
    for row in rows:
        row["exact_duplicate"] = sequence_counts[row["sequence"]] > 1


def summarise_groups(rows: list[dict]) -> list[dict]:
    summaries = []
    for group in sorted({row["group"] for row in rows}):
        selected = [row for row in rows if row["group"] == group]
        lengths = [row["fasta_length"] for row in selected]
        cysteine = [row["cysteine_pct"] for row in selected]
        summaries.append(
            {
                "group": group,
                "n": len(selected),
                "median_length": statistics.median(lengths),
                "min_length": min(lengths),
                "max_length": max(lengths),
                "median_cysteine_pct": statistics.median(cysteine),
                "min_cysteine_pct": min(cysteine),
                "max_cysteine_pct": max(cysteine),
                "secreted_n": sum(row["annotated_secreted"] for row in selected),
                "toxin_keyword_n": sum(row["toxin_keyword_present"] for row in selected),
                "length_match_n": sum(row["length_match"] for row in selected),
                "sequence_match_n": sum(row["sequence_match"] for row in selected),
                "invalid_residue_total": sum(row["invalid_residue_count"] for row in selected),
                "exact_duplicate_n": sum(row["exact_duplicate"] for row in selected),
            }
        )
    return summaries


def pairwise_similarity(rows: list[dict], k: int) -> list[dict]:
    output = []
    for group in sorted({row["group"] for row in rows}):
        selected = sorted(
            (row for row in rows if row["group"] == group), key=lambda row: row["entry"]
        )
        for first_index, first in enumerate(selected):
            for second in selected[first_index + 1 :]:
                output.append(
                    {
                        "group": group,
                        "entry_1": first["entry"],
                        "entry_2": second["entry"],
                        "k": k,
                        "kmer_jaccard": jaccard_similarity(
                            first["sequence"], second["sequence"], k
                        ),
                    }
                )
    return sorted(output, key=lambda row: row["kmer_jaccard"], reverse=True)


def write_csv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    if not rows:
        return
    columns = columns or list(rows[0])
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def format_value(value) -> str:
    if isinstance(value, bool):
        return "PASS" if value else "FAIL"
    if isinstance(value, float):
        return f"{value:.2f}"
    return str(value)


def make_svg(rows: list[dict], summaries: list[dict]) -> str:
    width, height = 920, 360
    colours = {"Known conotoxin": "#0f766e", "Non-toxin neuropeptide": "#d97706"}
    groups = [summary["group"] for summary in summaries]
    panels = [
        ("fasta_length", "Sequence length (amino acids)", 50, 420),
        ("cysteine_pct", "Cysteine content (%)", 500, 870),
    ]
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#1f2937}.axis{stroke:#64748b;stroke-width:1}.grid{stroke:#e2e8f0;stroke-width:1}.median{stroke:#111827;stroke-width:3}</style>',
    ]
    for feature, title, left, right in panels:
        values = [float(row[feature]) for row in rows]
        low = 0.0
        high = max(values) * 1.12 if values else 1.0
        top, bottom = 45, 285
        parts.append(f'<text x="{(left + right) / 2}" y="24" text-anchor="middle" font-size="16" font-weight="700">{html.escape(title)}</text>')
        for tick_index in range(6):
            tick_value = low + (high - low) * tick_index / 5
            y = bottom - (bottom - top) * tick_index / 5
            parts.append(f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}"/>')
            parts.append(f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" font-size="11">{tick_value:.1f}</text>')
        parts.append(f'<line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}"/>')
        parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}"/>')
        for group_index, group in enumerate(groups):
            x = left + (right - left) * (0.32 + 0.36 * group_index)
            selected = [row for row in rows if row["group"] == group]
            for point_index, row in enumerate(selected):
                value = float(row[feature])
                y = bottom - (bottom - top) * (value - low) / (high - low)
                jitter = ((point_index % 5) - 2) * 6 + ((point_index // 5) * 3)
                colour = colours.get(group, "#475569")
                parts.append(f'<circle cx="{x + jitter:.1f}" cy="{y:.1f}" r="5" fill="{colour}" fill-opacity="0.82" stroke="#ffffff" stroke-width="1"/>')
            median = statistics.median(float(row[feature]) for row in selected)
            median_y = bottom - (bottom - top) * (median - low) / (high - low)
            parts.append(f'<line class="median" x1="{x - 30:.1f}" y1="{median_y:.1f}" x2="{x + 30:.1f}" y2="{median_y:.1f}"/>')
            short_label = "Known conotoxin" if group == "Known conotoxin" else "Non-toxin control"
            parts.append(f'<text x="{x:.1f}" y="{bottom + 22}" text-anchor="middle" font-size="12">{html.escape(short_label)}</text>')
    parts.append('<text x="460" y="342" text-anchor="middle" font-size="11" fill="#475569">Each point is one reviewed UniProtKB sequence; black line = median. Descriptive only, not a classifier.</text>')
    parts.append("</svg>")
    return "".join(parts)


def html_table(rows: list[dict], columns: list[str]) -> str:
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body = []
    for row in rows:
        cells = "".join(
            f"<td>{html.escape(format_value(row.get(column, '')))}</td>" for column in columns
        )
        body.append(f"<tr>{cells}</tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def make_html_report(
    config: dict,
    rows: list[dict],
    summaries: list[dict],
    dataset_qc: list[dict],
    similarities: list[dict],
    svg: str,
    run_time: str,
) -> str:
    all_lengths_match = all(row["length_match"] for row in rows)
    all_sequences_match = all(row["sequence_match"] for row in rows)
    invalid_total = sum(row["invalid_residue_count"] for row in rows)
    missing_field_total = sum(row["missing_required_field_count"] for row in rows)
    duplicate_total = sum(row["exact_duplicate"] for row in rows)
    unmatched_total = sum(
        len(item["missing_in_metadata"]) + len(item["missing_in_fasta"])
        for item in dataset_qc
    )
    summary_rows = []
    for summary in summaries:
        summary_rows.append(
            {
                "Group": summary["group"],
                "n": summary["n"],
                "Length median (range)": f"{summary['median_length']:.0f} ({summary['min_length']}-{summary['max_length']})",
                "Cysteine median (range)": f"{summary['median_cysteine_pct']:.2f}% ({summary['min_cysteine_pct']:.2f}-{summary['max_cysteine_pct']:.2f}%)",
                "Secreted": f"{summary['secreted_n']}/{summary['n']}",
                "Toxin keyword": f"{summary['toxin_keyword_n']}/{summary['n']}",
            }
        )
    top_pairs = [
        {
            "Group": row["group"],
            "Entry 1": row["entry_1"],
            "Entry 2": row["entry_2"],
            "3-mer Jaccard": f"{row['kmer_jaccard']:.3f}",
        }
        for row in similarities[:5]
    ]
    checks = [
        ("FASTA IDs matched metadata IDs", unmatched_total == 0),
        ("FASTA lengths matched metadata", all_lengths_match),
        ("FASTA sequences matched TSV sequences", all_sequences_match),
        ("No missing required metadata fields", missing_field_total == 0),
        ("No invalid amino-acid symbols", invalid_total == 0),
        ("No exact duplicate sequences", duplicate_total == 0),
    ]
    check_html = "".join(
        f'<div class="check {"pass" if passed else "fail"}"><span>{"PASS" if passed else "CHECK"}</span>{html.escape(label)}</div>'
        for label, passed in checks
    )
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Protein Sequence QC Demo</title>
<style>
body{{font-family:Arial,sans-serif;margin:0;background:#f4f7fb;color:#172033;line-height:1.45}}
main{{max-width:1040px;margin:24px auto;padding:0 20px 50px}}
.hero{{background:linear-gradient(135deg,#0f4c5c,#0f766e);color:white;padding:26px 30px;border-radius:16px}}
.hero h1{{margin:0 0 8px;font-size:28px}} .hero p{{margin:5px 0;max-width:850px}}
.tag{{display:inline-block;margin-top:10px;padding:5px 9px;background:#ffffff22;border:1px solid #ffffff55;border-radius:999px;font-size:13px}}
section{{background:white;margin-top:16px;padding:22px 24px;border-radius:14px;box-shadow:0 3px 12px #1f29370d}}
h2{{font-size:19px;margin:0 0 14px;color:#0f4c5c}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:10px}}
.check{{padding:10px 12px;background:#f8fafc;border-left:4px solid #64748b;border-radius:6px}}
.check span{{font-size:11px;font-weight:bold;margin-right:8px;padding:3px 5px;border-radius:4px}}
.check.pass{{border-color:#15803d}} .check.pass span{{background:#dcfce7;color:#166534}}
.check.fail{{border-color:#b91c1c}} .check.fail span{{background:#fee2e2;color:#991b1b}}
table{{width:100%;border-collapse:collapse;font-size:13px}} th{{text-align:left;background:#e8f3f3;color:#164e63}}
th,td{{padding:9px 10px;border-bottom:1px solid #e2e8f0;vertical-align:top}}
.note{{background:#fff7ed;border-left:4px solid #ea580c;padding:12px 14px;border-radius:6px}}
.flow{{font-family:Consolas,monospace;background:#0f172a;color:#e2e8f0;padding:14px;border-radius:8px;white-space:pre-wrap}}
.small{{font-size:12px;color:#64748b}} svg{{max-width:100%;height:auto}}
</style>
</head>
<body><main>
<div class="hero">
  <h1>{html.escape(config['project_title'])}</h1>
  <p>{html.escape(config['purpose'])}</p>
  <div class="tag">Python standard library | offline | reproducible | auditable</div>
</div>
<section><h2>Workflow</h2><div class="flow">FASTA + TSV -&gt; ID matching -&gt; sequence QC -&gt; interpretable features -&gt; descriptive benchmark -&gt; traceable outputs</div></section>
<section><h2>Quality-control checks</h2><div class="grid">{check_html}</div></section>
<section><h2>Small benchmark summary</h2>{html_table(summary_rows, ['Group','n','Length median (range)','Cysteine median (range)','Secreted','Toxin keyword'])}</section>
<section><h2>Descriptive feature comparison</h2>{svg}</section>
<section><h2>Redundancy flag</h2>
<p>The 3-mer Jaccard value is used only to flag related records. It is not an alignment or a homology estimate.</p>
{html_table(top_pairs, ['Group','Entry 1','Entry 2','3-mer Jaccard'])}</section>
<section><h2>Scientific boundary</h2><div class="note">{html.escape(config['boundary'])}</div>
<p>This reference set mixes precursor proteins with fragment or mature-peptide records. Their lengths and signal-peptide status are therefore not directly comparable. A project-stage analysis should first standardise the sequence level, then use supervisor-approved species, genomic data type, toxin reference set, signal-peptide and transmembrane tools, homology and domain analysis, and an experimental validation plan.</p></section>
<section><h2>Reproducible outputs</h2><p>Source manifest with SHA-256 checksums, ID-match report, sequence feature table, group summary, pairwise redundancy table, QC report and run log.</p><p class="small">Generated {html.escape(run_time)} | Python {html.escape(platform.python_version())}</p></section>
</main></body></html>"""


def main() -> int:
    args = parse_args()
    config_path = Path(args.config).resolve()
    project_dir = config_path.parent
    output_dir = (project_dir / args.output).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_time = dt.datetime.now().astimezone().isoformat(timespec="seconds")

    print("[1/6] Reading configuration")
    config = read_config(config_path)

    all_rows: list[dict] = []
    all_dataset_qc: list[dict] = []
    all_id_rows: list[dict] = []
    all_manifest_rows: list[dict] = []

    print("[2/6] Reading FASTA and TSV files")
    for dataset in config["datasets"]:
        rows, dataset_qc, mixed_rows = analyse_dataset(dataset, project_dir)
        all_rows.extend(rows)
        all_dataset_qc.append(dataset_qc)
        for row in mixed_rows:
            if "file_type" in row:
                all_manifest_rows.append(row)
            else:
                all_id_rows.append(row)

    print("[3/6] Matching IDs and checking sequence quality")
    mark_global_duplicates(all_rows)
    summaries = summarise_groups(all_rows)

    print("[4/6] Calculating interpretable sequence features")
    k = int(config.get("kmer_size", 3))
    similarities = pairwise_similarity(all_rows, k)

    print("[5/6] Writing tables, manifest, checksums and log")
    feature_columns = [
        "entry", "entry_name", "protein_name", "organism", "group", "benchmark_role",
        "fasta_length", "metadata_length", "length_match", "sequence_match",
        "invalid_residues", "invalid_residue_count", "missing_required_fields",
        "missing_required_field_count", "cysteine_count", "cysteine_pct",
        "basic_pct", "acidic_pct", "hydrophobic_pct", "aromatic_pct",
        "annotated_secreted", "toxin_keyword_present", "exact_duplicate",
    ]
    write_csv(output_dir / "sequence_features.csv", all_rows, feature_columns)
    write_csv(output_dir / "group_summary.csv", summaries)
    write_csv(output_dir / "id_match_report.csv", all_id_rows)
    write_csv(output_dir / "source_manifest.csv", all_manifest_rows)
    write_csv(output_dir / "within_group_3mer_similarity.csv", similarities)

    svg = make_svg(all_rows, summaries)
    (output_dir / "comparison_plot.svg").write_text(svg, encoding="utf-8")

    unmatched = sum(
        len(item["missing_in_metadata"]) + len(item["missing_in_fasta"])
        for item in all_dataset_qc
    )
    qc_lines = [
        "PROTEIN SEQUENCE QC DEMO",
        f"Run time: {run_time}",
        f"Records analysed: {len(all_rows)}",
        f"Unmatched FASTA/TSV IDs: {unmatched}",
        f"Length mismatches: {sum(not row['length_match'] for row in all_rows)}",
        f"FASTA/TSV sequence mismatches: {sum(not row['sequence_match'] for row in all_rows)}",
        f"Missing required metadata fields: {sum(row['missing_required_field_count'] for row in all_rows)}",
        f"Invalid amino-acid symbols: {sum(row['invalid_residue_count'] for row in all_rows)}",
        f"Records involved in exact duplicates: {sum(row['exact_duplicate'] for row in all_rows)}",
        "",
        "BOUNDARY",
        config["boundary"],
    ]
    (output_dir / "qc_report.txt").write_text("\n".join(qc_lines) + "\n", encoding="utf-8")
    log_lines = [
        f"run_time={run_time}",
        f"python_version={platform.python_version()}",
        f"python_executable={Path(sys.executable).name}",
        f"platform={platform.platform()}",
        f"config_sha256={file_sha256(config_path)}",
        f"records={len(all_rows)}",
        f"kmer_size={k}",
    ]
    (output_dir / "run_log.txt").write_text("\n".join(log_lines) + "\n", encoding="utf-8")

    print("[6/6] Building offline HTML report")
    report = make_html_report(
        config, all_rows, summaries, all_dataset_qc, similarities, svg, run_time
    )
    report_path = output_dir / "protein_sequence_qc_report.html"
    report_path.write_text(report, encoding="utf-8")

    print("\nAnalysis complete")
    print(f"Records: {len(all_rows)}")
    print(f"QC report: {output_dir / 'qc_report.txt'}")
    print(f"HTML report: {report_path}")
    print("Boundary: descriptive QC benchmark, not a toxin classifier.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
