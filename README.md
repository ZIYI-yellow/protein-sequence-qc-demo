# Protein Sequence QC Demo

This private learning repository contains a small, reproducible demonstration of protein-sequence data quality control in Python. It uses 20 reviewed UniProtKB records: 10 known conotoxins and 10 non-toxin neuropeptides.

## What the workflow checks

- FASTA and metadata identifiers match;
- records and metadata are not missing;
- sequences contain valid amino-acid symbols;
- exact duplicate sequences are flagged;
- simple, interpretable sequence features are calculated;
- input provenance and file checksums are recorded.

## Run the clean demonstration

On Windows with Python 3 installed, open PowerShell in this folder and run:

```powershell
.\run_demo.ps1 -OpenReport
```

The report is written to `outputs/protein_sequence_qc_report.html`.

## Run the controlled-error demonstration

```powershell
.\run_error_demo.ps1 -OpenReport
```

This creates disposable copies containing four documented errors and checks whether the same workflow detects them. The source files in `data/` are not changed.

## Scientific boundary

This is a data-readiness and quality-control benchmark, not a toxin classifier. It does not claim discovery of a novel venom peptide or therapeutic lead. The small reference groups are not fully matched by species or protein family, and the records mix precursor proteins with fragment or mature-peptide entries. They therefore cannot support biological efficacy claims or model-performance claims.

## Current status

The repository is private while I reproduce the workflow, review each output, and document what I can explain independently. Public release will follow only after provenance, documentation, and scientific boundaries have been checked.
