# Modified RNA Site-Calling Pipeline

Author: Tang YH

This project provides a lightweight RNA modification analysis pipeline for multi-stage alignment, optional UMI deduplication, modification-site calling, motif extraction, and standalone RSEM quantification.

## Project Layout

```text
Modified_RNA_Site-Calling_Pipeline/
├── rnamod
├── environment.yml
├── config/
├── pipelines/
├── utils/
├── lib/
└── docs/
```

## Main Components

- `rnamod`: main command-line entry point
- `pipelines/star_pipeline.py`: multi-stage alignment and downstream processing
- `pipelines/rsem_pipeline.py`: standalone RSEM workflow starting from FASTQ
- `utils/extract_modification_sites.py`: BaseCounter-based site calling
- `utils/extract_motif.py`: motif extraction around candidate sites
- `utils/alignment_statistics.py`: alignment, UMI, and coverage summary generation
- `lib/config_parser.py`: configuration validation

## Environment Setup

```bash
conda env create -f environment.yml
conda activate rnamod-pipeline
```

## Quick Start

Run all samples:

```bash
./rnamod -c config/example.yaml
```

Run a single sample:

```bash
./rnamod -c config/example.yaml -s sample1
```

Run multiple samples in parallel:

```bash
./rnamod -c config/example.yaml -j 4
```

Dry run only:

```bash
./rnamod -c config/example.yaml --dry-run
```

## Core Features

- Multi-stage alignment with `bowtie2` and `STAR`
- Optional UMI deduplication with `umi_tools` or `umiCollapse`
- Modification-site extraction for `puMseq`, `Gloriseq`, and `Limeseq`
- Motif extraction from the reference FASTA
- Optional BAM-to-FASTQ export
- Sample-level alignment and coverage summaries
- Standalone RSEM quantification workflow

## Notes

- The pipeline is designed to be runnable from any working directory as long as the `rnamod` path is correct.
- Configuration examples are available in `config/example.yaml` and `config/rsem_example.yaml`.
- Detailed usage notes are available in `docs/`.
