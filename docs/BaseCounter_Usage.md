# BaseCounter Usage

## Overview

`BaseCounter` is a pileup-based base counting tool for RNA modification analysis. It supports `puMseq`, `Gloriseq`, and `Limeseq`, and can be used either through the command line or as a Python module.

## Key Features

- Fixed `orient_to_reference=False`
- Configurable sequencing method handling
- Optional strand-separated or merged output
- Minimum depth and minimum frequency filtering
- Streaming output mode for lower memory usage

## Sequencing Modes

### `puMseq`

- Focuses on reference `T` sites
- Useful for detecting T-to-C signatures associated with pseudouridine-related workflows

### `Gloriseq`

- Focuses on reference `A` sites
- Useful for detecting A-to-G signatures in Glori-seq style workflows

### `Limeseq`

- Reports all four reference bases
- Useful for broader modification screening

## Command-Line Examples

```bash
python utils/extract_modification_sites.py   --bam input.bam   --fasta reference.fa   --mod-type puMseq   --min-depth 10   --min-freq 0.05   --output pumseq_sites.tsv.gz
```

```bash
python utils/extract_modification_sites.py   --bam input.bam   --fasta reference.fa   --mod-type Limeseq   --min-depth 10   --merge-strands   --output limeseq_sites.tsv.gz
```

```bash
python utils/extract_modification_sites.py   --bam input.bam   --fasta reference.fa   --mod-type custom   --custom-ref-bases A,T   --min-depth 15   --output custom_AT.tsv.gz
```

## Python API Example

```python
from utils.extract_modification_sites import BaseCounter, SeqMethod

counter = BaseCounter(
    bam_path="input.bam",
    ref_fasta_path="reference.fa",
    seq_method=SeqMethod.PUMSEQ,
    min_depth=10,
    min_freq=0.05,
)

df = counter.run()
counter.save_to_csv(df, "output.tsv.gz", compression="gzip")
```

## Important Arguments

- `--bam`: input BAM path
- `--fasta`: reference FASTA path
- `--mod-type`: `puMseq`, `Gloriseq`, `Limeseq`, or `custom`
- `--output`: output table path
- `--min-depth`: minimum depth threshold
- `--min-freq`: minimum non-reference frequency threshold
- `--merge-strands`: merge `+` and `-` strands into one row per site
- `--custom-ref-bases`: comma-separated base set for `custom` mode

## Output Columns

Strand-separated mode:

```text
chrom  pos  ref  strand  depth  A  C  G  T  N
```

Merged mode:

```text
chrom  pos  ref  depth  A  C  G  T  N
```

## Integration

The main pipeline calls `BaseCounter` automatically from `pipelines/star_pipeline.py` based on the `modification` section of the YAML configuration.
