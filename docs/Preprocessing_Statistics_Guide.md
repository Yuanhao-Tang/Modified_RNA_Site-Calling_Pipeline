# Preprocessing Statistics Guide

## Overview

This document explains the tables produced by the `SampleStatistics` module. These files help assess preprocessing quality after sample merging and genomic annotation.

## Generated Tables

`SampleStatistics` produces five TSV tables:

1. `sample_summary.tsv`
2. `chromosome_stats.tsv`
3. `region_stats.tsv`
4. `gene_type_stats.tsv`
5. `mod_rate_bins.tsv`

## How To Use These Tables

### `sample_summary.tsv`

Use this table for a quick sample-level overview.

- `valid_sites / total_sites` gives a rough sense of completeness.
- A high missing rate may suggest low depth or overly strict filtering.
- Differences between `mean_mod_rate` and `median_mod_rate` can reveal whether modification signals are concentrated in a small subset of sites.

### `chromosome_stats.tsv`

Use this table to inspect chromosome-level patterns.

- `percentage` indicates how many sites fall on each chromosome.
- Large chromosome-specific shifts in `mean_mod_rate` may reflect biological or technical effects.

### `region_stats.tsv`

Use this table to compare genomic regions such as `exon`, `intron`, and `intergenic`.

- Compare `percentage` across regions to understand where sites are enriched.
- Compare `mean_mod_rate` and `mean_depth` to identify region-specific biases or signal enrichment.

### `gene_type_stats.tsv`

Use this table to compare gene-type categories.

- Common categories include `protein_coding`, `lncRNA`, `miRNA`, and `processed_pseudogene`.
- Differences in `mean_mod_rate` across gene types can guide downstream interpretation.

### `mod_rate_bins.tsv`

Use this table to inspect the distribution of modification rates.

- Low-rate bins often contain most sites.
- A sample with an unusually high fraction of high-rate sites may deserve additional review.

## Reading The Files

```python
import pandas as pd

sample_summary = pd.read_csv("preprocessing_stats_sample_summary.tsv", sep="	")
chrom_stats = pd.read_csv("preprocessing_stats_chromosome_stats.tsv", sep="	")
region_stats = pd.read_csv("preprocessing_stats_region_stats.tsv", sep="	")
gene_type_stats = pd.read_csv("preprocessing_stats_gene_type_stats.tsv", sep="	")
mod_bins = pd.read_csv("preprocessing_stats_mod_rate_bins.tsv", sep="	")
```
