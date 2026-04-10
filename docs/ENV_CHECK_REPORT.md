# Environment Check Report

## Recommended Environment

The recommended environment name is `rnamod-pipeline`, created from `environment.yml`.

## Expected Core Dependencies

### Python Packages

- `pyyaml`
- `pandas`
- `pysam`

### External Tools

- `bowtie2`
- `STAR`
- `samtools`
- `umi_tools`
- `RSEM`

## Suggested Setup

```bash
cd /home/tangyh/script/Modified_RNA_Site-Calling_Pipeline
conda env create -f environment.yml
conda activate rnamod-pipeline
```

If `environment.yml` changes later:

```bash
conda env update -f environment.yml --prune
```

## Practical Notes

- Build all required genome indexes before running the pipeline.
- Build an RSEM index if you plan to use the standalone RSEM workflow.
- If UMI deduplication is enabled, preprocess FASTQ files with `umi_tools extract` before starting the pipeline.

## Conclusion

A complete environment based on `environment.yml` should be sufficient for the main alignment workflow and the standalone RSEM workflow.
