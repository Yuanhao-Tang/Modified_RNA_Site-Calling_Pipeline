#!/usr/bin/env python3
"""
Standalone RSEM expression quantification workflow starting from FASTQ input.

Assumptions:
  * config['samples'][sample]['fastq'] provides one or two FASTQ paths
  * config['rsem'] contains at least:
        - work_dir : output root directory
        - index    : RSEM reference prefix
        - params   : optional extra RSEM arguments
        - threads  : optional override for the global thread count

This module builds an `rsem-calculate-expression` command from the configuration
and writes results to `work_dir/<sample>/rsem/`.
"""

from pathlib import Path
import sys
import subprocess
from typing import Iterable, List


def run_rsem_pipeline(sample_name: str, config: dict, dry_run: bool) -> None:
    """Run the RSEM workflow starting from FASTQ files."""
    rsem_config = config.get('rsem', {})
    if not rsem_config:
        raise ValueError("Missing 'rsem' section in the configuration")

    work_dir = Path(
        rsem_config.get(
            'work_dir',
            config.get('global', {}).get('work_dir', './rsem_workdir')
        )
    )
    sample_dir = work_dir / sample_name
    result_dir = sample_dir / 'rsem'
    log_dir = sample_dir / 'logs'

    if not dry_run:
        for d in [result_dir, log_dir]:
            d.mkdir(parents=True, exist_ok=True)

    fastqs = resolve_fastqs(config['samples'][sample_name])
    threads = int(
        rsem_config.get(
            'threads',
            config.get('global', {}).get('threads', 4)
        )
    )
    params = rsem_config.get('params', '')
    reference_prefix = rsem_config.get('index')
    if not reference_prefix:
        raise ValueError("Configuration key 'rsem.index' is empty")

    output_prefix = result_dir / sample_name
    command = build_rsem_command(
        fastqs=fastqs,
        reference_prefix=reference_prefix,
        output_prefix=output_prefix,
        threads=threads,
        extra_params=params,
    )

    log_file = log_dir / f"rsem_{sample_name}.log"
    run_cmd(f"{command} > {log_file} 2>&1", dry_run=dry_run)

    print(f"[OK] RSEM workflow finished for sample: {sample_name}")


def resolve_fastqs(sample_config: dict) -> List[str]:
    """Return FASTQ paths from sample configuration for single-end or paired-end input."""
    fastq_entry = sample_config.get('fastq')
    if fastq_entry is None:
        raise ValueError("Sample configuration is missing the 'fastq' field")

    if isinstance(fastq_entry, str):
        fastqs = [fastq_entry]
    elif isinstance(fastq_entry, Iterable):
        fastqs = list(fastq_entry)
    else:
        raise TypeError("The 'fastq' field must be a string or an iterable")

    if not fastqs:
        raise ValueError("The FASTQ list is empty")

    return fastqs


def build_rsem_command(
    fastqs: List[str],
    reference_prefix: str,
    output_prefix: Path,
    threads: int,
    extra_params: str = "",
) -> str:
    """Build the `rsem-calculate-expression` command line."""
    parts = [
        "rsem-calculate-expression",
        f"--num-threads {threads}",
    ]

    fastqs_str = ""
    if len(fastqs) == 1:
        if fastqs[0].endswith('.gz'):
            parts.append('--gzipped-read-file')
        fastqs_str = fastqs[0]
    elif len(fastqs) == 2:
        parts.append('--paired-end')
        if all(f.endswith('.gz') for f in fastqs):
            parts.append('--gzipped-read-file')
        fastqs_str = ' '.join(fastqs)
    else:
        raise ValueError("Only one or two FASTQ inputs are currently supported")

    if extra_params:
        parts.append(extra_params.strip())

    parts.extend([fastqs_str, reference_prefix, str(output_prefix)])
    return ' '.join(parts)


def run_cmd(cmd: str, dry_run: bool) -> None:
    """Run a command and support dry-run mode."""
    if dry_run:
        print(f"[DRY-RUN] {cmd}")
        return

    result = subprocess.run(cmd, shell=True, executable='/bin/bash')
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}", file=sys.stderr)
        sys.exit(1)


__all__ = [
    'run_rsem_pipeline',
    'resolve_fastqs',
    'build_rsem_command',
    'run_cmd',
]
