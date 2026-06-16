#!/usr/bin/env python3
"""
Configuration validation utilities.
"""
import sys
from pathlib import Path


def validate_config(config, allow_missing_files=False):
    """Validate the configuration structure and required file paths.

    When `allow_missing_files=True`, file-existence checks are relaxed for dry-run mode.
    """
    required_sections = ['global', 'samples', 'alignment', 'modification']
    for section in required_sections:
        if section not in config:
            print(f"[ERROR] Missing required configuration section: '{section}'", file=sys.stderr)
            sys.exit(1)

    work_dir = Path(config['global']['work_dir'])
    if not work_dir.exists():
        if allow_missing_files:
            print(f"[INFO] Working directory does not exist: {work_dir} (dry-run, skipping creation)")
        else:
            print(f"[INFO] Working directory does not exist and will be created: {work_dir}")
            work_dir.mkdir(parents=True, exist_ok=True)

    for sample_name, sample_info in config['samples'].items():
        fastq_paths = []
        if 'fastq' in sample_info:
            fastq_value = sample_info['fastq']
            if isinstance(fastq_value, (list, tuple)):
                fastq_paths.extend([Path(p) for p in fastq_value])
            else:
                fastq_paths.append(Path(fastq_value))
        else:
            fq1 = sample_info.get('fastq_1')
            fq2 = sample_info.get('fastq_2')
            if fq1:
                fastq_paths.append(Path(fq1))
            if fq2:
                fastq_paths.append(Path(fq2))
        if not fastq_paths:
            print(f"[ERROR] Sample '{sample_name}' is missing FASTQ configuration ('fastq' or fastq_1/fastq_2)", file=sys.stderr)
            sys.exit(1)
        if not allow_missing_files:
            for fastq_path in fastq_paths:
                if not fastq_path.exists():
                    print(f"[ERROR] FASTQ not found: {fastq_path}", file=sys.stderr)
                    sys.exit(1)

    if 'stages' not in config['alignment']:
        print('[ERROR] alignment.stages is not defined', file=sys.stderr)
        sys.exit(1)

    for stage in config['alignment']['stages']:
        stage_name = stage.get('name', '<unknown>')
        if 'index' not in stage:
            print(f"[ERROR] Stage '{stage_name}' is missing 'index'", file=sys.stderr)
            sys.exit(1)
        if 'fasta' not in stage:
            print(f"[ERROR] Stage '{stage_name}' is missing 'fasta'", file=sys.stderr)
            sys.exit(1)
        if not allow_missing_files:
            if not Path(stage['fasta']).exists():
                print(f"[ERROR] Reference FASTA not found: {stage['fasta']}", file=sys.stderr)
                sys.exit(1)

    mod_type = config['modification']['type']
    if mod_type not in ['puMseq', 'Gloriseq', 'Limeseq']:
        print(f"[ERROR] Unknown modification type: {mod_type}", file=sys.stderr)
        sys.exit(1)

    if config.get('umi_dedup', {}).get('enable', False):
        print('[INFO] UMI deduplication is enabled; make sure FASTQ files were preprocessed with umi_tools extract')

    print('[OK] Configuration validation passed')


__all__ = ['validate_config']
