#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys


def get_project_root():
    """Return the absolute project root path based on the script location."""
    return Path(__file__).parent.parent.absolute()


def run_star_pipeline(sample_name, config, dry_run):
    """
    Main STAR-based workflow including multi-stage alignment, UMI deduplication,
    modification-site calling, and motif annotation.
    """
    # 0. Initialize directories
    work_dir = Path(config['global']['work_dir'])
    sample_dir = work_dir / sample_name
    align_dir = sample_dir / 'alignment'
    mod_dir = sample_dir / 'modification'
    log_dir = sample_dir / 'logs'
    if not dry_run:
        for d in [align_dir, mod_dir, log_dir]:
            d.mkdir(parents=True, exist_ok=True)

    print('### Step 1: Multi-stage alignment ###')
    stage_results = {}
    current_fastq = config['samples'][sample_name]['fastq']

    for stage in config['alignment']['stages']:
        stage_name = stage['name']
        aligner = stage['aligner']
        stage_dir = align_dir / stage_name
        if not dry_run:
            stage_dir.mkdir(parents=True, exist_ok=True)

        print(f"  {stage_name} alignment ({aligner})")

        if aligner == 'bowtie2':
            result = run_bowtie2(current_fastq, stage, stage_dir, log_dir, dry_run)
        elif aligner == 'star':
            result = run_star(current_fastq, stage, stage_dir, log_dir, config, dry_run)
        else:
            raise ValueError(f"Unknown aligner: {aligner}")

        stage_results[stage_name] = {
            'bam': result['bam'],
            'unmapped': result.get('unmapped')
        }

        # Pass unmapped reads to the next stage.
        if stage.get('keep_unmapped', False):
            unmapped_path = result.get('unmapped')
            if unmapped_path and (Path(unmapped_path).exists() or dry_run):
                current_fastq = result['unmapped']
            else:
                print('  [WARN] Unmapped output was not found; stopping downstream alignment stages')
                break
        else:
            break

    print('### Step 2: Per-stage post-processing ###')
    for stage_name, result in stage_results.items():
        bam = result['bam']
        stage_config = get_stage_config(config, stage_name)

        # 2a UMI deduplication (optional)
        if config.get('umi_dedup', {}).get('enable', False):
            print('  - UMI deduplication')
            bam = umi_dedup(bam, config, log_dir, dry_run)
            result['bam'] = bam

        # 2b Site calling
        print('  - Extract modification sites')
        sites_table = extract_modification_sites(
            bam=bam,
            stage_name=stage_name,
            stage_config=stage_config,
            config=config,
            mod_dir=mod_dir,
            log_dir=log_dir,
            dry_run=dry_run
        )

        # 2c Motif annotation (optional)
        if config.get('motif', {}).get('enable', True):
            print('  - Add motif annotation')
            add_motif(
                sites_table=sites_table,
                stage_config=stage_config,
                config=config,
                dry_run=dry_run
            )

        # 2d BAM -> FASTQ export (optional)
        if config.get('fastq_export', {}).get('enable', False):
            print('  - Export FASTQ')
            export_fastq(
                bam=bam,
                sample_dir=sample_dir,
                stage_name=stage_name,
                config=config,
                log_dir=log_dir,
                dry_run=dry_run
            )

    # 3. Generate alignment statistics report
    if not dry_run:
        print('### Generate alignment statistics report ###')
        try:
            # Add the project root to sys.path so utils can be imported.
            project_root = get_project_root()
            if str(project_root) not in sys.path:
                sys.path.insert(0, str(project_root))

            from utils.alignment_statistics import generate_alignment_report, format_report_text
            report_dir = sample_dir / 'reports'
            report_dir.mkdir(exist_ok=True)
            report_file = report_dir / 'alignment_report.json'
            report_text_file = report_dir / 'alignment_report.txt'

            report = generate_alignment_report(sample_name, sample_dir, config, report_file)
            report_text = format_report_text(report)

            with open(report_text_file, 'w', encoding='utf-8') as f:
                f.write(report_text)

            print(f"  Report written: {report_text_file}")
        except Exception as e:
            print(f"  [WARN] Failed to generate statistics report: {e}")
            import traceback
            if config['global'].get('verbose', False):
                traceback.print_exc()

    # 4. Clean intermediate files (optional)
    if not config['global'].get('keep_intermediate', False):
        print('### Clean intermediate files ###')
        cleanup(sample_dir, dry_run)

    print(f"[OK] Sample {sample_name} completed")


def run_bowtie2(fastq, stage, stage_dir, log_dir, dry_run):
    bam_out = stage_dir / 'aligned.bam'
    keep_unmapped = stage.get('keep_unmapped', False)
    unmapped_fq = stage_dir / 'unmapped.fq.gz' if keep_unmapped else None
    threads = int(stage.get('threads') or 8)
    params = stage.get('params', {}).get('bowtie2', '')

    un_opt = f"--un-gz {unmapped_fq} " if keep_unmapped else ""
    cmd = (
        f"bowtie2 -p {threads} {params} "
        f"-x {stage['index']} "
        f"-U {fastq} "
        f"{un_opt}"
        f"2> {log_dir}/bowtie2_{stage['name']}.log | "
        f"samtools sort -@ {threads} -O BAM -o {bam_out} && "
        f"samtools index {bam_out}"
    )
    run_cmd(cmd, dry_run)
    return {'bam': bam_out, 'unmapped': str(unmapped_fq) if keep_unmapped else None}


def run_star(fastq, stage, stage_dir, log_dir, config, dry_run):
    threads = int(config['global'].get('threads', 8))
    bam_out = stage_dir / 'Aligned.sortedByCoord.out.bam'
    keep_unmapped = stage.get('keep_unmapped', False)
    unmapped_fq = stage_dir / 'Unmapped.out.mate1' if keep_unmapped else None
    unmapped_fq_gz = stage_dir / 'Unmapped.out.mate1.gz' if keep_unmapped else None
    params = stage.get('params', {}).get('star', '')

    unmapped_opt = '--outReadsUnmapped Fastx' if keep_unmapped else ''
    cmd = (
        f"STAR --runThreadN {threads} "
        f"--genomeDir {stage['index']} "
        f"--readFilesIn {fastq} "
        f"--readFilesCommand zcat "
        f"{params} "
        f"--outSAMtype BAM SortedByCoordinate "
        f"--outFileNamePrefix {stage_dir}/ "
        f"{unmapped_opt} "
        f"> {log_dir}/star_{stage['name']}.log 2>&1"
    )
    run_cmd(cmd, dry_run)
    run_cmd(f"samtools index {bam_out}", dry_run)
    # Compress unmapped FASTQ output.
    if keep_unmapped and (unmapped_fq.exists() or dry_run):
        run_cmd(f"gzip {unmapped_fq}", dry_run)
    return {'bam': bam_out, 'unmapped': str(unmapped_fq_gz) if keep_unmapped else None}


def umi_dedup(bam, config, log_dir, dry_run):
    dedup_config = config['umi_dedup']
    tool = dedup_config.get('tool', 'umi_tools')
    dedup_bam = bam.parent / f"{bam.stem}.dedup.bam"

    if tool == 'umi_tools':
        cmd = (
            f"umi_tools dedup "
            f"--method {dedup_config['method']} "
            f"--edit-distance-threshold {dedup_config['edit_distance']} "
            f"-I {bam} -S {dedup_bam} "
            f"> {log_dir}/umi_dedup_{bam.parent.name}.log 2>&1 && "
            f"samtools index {dedup_bam}"
        )
    elif tool == 'umiCollapse':
        cmd = (
            f"umiCollapse -i {bam} -o {dedup_bam} "
            f"--edit-distance {dedup_config['edit_distance']} "
            f"> {log_dir}/umi_dedup_{bam.parent.name}.log 2>&1 && "
            f"samtools index {dedup_bam}"
        )
    else:
        raise ValueError(f"Unknown UMI deduplication tool: {tool}")

    run_cmd(cmd, dry_run)
    return dedup_bam


def extract_modification_sites(bam, stage_name, stage_config, config, mod_dir, log_dir, dry_run):
    sites_table = mod_dir / f"{stage_name}.sites.tsv"
    mod_type = config['modification']['type']
    min_depth = config['modification']['min_depth']
    min_freq = config['modification'].get('min_modification_freq', 0.05)

    # Use an absolute project-root-relative script path.
    project_root = get_project_root()
    script_path = project_root / 'utils' / 'extract_modification_sites.py'

    cmd = (
        f"python {script_path} "
        f"--bam {bam} "
        f"--fasta {stage_config['fasta']} "
        f"--mod-type {mod_type} "
        f"--min-depth {min_depth} "
        f"--min-freq {min_freq} "
        f"--output {sites_table} "
        f"> {log_dir}/extract_sites_{stage_name}.log 2>&1"
    )
    run_cmd(cmd, dry_run)
    return sites_table


def add_motif(sites_table, stage_config, config, dry_run):
    fasta = stage_config['fasta']
    upstream = config['motif']['upstream']
    downstream = config['motif']['downstream']

    # Use an absolute project-root-relative script path.
    project_root = get_project_root()
    script_path = project_root / 'utils' / 'extract_motif.py'

    cmd = (
        f"python {script_path} "
        f"--input {sites_table} "
        f"--fasta {fasta} "
        f"--upstream {upstream} "
        f"--downstream {downstream} "
        f"--inplace"
    )
    run_cmd(cmd, dry_run)


def export_fastq(bam, sample_dir, stage_name, config, log_dir, dry_run):
    fastq_cfg = config.get('fastq_export', {})
    output_subdir = fastq_cfg.get('output_dir', 'fastq')
    fastq_dir = sample_dir / output_subdir
    if not dry_run:
        fastq_dir.mkdir(parents=True, exist_ok=True)

    threads = int(config.get('global', {}).get('threads', 8))
    output_path = fastq_dir / f"{stage_name}.fastq.gz"
    log_file = log_dir / f"fastq_{stage_name}.log"

    cmd = (
        f"samtools fastq "
        f"-@ {threads} "
        f"{bam} | gzip > {output_path} 2> {log_file}"
    )
    run_cmd(cmd, dry_run)


def get_stage_config(config, stage_name):
    for stage in config['alignment']['stages']:
        if stage['name'] == stage_name:
            return stage
    raise ValueError(f"Stage not found: {stage_name}")


def cleanup(sample_dir, dry_run):
    align_dir = sample_dir / 'alignment'
    if not align_dir.exists():
        if dry_run:
            print('  [DRY-RUN] Clean intermediate files (directory does not exist, skipping)')
        return

    # Remove unmapped FASTQ files.
    for unmapped in align_dir.rglob('unmapped*.fq.gz'):
        if dry_run:
            print(f"  [DRY-RUN] Remove: {unmapped}")
        else:
            unmapped.unlink(missing_ok=True)
            print(f"  Remove: {unmapped.name}")

    # Remove pre-dedup BAM files if deduplicated BAM files exist.
    for bam in align_dir.rglob('aligned.bam'):
        dedup_bam = bam.parent / f"{bam.stem}.dedup.bam"
        if dedup_bam.exists() or dry_run:
            if dry_run:
                print(f"  [DRY-RUN] Remove: {bam.name}")
            else:
                bam.unlink(missing_ok=True)
                (bam.parent / f"{bam.name}.bai").unlink(missing_ok=True)
                print(f"  Remove: {bam.name} (deduplicated)")


def run_cmd(cmd, dry_run):
    if dry_run:
        print(f"[DRY-RUN] {cmd}")
        return
    result = subprocess.run(cmd, shell=True, executable='/bin/bash')
    if result.returncode != 0:
        print(f"[ERROR] Command failed: {cmd}", file=sys.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd)


__all__ = [
    'run_star_pipeline', 'run_bowtie2', 'run_star', 'umi_dedup',
    'extract_modification_sites', 'add_motif', 'export_fastq',
    'get_stage_config', 'cleanup', 'run_cmd'
]
