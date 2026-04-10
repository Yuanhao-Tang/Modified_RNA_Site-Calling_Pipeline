#!/usr/bin/env python3
"""
Alignment statistics report generation module
Extract statistics from bowtie2, STAR, UMI deduplication logs, and BAM files
"""

import re
import json
from pathlib import Path
from typing import Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

try:
    import pysam
except ImportError:
    pysam = None


def parse_bowtie2_log(log_file: Path) -> Dict:
    """Parse a bowtie2 log file and extract alignment statistics"""
    stats = {
        'total_reads': 0,
        'aligned_once': 0,
        'aligned_more_than_once': 0,
        'unaligned': 0,
        'overall_alignment_rate': 0.0
    }
    
    if not log_file.exists():
        return stats
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Extract total read count
    match = re.search(r'(\d+) reads;', content)
    if match:
        stats['total_reads'] = int(match.group(1))
    
    # Extract uniquely aligned reads
    match = re.search(r'(\d+) \([\d.]+%\) aligned exactly 1 time', content)
    if match:
        stats['aligned_once'] = int(match.group(1))
    
    # Extract multi-mapped reads
    match = re.search(r'(\d+) \([\d.]+%\) aligned >1 times', content)
    if match:
        stats['aligned_more_than_once'] = int(match.group(1))
    
    # Extract unaligned reads
    match = re.search(r'(\d+) \([\d.]+%\) aligned 0 times', content)
    if match:
        stats['unaligned'] = int(match.group(1))
    
    # Extract overall alignment rate
    match = re.search(r'overall alignment rate ([\d.]+)%', content)
    if match:
        stats['overall_alignment_rate'] = float(match.group(1))
    
    # Compute unique alignment rate
    if stats['total_reads'] > 0:
        stats['unique_alignment_rate'] = (stats['aligned_once'] / stats['total_reads']) * 100
        stats['multi_alignment_rate'] = (stats['aligned_more_than_once'] / stats['total_reads']) * 100
    else:
        stats['unique_alignment_rate'] = 0.0
        stats['multi_alignment_rate'] = 0.0
    
    return stats


def parse_star_log(log_file: Path) -> Dict:
    """Parse a STAR log file and extract alignment statistics"""
    stats = {
        'total_reads': 0,
        'uniquely_mapped': 0,
        'uniquely_mapped_pct': 0.0,
        'multi_mapped': 0,
        'multi_mapped_pct': 0.0,
        'too_many_mapped': 0,
        'too_many_mapped_pct': 0.0,
        'unmapped_too_many_mismatches': 0,
        'unmapped_too_short': 0,
        'unmapped_other': 0,
        'mismatch_rate': 0.0,
        'deletion_rate': 0.0,
        'insertion_rate': 0.0
    }
    
    if not log_file.exists():
        return stats
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Extract total read count
    match = re.search(r'Number of input reads\s+\|\s+(\d+)', content)
    if match:
        stats['total_reads'] = int(match.group(1))
    
    # Extract uniquely aligned reads
    match = re.search(r'Uniquely mapped reads number\s+\|\s+(\d+)', content)
    if match:
        stats['uniquely_mapped'] = int(match.group(1))
    
    match = re.search(r'Uniquely mapped reads %\s+\|\s+([\d.]+)%', content)
    if match:
        stats['uniquely_mapped_pct'] = float(match.group(1))
    
    # Extract multi-mapped reads
    match = re.search(r'Number of reads mapped to multiple loci\s+\|\s+(\d+)', content)
    if match:
        stats['multi_mapped'] = int(match.group(1))
    
    match = re.search(r'% of reads mapped to multiple loci\s+\|\s+([\d.]+)%', content)
    if match:
        stats['multi_mapped_pct'] = float(match.group(1))
    
    # Extract reads mapped to too many loci
    match = re.search(r'Number of reads mapped to too many loci\s+\|\s+(\d+)', content)
    if match:
        stats['too_many_mapped'] = int(match.group(1))
    
    match = re.search(r'% of reads mapped to too many loci\s+\|\s+([\d.]+)%', content)
    if match:
        stats['too_many_mapped_pct'] = float(match.group(1))
    
    # Extract reasons for unmapped reads
    match = re.search(r'Number of reads unmapped: too many mismatches\s+\|\s+(\d+)', content)
    if match:
        stats['unmapped_too_many_mismatches'] = int(match.group(1))
    
    match = re.search(r'Number of reads unmapped: too short\s+\|\s+(\d+)', content)
    if match:
        stats['unmapped_too_short'] = int(match.group(1))
    
    match = re.search(r'Number of reads unmapped: other\s+\|\s+(\d+)', content)
    if match:
        stats['unmapped_other'] = int(match.group(1))
    
    # Extract error rates
    match = re.search(r'Mismatch rate per base, %\s+\|\s+([\d.]+)%', content)
    if match:
        stats['mismatch_rate'] = float(match.group(1))
    
    match = re.search(r'Deletion rate per base\s+\|\s+([\d.]+)', content)
    if match:
        stats['deletion_rate'] = float(match.group(1))
    
    match = re.search(r'Insertion rate per base\s+\|\s+([\d.]+)', content)
    if match:
        stats['insertion_rate'] = float(match.group(1))
    
    # Compute total alignment rate
    if stats['total_reads'] > 0:
        total_mapped = stats['uniquely_mapped'] + stats['multi_mapped'] + stats['too_many_mapped']
        stats['total_alignment_rate'] = (total_mapped / stats['total_reads']) * 100
    else:
        stats['total_alignment_rate'] = 0.0
    
    return stats


def parse_umi_dedup_log(log_file: Path) -> Dict:
    """Parse a UMI deduplication log file and extract deduplication statistics"""
    stats = {
        'input_reads': 0,
        'output_reads': 0,
        'deduplication_rate': 0.0,
        'mean_umi_per_position': 0.0,
        'max_umi_per_position': 0
    }
    
    if not log_file.exists():
        return stats
    
    with open(log_file, 'r') as f:
        content = f.read()
    
    # Actual umi_tools output formats
    # Format 1: "INFO Reads: Input Reads: 76487311" (new format)
    match = re.search(r'Reads:\s+Input Reads:\s+(\d+)', content)
    if match:
        stats['input_reads'] = int(match.group(1))
    else:
        # Format 2: "Input: 123 reads" (legacy-compatible format)
        match = re.search(r'Input:\s+(\d+) reads', content)
        if match:
            stats['input_reads'] = int(match.group(1))
    
    # Format 1: "INFO Number of reads out: 2344694" (new format)
    match = re.search(r'Number of reads out:\s+(\d+)', content)
    if match:
        stats['output_reads'] = int(match.group(1))
    else:
        # Format 2: "Output: 123 reads" (legacy-compatible format)
        match = re.search(r'Output:\s+(\d+) reads', content)
        if match:
            stats['output_reads'] = int(match.group(1))
    
    # Compute deduplication rate
    if stats['input_reads'] > 0:
        stats['deduplication_rate'] = ((stats['input_reads'] - stats['output_reads']) / stats['input_reads']) * 100
    
    # Extract mean UMI count per position
    match = re.search(r'Mean number of unique UMIs per position:\s+([\d.]+)', content)
    if match:
        stats['mean_umi_per_position'] = float(match.group(1))
    
    # Extract maximum UMI count per position
    match = re.search(r'Max\. number of unique UMIs per position:\s+(\d+)', content)
    if match:
        stats['max_umi_per_position'] = int(match.group(1))
    
    return stats


def calculate_bam_coverage(bam_file: Path, min_mapq: int = 0) -> Dict:
    """Calculate BAM coverage statistics"""
    stats = {
        'total_reads': 0,
        'mapped_reads': 0,
        'unmapped_reads': 0,
        'total_bases': 0,
        'mapped_bases': 0,
        'mean_coverage': 0.0,
        'mean_read_length': 0.0
    }
    
    if pysam is None:
        print("[WARN] pysam is not installed, unable to calculate BAM coverage")
        return stats
    
    if not bam_file.exists():
        return stats
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            total_bases = 0
            mapped_bases = 0
            read_lengths = []
            
            for read in bam:
                stats['total_reads'] += 1
                
                if read.is_unmapped:
                    stats['unmapped_reads'] += 1
                else:
                    if read.mapping_quality >= min_mapq:
                        stats['mapped_reads'] += 1
                        read_len = read.query_length
                        if read_len:
                            read_lengths.append(read_len)
                            total_bases += read_len
                            mapped_bases += read_len
            
            stats['total_bases'] = total_bases
            stats['mapped_bases'] = mapped_bases
            
            if read_lengths:
                stats['mean_read_length'] = sum(read_lengths) / len(read_lengths)
            
            # Get the total reference length
            ref_length = sum(bam.lengths)
            if ref_length > 0:
                stats['mean_coverage'] = (mapped_bases / ref_length) if mapped_bases > 0 else 0.0
                stats['reference_length'] = ref_length
    except Exception as e:
        print(f"[WARN] Unable to calculate BAM coverage {bam_file}: {e}")
    
    return stats


def calculate_chromosome_coverage(bam_file: Path, min_mapq: int = 0) -> Dict:
    """Calculate non-overlapping covered bases and mean covered-region depth for each chromosome
    
    Use interval merging to efficiently calculate non-overlapping covered regions
    
    Return format:
    {
        'chromosome_stats': {
            'chr1': {
                'covered_bases': 1000000,  # non-overlapping covered bases
                'mean_coverage': 5.2,      # mean coverage across covered regions
                'total_bases': 2000000     # total aligned bases for this chromosome
            },
            ...
        }
    }
    """
    from collections import defaultdict
    
    if pysam is None:
        print("[WARN] pysam is not installed, unable to calculate chromosome coverage")
        return {'chromosome_stats': {}}
    
    if not bam_file.exists():
        return {'chromosome_stats': {}}
    
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            # First pass: collect covered intervals
            chrom_intervals = defaultdict(list)  # {chrom: [(start, end), ...]}
            chrom_total_bases = defaultdict(int)   # {chrom: total_bases}
            
            for read in bam:
                if read.is_unmapped or read.mapping_quality < min_mapq:
                    continue
                
                chrom = bam.get_reference_name(read.reference_id)
                if not chrom:
                    continue
                
                start = read.reference_start
                end = read.reference_end
                
                if start is not None and end is not None and end > start:
                    chrom_intervals[chrom].append((start, end))
                    chrom_total_bases[chrom] += (end - start)
            
            # Second pass: merge intervals and calculate coverage
            result = {}
            for chrom, intervals in chrom_intervals.items():
                if not intervals:
                    continue
                
                # Merge overlapping intervals
                intervals.sort()
                merged = []
                for start, end in intervals:
                    if not merged or merged[-1][1] < start:
                        merged.append([start, end])
                    else:
                        merged[-1][1] = max(merged[-1][1], end)
                
                # Compute non-overlapping covered bases
                covered_bases = sum(end - start for start, end in merged)
                
                # Compute mean coverage across covered regions
                # Use pileup to measure depth at covered positions
                mean_coverage = 0.0
                coverage_sum = 0
                coverage_count = 0
                
                try:
                    # Calculate coverage for each merged interval
                    for start, end in merged:
                        for pileupcolumn in bam.pileup(chrom, start, end, stepper='all'):
                            pos = pileupcolumn.pos
                            if start <= pos < end:
                                coverage_sum += pileupcolumn.n
                                coverage_count += 1
                    
                    if coverage_count > 0:
                        mean_coverage = coverage_sum / coverage_count
                    else:
                        # If pileup fails, approximate depth with total_bases / covered_bases
                        mean_coverage = chrom_total_bases[chrom] / covered_bases if covered_bases > 0 else 0.0
                except Exception as e:
                    # If pileup fails, approximate depth with total_bases / covered_bases
                    mean_coverage = chrom_total_bases[chrom] / covered_bases if covered_bases > 0 else 0.0
                
                result[chrom] = {
                    'covered_bases': covered_bases,
                    'mean_coverage': round(mean_coverage, 2),
                    'total_bases': chrom_total_bases[chrom]
                }
            
            return {'chromosome_stats': result}
    except Exception as e:
        print(f"[WARN] Unable to calculate chromosome coverage {bam_file}: {e}")
        import traceback
        traceback.print_exc()
        return {'chromosome_stats': {}}


def generate_alignment_report(
    sample_name: str,
    sample_dir: Path,
    config: Dict,
    output_file: Optional[Path] = None
) -> Dict:
    """Generate a complete alignment statistics report"""
    report = {
        'sample_name': sample_name,
        'stages': {}
    }
    
    align_dir = sample_dir / 'alignment'
    log_dir = sample_dir / 'logs'
    
    if not align_dir.exists():
        return report
    
    # Iterate over each stage
    for stage_dir in sorted(align_dir.iterdir()):
        if not stage_dir.is_dir():
            continue
        
        stage_name = stage_dir.name
        stage_report = {
            'stage_name': stage_name,
            'aligner': None,
            'alignment_stats': {},
            'umi_dedup_stats': {},
            'coverage_stats': {}
        }
        
        # Determine aligner type
        if (stage_dir / 'aligned.bam').exists():
            stage_report['aligner'] = 'bowtie2'
            bam_file = stage_dir / 'aligned.bam'
            log_file = log_dir / f'bowtie2_{stage_name}.log'
            stage_report['alignment_stats'] = parse_bowtie2_log(log_file)
        elif (stage_dir / 'Aligned.sortedByCoord.out.bam').exists():
            stage_report['aligner'] = 'star'
            bam_file = stage_dir / 'Aligned.sortedByCoord.out.bam'
            # STAR detailed statistics are stored in Log.final.out rather than star_{stage_name}.log
            log_file = stage_dir / 'Log.final.out'
            if not log_file.exists():
                # If Log.final.out is missing, fall back to the log file under logs/
                log_file = log_dir / f'star_{stage_name}.log'
            stage_report['alignment_stats'] = parse_star_log(log_file)
        else:
            continue
        
        # Check whether a deduplicated BAM exists and select the BAM used for coverage calculation
        final_bam = bam_file  # Use the original BAM by default
        if config.get('umi_dedup', {}).get('enable', False):
            if stage_report['aligner'] == 'bowtie2':
                dedup_bam = stage_dir / 'aligned.dedup.bam'
            elif stage_report['aligner'] == 'star':
                dedup_bam = stage_dir / 'Aligned.sortedByCoord.out.dedup.bam'
            else:
                dedup_bam = None
            
            if dedup_bam and dedup_bam.exists():
                umi_log = log_dir / f'umi_dedup_{stage_name}.log'
                stage_report['umi_dedup_stats'] = parse_umi_dedup_log(umi_log)
                final_bam = dedup_bam  # Use the deduplicated BAM
        
        # Calculate coverage statistics using the selected final BAM
        stage_report['coverage_stats'] = calculate_bam_coverage(final_bam)
        
        # Calculate chromosome coverage statistics using the selected final BAM
        chrom_coverage = calculate_chromosome_coverage(final_bam)
        stage_report['chromosome_coverage'] = chrom_coverage.get('chromosome_stats', {})
        
        report['stages'][stage_name] = stage_report
    
    # Save the report
    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2, ensure_ascii=False)
    
    return report


def format_report_text(report: Dict) -> str:
    """Format the report into a human-readable text summary"""
    lines = []
    lines.append("=" * 80)
    lines.append(f"Sample alignment statistics report: {report['sample_name']}")
    lines.append("=" * 80)
    lines.append("")
    
    for stage_name, stage_data in report['stages'].items():
        lines.append(f"### Stage: {stage_name} ({stage_data['aligner']})")
        lines.append("-" * 80)
        
        # Alignment statistics
        align_stats = stage_data.get('alignment_stats', {})
        if align_stats:
            lines.append("\n[Alignment statistics]")
            if stage_data['aligner'] == 'bowtie2':
                lines.append(f"  Total reads: {align_stats.get('total_reads', 0):,}")
                lines.append(f"  Uniquely aligned: {align_stats.get('aligned_once', 0):,} ({align_stats.get('unique_alignment_rate', 0):.2f}%)")
                lines.append(f"  Multi-mapped: {align_stats.get('aligned_more_than_once', 0):,} ({align_stats.get('multi_alignment_rate', 0):.2f}%)")
                lines.append(f"  Unaligned: {align_stats.get('unaligned', 0):,}")
                lines.append(f"  Overall alignment rate: {align_stats.get('overall_alignment_rate', 0):.2f}%")
            elif stage_data['aligner'] == 'star':
                lines.append(f"  Total reads: {align_stats.get('total_reads', 0):,}")
                lines.append(f"  Uniquely aligned: {align_stats.get('uniquely_mapped', 0):,} ({align_stats.get('uniquely_mapped_pct', 0):.2f}%)")
                lines.append(f"  Multi-mapped: {align_stats.get('multi_mapped', 0):,} ({align_stats.get('multi_mapped_pct', 0):.2f}%)")
                lines.append(f"  Mapped to too many loci: {align_stats.get('too_many_mapped', 0):,} ({align_stats.get('too_many_mapped_pct', 0):.2f}%)")
                lines.append(f"  Overall alignment rate: {align_stats.get('total_alignment_rate', 0):.2f}%")
                lines.append(f"  Unaligned (too many mismatches): {align_stats.get('unmapped_too_many_mismatches', 0):,}")
                lines.append(f"  Unaligned (too short): {align_stats.get('unmapped_too_short', 0):,}")
                lines.append(f"  Unaligned (other): {align_stats.get('unmapped_other', 0):,}")
                lines.append(f"  Mismatch rate: {align_stats.get('mismatch_rate', 0):.2f}%")
        
        # UMI deduplication statistics
        umi_stats = stage_data.get('umi_dedup_stats', {})
        if umi_stats and umi_stats.get('input_reads', 0) > 0:
            lines.append("\n[UMI deduplication statistics]")
            lines.append(f"  Input reads: {umi_stats.get('input_reads', 0):,}")
            lines.append(f"  Output reads: {umi_stats.get('output_reads', 0):,}")
            lines.append(f"  Deduplication rate: {umi_stats.get('deduplication_rate', 0):.2f}%")
            lines.append(f"  Mean UMIs per position: {umi_stats.get('mean_umi_per_position', 0):.2f}")
            lines.append(f"  Maximum UMIs per position: {umi_stats.get('max_umi_per_position', 0)}")
        
        # Coverage statistics
        cov_stats = stage_data.get('coverage_stats', {})
        if cov_stats and cov_stats.get('total_reads', 0) > 0:
            lines.append("\n[Coverage statistics]")
            # Check whether coverage is based on a deduplicated BAM
            if umi_stats and umi_stats.get('input_reads', 0) > 0:
                lines.append("  (based on the UMI-deduplicated BAM)")
            else:
                lines.append("  (based on the original aligned BAM)")
            lines.append(f"  Total reads: {cov_stats.get('total_reads', 0):,}")
            lines.append(f"  Mapped reads: {cov_stats.get('mapped_reads', 0):,}")
            lines.append(f"  Unmapped reads: {cov_stats.get('unmapped_reads', 0):,}")
            lines.append(f"  Mean read length: {cov_stats.get('mean_read_length', 0):.2f} bp")
            if 'reference_length' in cov_stats:
                lines.append(f"  Reference length: {cov_stats.get('reference_length', 0):,} bp")
            lines.append(f"  Mean coverage: {cov_stats.get('mean_coverage', 0):.2f}x")
        
        lines.append("")
    
    return "\n".join(lines)


def generate_coverage_summary_table(work_dir: Path, config: Dict, output_file: Optional[Path] = None):
    """Generate a chromosome-level coverage summary table across all samples
    
    Parameters:
        work_dir: Working directory containing per-sample subdirectories
        config: Configuration dictionary
        output_file: Optional output CSV path
    
    Returns:
        pandas DataFrame，A pandas DataFrame containing chromosome-level coverage statistics for all samples
    """
    import pandas as pd
    from collections import defaultdict
    
    # Collect all rows
    all_data = []
    
    # Iterate over all samples
    for sample_name in config.get('samples', {}).keys():
        sample_dir = work_dir / sample_name
        if not sample_dir.exists():
            continue
        
        # Load or generate the report
        report_file = sample_dir / 'reports' / 'alignment_report.json'
        if report_file.exists():
            with open(report_file, 'r') as f:
                report = json.load(f)
        else:
            # Generate the report if it does not exist
            report = generate_alignment_report(sample_name, sample_dir, config)
        
        # Extract chromosome coverage for each stage
        for stage_name, stage_data in report.get('stages', {}).items():
            chrom_coverage = stage_data.get('chromosome_coverage', {})
            
            for chrom, stats in chrom_coverage.items():
                all_data.append({
                    'sample': sample_name,
                    'stage': stage_name,
                    'chromosome': chrom,
                    'covered_bases': stats.get('covered_bases', 0),
                    'mean_coverage': stats.get('mean_coverage', 0.0),
                    'total_bases': stats.get('total_bases', 0)
                })
    
    # Build the DataFrame
    if not all_data:
        # Return an empty DataFrame
        df = pd.DataFrame(columns=['sample', 'stage', 'chromosome', 'covered_bases', 'mean_coverage', 'total_bases'])
    else:
        df = pd.DataFrame(all_data)
        
        # Sort by sample, stage, and chromosome
        df = df.sort_values(['sample', 'stage', 'chromosome'])
    
    # Save to CSV
    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_file, index=False, encoding='utf-8')
        print(f"[INFO] Coverage summary table written: {output_file}")
    
    return df


def generate_umi_summary_table(work_dir: Path, config: Dict, output_file: Optional[Path] = None):
    """Generate a UMI deduplication summary table across all samples
    
    Parameters:
        work_dir: Working directory containing per-sample subdirectories
        config: Configuration dictionary
        output_file: Optional output CSV path
    
    Returns:
        pandas DataFrame，A pandas DataFrame containing UMI deduplication statistics for all samples and stages
    """
    import pandas as pd
    
    # Collect all rows
    all_data = []
    
    # Iterate over all samples
    for sample_name in config.get('samples', {}).keys():
        sample_dir = work_dir / sample_name
        if not sample_dir.exists():
            continue
        
        # Load or generate the report
        report_file = sample_dir / 'reports' / 'alignment_report.json'
        if report_file.exists():
            with open(report_file, 'r') as f:
                report = json.load(f)
        else:
            # Generate the report if it does not exist
            report = generate_alignment_report(sample_name, sample_dir, config)
        
        # Extract UMI deduplication statistics for each stage
        for stage_name, stage_data in report.get('stages', {}).items():
            umi_stats = stage_data.get('umi_dedup_stats', {})
            
            if umi_stats and umi_stats.get('input_reads', 0) > 0:
                # Stage has UMI deduplication statistics
                all_data.append({
                    'sample': sample_name,
                    'stage': stage_name,
                    'umi_input_reads': umi_stats.get('input_reads', 0),
                    'umi_output_reads': umi_stats.get('output_reads', 0),
                    'umi_dedup_rate': round(umi_stats.get('deduplication_rate', 0.0), 2),
                    'umi_mean_per_position': round(umi_stats.get('mean_umi_per_position', 0.0), 2),
                    'umi_max_per_position': umi_stats.get('max_umi_per_position', 0)
                })
            else:
                # No UMI deduplication statistics are available
                all_data.append({
                    'sample': sample_name,
                    'stage': stage_name,
                    'umi_input_reads': None,
                    'umi_output_reads': None,
                    'umi_dedup_rate': None,
                    'umi_mean_per_position': None,
                    'umi_max_per_position': None
                })
    
    # Build the DataFrame
    if not all_data:
        # Return an empty DataFrame
        df = pd.DataFrame(columns=['sample', 'stage', 'umi_input_reads', 'umi_output_reads', 
                                    'umi_dedup_rate', 'umi_mean_per_position', 'umi_max_per_position'])
    else:
        df = pd.DataFrame(all_data)
        
        # Sort by sample and stage
        df = df.sort_values(['sample', 'stage'])
    
    # Save to CSV
    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_file, index=False, encoding='utf-8')
        print(f"[INFO] UMI summary table written: {output_file}")
    
    return df


def generate_sample_summary_table(work_dir: Path, config: Dict, output_file: Optional[Path] = None):
    """Generate a complete summary table across all samples
    
    Includes alignment, UMI, and coverage statistics for each sample and stage
    
    Parameters:
        work_dir: Working directory containing per-sample subdirectories
        config: Configuration dictionary
        output_file: Optional output CSV path
    
    Returns:
        pandas DataFrame，A pandas DataFrame containing the complete per-stage statistics for all samples
    """
    import pandas as pd
    
    # Collect all rows
    all_data = []
    
    # Iterate over all samples
    for sample_name in config.get('samples', {}).keys():
        sample_dir = work_dir / sample_name
        if not sample_dir.exists():
            continue
        
        # Load or generate the report
        report_file = sample_dir / 'reports' / 'alignment_report.json'
        if report_file.exists():
            with open(report_file, 'r') as f:
                report = json.load(f)
        else:
            # Generate the report if it does not exist
            report = generate_alignment_report(sample_name, sample_dir, config)
        
        # Extract all statistics for each stage
        for stage_name, stage_data in report.get('stages', {}).items():
            align_stats = stage_data.get('alignment_stats', {})
            umi_stats = stage_data.get('umi_dedup_stats', {})
            cov_stats = stage_data.get('coverage_stats', {})
            
            # Build the output row
            row = {
                'sample': sample_name,
                'stage': stage_name,
                'aligner': stage_data.get('aligner', '')
            }
            
            # Alignment statistics by aligner type
            if stage_data.get('aligner') == 'bowtie2':
                row['total_reads'] = align_stats.get('total_reads', 0)
                row['unique_aligned'] = align_stats.get('aligned_once', 0)
                row['unique_aligned_pct'] = round(align_stats.get('unique_alignment_rate', 0.0), 2)
                row['multi_aligned'] = align_stats.get('aligned_more_than_once', 0)
                row['multi_aligned_pct'] = round(align_stats.get('multi_alignment_rate', 0.0), 2)
                row['unaligned'] = align_stats.get('unaligned', 0)
                row['alignment_rate'] = round(align_stats.get('overall_alignment_rate', 0.0), 2)
                row['too_many_mapped'] = None
                row['too_many_mapped_pct'] = None
                row['mismatch_rate'] = None
            elif stage_data.get('aligner') == 'star':
                row['total_reads'] = align_stats.get('total_reads', 0)
                row['unique_aligned'] = align_stats.get('uniquely_mapped', 0)
                row['unique_aligned_pct'] = round(align_stats.get('uniquely_mapped_pct', 0.0), 2)
                row['multi_aligned'] = align_stats.get('multi_mapped', 0)
                row['multi_aligned_pct'] = round(align_stats.get('multi_mapped_pct', 0.0), 2)
                row['unaligned'] = (align_stats.get('unmapped_too_many_mismatches', 0) + 
                                   align_stats.get('unmapped_too_short', 0) + 
                                   align_stats.get('unmapped_other', 0))
                row['alignment_rate'] = round(align_stats.get('total_alignment_rate', 0.0), 2)
                row['too_many_mapped'] = align_stats.get('too_many_mapped', 0)
                row['too_many_mapped_pct'] = round(align_stats.get('too_many_mapped_pct', 0.0), 2)
                row['mismatch_rate'] = round(align_stats.get('mismatch_rate', 0.0), 2)
            else:
                row['total_reads'] = None
                row['unique_aligned'] = None
                row['unique_aligned_pct'] = None
                row['multi_aligned'] = None
                row['multi_aligned_pct'] = None
                row['unaligned'] = None
                row['alignment_rate'] = None
                row['too_many_mapped'] = None
                row['too_many_mapped_pct'] = None
                row['mismatch_rate'] = None
            
            # UMI statistics
            if umi_stats and umi_stats.get('input_reads', 0) > 0:
                row['umi_input_reads'] = umi_stats.get('input_reads', 0)
                row['umi_output_reads'] = umi_stats.get('output_reads', 0)
                row['umi_dedup_rate'] = round(umi_stats.get('deduplication_rate', 0.0), 2)
                row['umi_mean_per_position'] = round(umi_stats.get('mean_umi_per_position', 0.0), 2)
                row['umi_max_per_position'] = umi_stats.get('max_umi_per_position', 0)
            else:
                row['umi_input_reads'] = None
                row['umi_output_reads'] = None
                row['umi_dedup_rate'] = None
                row['umi_mean_per_position'] = None
                row['umi_max_per_position'] = None
            
            # Coverage statistics
            if cov_stats and cov_stats.get('total_reads', 0) > 0:
                row['coverage_total_reads'] = cov_stats.get('total_reads', 0)
                row['coverage_mapped_reads'] = cov_stats.get('mapped_reads', 0)
                row['coverage_unmapped_reads'] = cov_stats.get('unmapped_reads', 0)
                row['coverage_mean_read_length'] = round(cov_stats.get('mean_read_length', 0.0), 2)
                row['coverage_mean_coverage'] = round(cov_stats.get('mean_coverage', 0.0), 2)
            else:
                row['coverage_total_reads'] = None
                row['coverage_mapped_reads'] = None
                row['coverage_unmapped_reads'] = None
                row['coverage_mean_read_length'] = None
                row['coverage_mean_coverage'] = None
            
            all_data.append(row)
    
    # Build the DataFrame
    if not all_data:
        # Return an empty DataFrame with all expected columns
        columns = ['sample', 'stage', 'aligner', 
                  'total_reads', 'unique_aligned', 'unique_aligned_pct', 
                  'multi_aligned', 'multi_aligned_pct', 'unaligned', 'alignment_rate',
                  'too_many_mapped', 'too_many_mapped_pct', 'mismatch_rate',
                  'umi_input_reads', 'umi_output_reads', 'umi_dedup_rate', 
                  'umi_mean_per_position', 'umi_max_per_position',
                  'coverage_total_reads', 'coverage_mapped_reads', 'coverage_unmapped_reads',
                  'coverage_mean_read_length', 'coverage_mean_coverage']
        df = pd.DataFrame(columns=columns)
    else:
        df = pd.DataFrame(all_data)
        
        # Sort by sample and stage
        df = df.sort_values(['sample', 'stage'])
    
    # Save to CSV
    if output_file:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_file, index=False, encoding='utf-8')
        print(f"[INFO] Sample summary table written: {output_file}")
    
    return df


if __name__ == '__main__':
    # Test block
    import sys
    if len(sys.argv) > 1:
        sample_dir = Path(sys.argv[1])
        sample_name = sample_dir.name
        config = {'umi_dedup': {'enable': True}}
        report = generate_alignment_report(sample_name, sample_dir, config)
        print(format_report_text(report))

