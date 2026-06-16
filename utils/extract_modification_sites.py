#!/usr/bin/env python3
"""
BaseCounter: BAM base-counting tool
For RNA modification detection (PUMseq/Gloriseq/Limeseq)

Supports two usage modes:
1. Import as a module: from utils.extract_modification_sites import BaseCounter
2. Run from the command line: python utils/extract_modification_sites.py --bam ... --fasta ...
"""

import os
import argparse
import pysam
import pandas as pd
from collections import defaultdict
from typing import List, Set, Optional, Dict
from enum import Enum


# ========= Sequencing method configuration =========
class SeqMethod(Enum):
    """Enumeration of sequencing methods."""
    PUMSEQ = "puMseq"      # Keep the case consistent with the configuration file
    GLORISEQ = "Gloriseq"
    LIMESEQ = "Limeseq"
    CUSTOM = "custom"


class SeqMethodConfig:
    """Configuration for a sequencing method."""
    def __init__(
        self,
        name: str,
        ref_bases: Set[str],
        description: str = ""
    ):
        self.name = name
        self.ref_bases = set(b.upper() for b in ref_bases)
        self.description = description

    def should_include(self, ref_base: str) -> bool:
        """Return whether this reference base should be included."""
        return ref_base.upper() in self.ref_bases


# ========= Predefined configurations =========
SEQMETHOD_CONFIGS = {
    SeqMethod.PUMSEQ: SeqMethodConfig(
        name="PUMseq",
        ref_bases={"T"},
        description="Only report reference T sites (detecting T>C conversion for pseudouridine-related analysis)"
    ),
    SeqMethod.GLORISEQ: SeqMethodConfig(
        name="Gloriseq",
        ref_bases={"A"},
        description="Only report reference A sites (detecting A>G edits for m6A-related analysis)"
    ),
    SeqMethod.LIMESEQ: SeqMethodConfig(
        name="Limeseq",
        ref_bases={"A", "C", "G", "T"},
        description="Report all four reference bases (broad modification screening)"
    ),
}


# ========= Main class =========
class BaseCounter:
    """Main base-counting class

    Note: orient_to_reference is fixed to False
          Bases on reverse-strand reads are kept in their original read orientation instead of being converted to the forward reference strand
    """

    def __init__(
        self,
        bam_path: str,
        ref_fasta_path: str,
        seq_method: SeqMethod = SeqMethod.LIMESEQ,
        custom_ref_bases: Optional[Set[str]] = None,
        # Filtering parameters
        min_mapq: int = 10,
        min_bq: int = 20,
        min_depth: int = 10,
        min_freq: float = 0.0,  # Minimum modification frequency
        max_depth: int = 1_000_000,
        ignore_overlaps: bool = True,
        include_duplicates: bool = False,
        include_secondary: bool = False,
        include_supp: bool = False,
        merge_strands: bool = False,  # Default False to keep strand-separated output, matching the main pipeline
        # Runtime parameters
        limit: Optional[int] = None,
        threads: int = 4,
        verbose: bool = True,
        streaming: bool = True  # Streaming output mode (default True to reduce memory usage)
    ):
        self.bam_path = bam_path
        self.ref_fasta_path = ref_fasta_path
        self.min_mapq = min_mapq
        self.min_bq = min_bq
        self.min_depth = min_depth
        self.min_freq = min_freq
        self.max_depth = max_depth
        self.ignore_overlaps = ignore_overlaps
        self.include_duplicates = include_duplicates
        self.include_secondary = include_secondary
        self.include_supp = include_supp
        self.orient_to_reference = False  # Fixed to False
        self.merge_strands = merge_strands
        self.limit = limit
        self.threads = threads
        self.verbose = verbose
        self.streaming = streaming

        # Set the sequencing method configuration
        if seq_method == SeqMethod.CUSTOM:
            if custom_ref_bases is None:
                raise ValueError("CUSTOM mode requires the custom_ref_bases argument")
            self.config = SeqMethodConfig("Custom", custom_ref_bases)
        else:
            self.config = SEQMETHOD_CONFIGS[seq_method]

        # Base complement map kept for possible future use
        self.comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}

        # Validate input files
        self._validate_files()

    def _get_rna_ref_base(self, ref_base: str, strand: str) -> str:
        """Return the RNA-oriented reference base for a genomic strand label.

        In strand-separated mode, `strand` is reported in RNA orientation:
        - `+`: RNA sequence follows the FASTA forward strand.
        - `-`: RNA sequence follows the FASTA reverse-complement strand.
        """
        ref_base = ref_base.upper()
        if strand == "-":
            return self.comp.get(ref_base, "N")
        return ref_base

    def _read_base_to_rna_base(self, read_base: str) -> str:
        """Convert a read-oriented cDNA base to RNA base.

        For reverse-stranded RNA-seq libraries, sequenced cDNA bases are complements
        of original RNA bases. Therefore RNA-base counts are obtained by complementing
        read-oriented bases before accumulation.
        """
        return self.comp.get(read_base.upper(), "N")

    def _get_rna_strand(self, aln) -> str:
        """Map alignment orientation to RNA strand.

        For reverse-stranded libraries:
        - reverse alignments (is_reverse=True) originate from RNA `+`
        - forward alignments (is_reverse=False) originate from RNA `-`
        """
        return "+" if aln.is_reverse else "-"

    def _get_output_ref_base(self, ref_base: str) -> str:
        """Return the reference base written to the output `ref` column.

        The output `ref` value always keeps the original FASTA base to stay aligned with the reference coordinate system.
        """
        return ref_base.upper()

    def _get_read_oriented_base_and_quality(self, aln, qpos):
        """Return the base and quality value in the original read orientation.

        `pysam` reports `query_sequence` in reference-oriented space for reverse-strand alignments.
        This converts the sequence back to the original read orientation so reverse-strand counts are not complemented before counting.
        """
        if qpos is None or aln.query_sequence is None:
            return None, None

        quals = aln.query_qualities
        if aln.is_reverse:
            read_len = len(aln.query_sequence)
            read_qpos = read_len - 1 - qpos
            base = aln.get_forward_sequence()[read_qpos].upper()
            q = (quals[read_qpos] if quals is not None else 255)
        else:
            base = aln.query_sequence[qpos].upper()
            q = (quals[qpos] if quals is not None else 255)

        if base not in ("A", "C", "G", "T", "N"):
            base = "N"
        return base, q

    def _validate_files(self):
        """Validate input files and create indexes when needed."""
        assert os.path.exists(self.bam_path), f"BAM not found: {self.bam_path}"
        assert os.path.exists(self.ref_fasta_path), f"Reference FASTA not found: {self.ref_fasta_path}"

        if not os.path.exists(self.bam_path + ".bai"):
            if self.verbose:
                print(f"[INFO] Indexing BAM: {self.bam_path}")
            pysam.index(self.bam_path)

        if not os.path.exists(self.ref_fasta_path + ".fai"):
            if self.verbose:
                print(f"[INFO] Indexing FASTA: {self.ref_fasta_path}")
            pysam.faidx(self.ref_fasta_path)

    def _get_contigs(self, bam):
        """Return contigs that have mapped coverage."""
        try:
            stats = bam.get_index_statistics()
            contigs = []
            for s in stats:
                if s.mapped > 0 and s.tid >= 0:
                    name = bam.get_reference_name(s.tid)
                    length = bam.get_reference_length(name)
                    contigs.append((name, length))
        except Exception:
            contigs = list(zip(bam.references, bam.lengths))
        return contigs

    def _should_skip_alignment(self, aln):
        """Return whether this alignment should be skipped."""
        if (not self.include_secondary) and aln.is_secondary:
            return True
        if (not self.include_supp) and aln.is_supplementary:
            return True
        if (not self.include_duplicates) and aln.is_duplicate:
            return True
        if aln.mapping_quality < self.min_mapq:
            return True
        return False

    def run(self, output_path: Optional[str] = None) -> pd.DataFrame:
        """Run base counting
        
        Args:
            output_path: optional output path. If provided with streaming=True, rows are written directly to disk
            
        Returns:
            DataFrame: if streaming=False or output_path is not provided, return a DataFrame;
                      if streaming=True and output_path is provided, return an empty summary DataFrame
        """
        if self.verbose:
            print(f"[INFO] ========== Starting BaseCounter ==========")
            print(f"[INFO] Sequencing method: {self.config.name}")
            print(f"[INFO] Description: {self.config.description}")
            print(f"[INFO] Filtered reference bases: {sorted(self.config.ref_bases)}")
            print(f"[INFO] Merge strands: {'yes' if self.merge_strands else 'no'}")
            print(f"[INFO] Base orientation conversion: no (orient_to_reference=False)")
            print(f"[INFO] Minimum depth: {self.min_depth}, Minimum modification frequency: {self.min_freq}")
            print(f"[INFO] Output mode: {'streaming' if (self.streaming and output_path) else 'in-memory'}")
            print(f"[INFO] BAM: {self.bam_path}")
            print(f"[INFO] Reference: {self.ref_fasta_path}")
            print("[INFO] " + "=" * 40)

        # Streaming mode: write rows directly to the output file
        if self.streaming and output_path:
            return self._run_streaming(output_path)
        else:
            return self._run_in_memory()
    
    def _run_streaming(self, output_path: str) -> pd.DataFrame:
        """Streaming mode: write rows directly to the output file to save memory."""
        import gzip
        
        # Determine whether compression is needed
        compression = output_path.endswith('.gz')
        sep = '\t' if output_path.endswith('.tsv') or output_path.endswith('.tsv.gz') else ','
        
        # Open the output file
        if compression:
            fout = gzip.open(output_path, 'wt', encoding='utf-8')
        else:
            fout = open(output_path, 'w', encoding='utf-8')
        
        try:
            # Write the header
            if self.merge_strands:
                header = sep.join(["chrom", "pos", "ref", "depth", "A", "C", "G", "T", "N"])
            else:
                header = sep.join(["chrom", "pos", "ref", "strand", "depth", "A", "C", "G", "T", "N"])
            fout.write(header + "\n")
            
            bam = pysam.AlignmentFile(self.bam_path, "rb", threads=self.threads)
            fasta = pysam.FastaFile(self.ref_fasta_path)
            
            contigs = self._get_contigs(bam)
            total_positions = 0
            filtered_positions = 0
            low_depth_positions = 0
            low_freq_positions = 0
            output_rows = 0
            
            for chrom, length in contigs:
                if self.verbose:
                    print(f"[INFO] Processing: {chrom} (length: {length:,})")
                
                for col in bam.pileup(
                    chrom, 0, length,
                    truncate=True,
                    stepper="samtools",
                    min_base_quality=self.min_bq,
                    max_depth=self.max_depth,
                    ignore_overlaps=self.ignore_overlaps
                ):
                    total_positions += 1
                    
                    # Fetch the reference base
                    pos0 = col.pos
                    try:
                        ref_base = fasta.fetch(chrom, pos0, pos0 + 1).upper()
                        if not ref_base:
                            ref_base = "N"
                    except Exception:
                        ref_base = "N"
                    
                    # Determine counting mode based on merge_strands
                    if self.merge_strands:
                        # Merged mode cannot distinguish RNA strand orientation, so reference-forward bases are used for filtering.
                        if not self.config.should_include(ref_base):
                            filtered_positions += 1
                            continue
                        counts = defaultdict(int)
                        depth = 0
                    else:
                        counts_plus = defaultdict(int)
                        counts_minus = defaultdict(int)
                        depth_plus = 0
                        depth_minus = 0
                    
                    # Iterate over all pileup reads
                    for pr in col.pileups:
                        if pr.is_refskip or pr.is_del:
                            continue
                        
                        aln = pr.alignment
                        if self._should_skip_alignment(aln):
                            continue
                        
                        qpos = pr.query_position
                        if qpos is None or aln.query_sequence is None:
                            continue
                        
                        base, q = self._get_read_oriented_base_and_quality(aln, qpos)
                        if base is None:
                            continue
                        if q is None or q < self.min_bq:
                            continue
                        
                        # Update counts according to the selected mode
                        if self.merge_strands:
                            counts[base] += 1
                            depth += 1
                        else:
                            rna_base = self._read_base_to_rna_base(base)
                            rna_strand = self._get_rna_strand(aln)
                            if rna_strand == "+":
                                counts_plus[rna_base] += 1
                                depth_plus += 1
                            else:
                                counts_minus[rna_base] += 1
                                depth_minus += 1
                    
                    # Write results according to the selected mode
                    if self.merge_strands:
                        if depth >= self.min_depth:
                            ref_count = counts[ref_base]
                            alt_freq = (depth - ref_count) / depth if depth > 0 else 0
                            if alt_freq >= self.min_freq:
                                row = sep.join([
                                    str(chrom), str(pos0 + 1), str(ref_base), str(depth),
                                    str(counts["A"]), str(counts["C"]), str(counts["G"]),
                                    str(counts["T"]), str(counts["N"])
                                ])
                                fout.write(row + "\n")
                                output_rows += 1
                            else:
                                low_freq_positions += 1
                        elif depth > 0:
                            low_depth_positions += 1
                    else:
                        # Positive strand
                        if depth_plus >= self.min_depth:
                            candidate_base_plus = self._get_rna_ref_base(ref_base, "+")
                            output_ref_plus = self._get_output_ref_base(ref_base)
                            if not self.config.should_include(candidate_base_plus):
                                filtered_positions += 1
                            else:
                                ref_count = counts_plus[candidate_base_plus]
                                alt_freq = (depth_plus - ref_count) / depth_plus if depth_plus > 0 else 0
                                if alt_freq >= self.min_freq:
                                    row = sep.join([
                                        str(chrom), str(pos0 + 1), str(output_ref_plus), "+", str(depth_plus),
                                        str(counts_plus["A"]), str(counts_plus["C"]),
                                        str(counts_plus["G"]), str(counts_plus["T"]), str(counts_plus["N"])
                                    ])
                                    fout.write(row + "\n")
                                    output_rows += 1
                                else:
                                    low_freq_positions += 1
                        elif depth_plus > 0:
                            low_depth_positions += 1
                        
                        # Negative strand
                        if depth_minus >= self.min_depth:
                            candidate_base_minus = self._get_rna_ref_base(ref_base, "-")
                            output_ref_minus = self._get_output_ref_base(ref_base)
                            if not self.config.should_include(candidate_base_minus):
                                filtered_positions += 1
                            else:
                                ref_count = counts_minus[candidate_base_minus]
                                alt_freq = (depth_minus - ref_count) / depth_minus if depth_minus > 0 else 0
                                if alt_freq >= self.min_freq:
                                    row = sep.join([
                                        str(chrom), str(pos0 + 1), str(output_ref_minus), "-", str(depth_minus),
                                        str(counts_minus["A"]), str(counts_minus["C"]),
                                        str(counts_minus["G"]), str(counts_minus["T"]), str(counts_minus["N"])
                                    ])
                                    fout.write(row + "\n")
                                    output_rows += 1
                                else:
                                    low_freq_positions += 1
                        elif depth_minus > 0:
                            low_depth_positions += 1
                    
                    # Apply row limit if requested
                    if self.limit is not None and output_rows >= self.limit:
                        if self.verbose:
                            print(f"[INFO] Reached row limit {self.limit}, stopping")
                        break
                
                if self.limit is not None and output_rows >= self.limit:
                    break
            
            bam.close()
            fasta.close()
            fout.close()
            
            if self.verbose:
                print("[INFO] " + "=" * 40)
                print(f"[INFO] Total scanned positions: {total_positions:,}")
                print(f"[INFO] Filtered by reference base: {filtered_positions:,}")
                print(f"[INFO] Filtered by low depth: {low_depth_positions:,} (depth < {self.min_depth})")
                print(f"[INFO] Filtered by low frequency: {low_freq_positions:,} (freq < {self.min_freq})")
                print(f"[INFO] Output rows: {output_rows:,}")
                file_size = os.path.getsize(output_path)
                print(f"[INFO] Output file size: {file_size:,} bytes ({file_size/1024/1024:.2f} MB)")
            
            # Return an empty summary DataFrame to keep the interface consistent
            return pd.DataFrame({
                'total_positions': [total_positions],
                'filtered_positions': [filtered_positions],
                'low_depth_positions': [low_depth_positions],
                'low_freq_positions': [low_freq_positions],
                'output_rows': [output_rows],
                'output_path': [output_path]
            })
            
        except Exception as e:
            fout.close()
            raise e
    
    def _run_in_memory(self) -> pd.DataFrame:
        """In-memory mode: collect all results in memory and return a DataFrame."""
        bam = pysam.AlignmentFile(self.bam_path, "rb", threads=self.threads)
        fasta = pysam.FastaFile(self.ref_fasta_path)

        contigs = self._get_contigs(bam)
        rows = []
        total_positions = 0
        filtered_positions = 0
        low_depth_positions = 0
        low_freq_positions = 0

        for chrom, length in contigs:
            if self.verbose:
                print(f"[INFO] Processing: {chrom} (length: {length:,})")

            for col in bam.pileup(
                chrom, 0, length,
                truncate=True,
                stepper="samtools",
                min_base_quality=self.min_bq,
                max_depth=self.max_depth,
                ignore_overlaps=self.ignore_overlaps
            ):
                total_positions += 1

                # Fetch the reference base
                pos0 = col.pos
                try:
                    ref_base = fasta.fetch(chrom, pos0, pos0 + 1).upper()
                    if not ref_base:
                        ref_base = "N"
                except Exception:
                    ref_base = "N"

                # Determine counting mode based on merge_strands
                if self.merge_strands:
                    # Merged mode cannot distinguish RNA strand orientation, so reference-forward bases are used for filtering.
                    if not self.config.should_include(ref_base):
                        filtered_positions += 1
                        continue
                    counts = defaultdict(int)
                    depth = 0
                else:
                    counts_plus = defaultdict(int)
                    counts_minus = defaultdict(int)
                    depth_plus = 0
                    depth_minus = 0

                # Iterate over all pileup reads
                for pr in col.pileups:
                    if pr.is_refskip or pr.is_del:
                        continue

                    aln = pr.alignment
                    if self._should_skip_alignment(aln):
                        continue

                    qpos = pr.query_position
                    if qpos is None or aln.query_sequence is None:
                        continue

                    base, q = self._get_read_oriented_base_and_quality(aln, qpos)
                    if base is None:
                        continue
                    if q is None or q < self.min_bq:
                        continue

                    # Update counts according to the selected mode
                    if self.merge_strands:
                        counts[base] += 1
                        depth += 1
                    else:
                        rna_base = self._read_base_to_rna_base(base)
                        rna_strand = self._get_rna_strand(aln)
                        if rna_strand == "+":
                            counts_plus[rna_base] += 1
                            depth_plus += 1
                        else:
                            counts_minus[rna_base] += 1
                            depth_minus += 1

                # Emit results according to the selected mode with depth and frequency filters
                if self.merge_strands:
                    # Merged mode: one row per site
                    if depth >= self.min_depth:
                        ref_count = counts[ref_base]
                        alt_freq = (depth - ref_count) / depth if depth > 0 else 0
                        if alt_freq >= self.min_freq:
                            rows.append((
                                chrom, pos0 + 1, ref_base, depth,
                                counts["A"], counts["C"], counts["G"],
                                counts["T"], counts["N"]
                            ))
                        else:
                            low_freq_positions += 1
                    elif depth > 0:
                        low_depth_positions += 1
                else:
                    # Strand-separated mode: up to two rows per site
                    if depth_plus == 0 and depth_minus == 0:
                        continue

                    # Positive-strand depth and frequency filtering
                    if depth_plus >= self.min_depth:
                        candidate_base_plus = self._get_rna_ref_base(ref_base, "+")
                        output_ref_plus = self._get_output_ref_base(ref_base)
                        if not self.config.should_include(candidate_base_plus):
                            filtered_positions += 1
                        else:
                            ref_count = counts_plus[candidate_base_plus]
                            alt_freq = (depth_plus - ref_count) / depth_plus if depth_plus > 0 else 0
                            if alt_freq >= self.min_freq:
                                rows.append((
                                    chrom, pos0 + 1, output_ref_plus, "+", depth_plus,
                                    counts_plus["A"], counts_plus["C"],
                                    counts_plus["G"], counts_plus["T"], counts_plus["N"]
                                ))
                            else:
                                low_freq_positions += 1
                    elif depth_plus > 0:
                        low_depth_positions += 1

                    # Negative-strand depth and frequency filtering
                    if depth_minus >= self.min_depth:
                        candidate_base_minus = self._get_rna_ref_base(ref_base, "-")
                        output_ref_minus = self._get_output_ref_base(ref_base)
                        if not self.config.should_include(candidate_base_minus):
                            filtered_positions += 1
                        else:
                            ref_count = counts_minus[candidate_base_minus]
                            alt_freq = (depth_minus - ref_count) / depth_minus if depth_minus > 0 else 0
                            if alt_freq >= self.min_freq:
                                rows.append((
                                    chrom, pos0 + 1, output_ref_minus, "-", depth_minus,
                                    counts_minus["A"], counts_minus["C"],
                                    counts_minus["G"], counts_minus["T"], counts_minus["N"]
                                ))
                            else:
                                low_freq_positions += 1
                    elif depth_minus > 0:
                        low_depth_positions += 1

                # Apply row limit if requested
                if self.limit is not None and len(rows) >= self.limit:
                    if self.verbose:
                        print(f"[INFO] Reached row limit {self.limit}, stopping")
                    break

            if self.limit is not None and len(rows) >= self.limit:
                break

        bam.close()
        fasta.close()

        if self.verbose:
            print("[INFO] " + "=" * 40)
            print(f"[INFO] Total scanned positions: {total_positions:,}")
            print(f"[INFO] Filtered by reference base: {filtered_positions:,}")
            print(f"[INFO] Filtered by low depth: {low_depth_positions:,} (depth < {self.min_depth})")
            print(f"[INFO] Filtered by low frequency: {low_freq_positions:,} (freq < {self.min_freq})")
            print(f"[INFO] Output rows: {len(rows):,}")

        # Build the DataFrame according to the selected output mode
        if self.merge_strands:
            df = pd.DataFrame(
                rows,
                columns=["chrom", "pos", "ref", "depth", "A", "C", "G", "T", "N"]
            )
        else:
            df = pd.DataFrame(
                rows,
                columns=["chrom", "pos", "ref", "strand", "depth", "A", "C", "G", "T", "N"]
            )

        return df

    def save_to_csv(self, df: pd.DataFrame, output_path: str, compression: Optional[str] = None):
        """Save results to CSV/TSV (mainly for in-memory mode).
        
        Note: if streaming mode is used, pass output_path directly to run() instead of calling this method,
        so this method is not needed.
        """
        if self.verbose:
            print(f"[INFO] Saving results to: {output_path}")
        
        # Determine the separator from the file extension
        sep = '\t' if output_path.endswith('.tsv') or output_path.endswith('.tsv.gz') else ','
        
        df.to_csv(output_path, sep=sep, index=False, compression=compression)
        
        if self.verbose:
            file_size = os.path.getsize(output_path)
            print(f"[INFO] Done! File size: {file_size:,} bytes ({file_size/1024/1024:.2f} MB)")


# ========================================================================
# Command-line interface (compatible with the main pipeline)
# ========================================================================
def main():
    """CLI entry point compatible with the main pipeline argument format."""
    parser = argparse.ArgumentParser(
        description='Extract modification sites from BAM (BaseCounter)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # PUMseq example (reference T sites only)
  python extract_modification_sites.py --bam in.bam --fasta ref.fa \\
    --mod-type puMseq --min-depth 10 --min-freq 0.05 --output out.tsv

  # Limeseq example (all reference bases)
  python extract_modification_sites.py --bam in.bam --fasta ref.fa \\
    --mod-type Limeseq --min-depth 10 --output out.tsv.gz

  # Merge strands
  python extract_modification_sites.py --bam in.bam --fasta ref.fa \\
    --mod-type puMseq --merge-strands --output out.tsv
        """
    )
    
    # Required arguments
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--fasta', required=True, help='Reference FASTA file')
    parser.add_argument('--mod-type', required=True, 
                        choices=['puMseq', 'Gloriseq', 'Limeseq', 'custom'],
                        help='Modification type')
    parser.add_argument('--custom-ref-bases', help='Custom reference base set, comma-separated (requires --mod-type custom)')
    parser.add_argument('--output', required=True, help='Output TSV file (supports .gz compression)')
    
    # Threshold parameters
    parser.add_argument('--min-depth', type=int, default=10, help='Minimum sequencing depth (default: 10)')
    parser.add_argument('--min-freq', type=float, default=0.0, help='Minimum modification frequency (default: 0.0)')
    
    # Filtering parameters
    parser.add_argument('--min-mapq', type=int, default=10, help='Minimum mapping quality (default: 10)')
    parser.add_argument('--min-bq', type=int, default=20, help='Minimum base quality (default: 20)')
    parser.add_argument('--max-depth', type=int, default=1_000_000, help='Maximum pileup depth (default: 1000000)')
    parser.add_argument('--include-duplicates', action='store_true', help='Include duplicate reads (default: False)')
    parser.add_argument('--include-secondary', action='store_true', help='Include secondary alignments (default: False)')
    parser.add_argument('--include-supplementary', action='store_true', help='Include supplementary alignments (default: False)')
    parser.add_argument('--ignore-overlaps', action='store_true', default=True, help='Ignore overlapping paired reads (default: True)')
    
    # Output mode
    parser.add_argument('--merge-strands', action='store_true', help='Merge strands (default: separate output)')
    
    # Debug parameters
    parser.add_argument('--limit', type=int, help='Limit output rows for testing (default: unlimited)')
    parser.add_argument('--threads', type=int, default=4, help='BAM read threads (default: 4)')
    parser.add_argument('--quiet', action='store_true', help='Quiet mode (suppress progress output)')
    
    args = parser.parse_args()

    # Map mod-type values to SeqMethod
    mod_type_map = {
        'puMseq': SeqMethod.PUMSEQ,
        'Gloriseq': SeqMethod.GLORISEQ,
        'Limeseq': SeqMethod.LIMESEQ,
        'custom': SeqMethod.CUSTOM,
    }

    custom_ref_bases = None
    if args.custom_ref_bases:
        custom_ref_bases = {
            b.strip().upper()
            for b in args.custom_ref_bases.split(',')
            if b.strip()
        }
        if not custom_ref_bases:
            parser.error('Parsed custom reference base set is empty; please check --custom-ref-bases')

    if args.mod_type == 'custom' and not custom_ref_bases:
        parser.error('--mod-type custom also requires --custom-ref-bases')
    if args.mod_type != 'custom' and custom_ref_bases:
        print('[WARN] Detected --custom-ref-bases; switching to custom mode automatically')
        args.mod_type = 'custom'

    counter = BaseCounter(
        bam_path=args.bam,
        ref_fasta_path=args.fasta,
        seq_method=mod_type_map[args.mod_type],
        custom_ref_bases=custom_ref_bases,
        min_mapq=args.min_mapq,
        min_bq=args.min_bq,
        min_depth=args.min_depth,
        min_freq=args.min_freq,
        max_depth=args.max_depth,
        ignore_overlaps=args.ignore_overlaps,
        include_duplicates=args.include_duplicates,
        include_secondary=args.include_secondary,
        include_supp=args.include_supplementary,
        merge_strands=args.merge_strands,
        limit=args.limit,
        threads=args.threads,
        verbose=not args.quiet
    )
    
    # Run counting in streaming mode to reduce memory usage
    df = counter.run(output_path=args.output)
    
    print(f"[OK] Done! Output file: {args.output}")


# ========================================================================
# Standalone example block for direct execution
# ========================================================================
if __name__ == "__main__":
    import sys
    
    # If command-line arguments are provided, run CLI mode
    if len(sys.argv) > 1:
        main()
        sys.exit(0)

    # Otherwise run the example block
    print("=" * 70)
    print("BaseCounter - standalone example")
    print("Tip: use --help to view command-line options")
    print("=" * 70)

    # Example configuration (replace with real paths)
    BAM_PATH = "/home/tangyh/project/new_pu/work_dir2/D-D/ProcessedData/sncRNA_alignment/aligned.bam"
    REF_FASTA_PATH = "/home/tangyh/reference/db/reference_sequence/Homo_sapiens.GRCh38.sncRNA.fa"
    OUTPUT_DIR = "./results"

    # Check whether the example files exist
    if not os.path.exists(BAM_PATH):
        print(f"\n[ERROR] Example BAM file does not exist: {BAM_PATH}")
        print("Please update BAM_PATH and REF_FASTA_PATH to real paths")
        print("\nOr use command-line mode:")
        print("  python extract_modification_sites.py --bam <bam> --fasta <fasta> --mod-type puMseq --output out.tsv")
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Example 1: PUMseq analysis
    print("\n" + "=" * 70)
    print("[Example] PUMseq analysis - strand-separated output with minimum depth 20")
    print("=" * 70)

    counter = BaseCounter(
        bam_path=BAM_PATH,
        ref_fasta_path=REF_FASTA_PATH,
        seq_method=SeqMethod.PUMSEQ,
        min_mapq=20,
        min_bq=20,
        min_depth=20,
        min_freq=0.05,
        merge_strands=False,
        limit=10000  # Limit to 10000 rows for quick testing
    )

    # Use streaming output (recommended to save memory)
    output_file = os.path.join(OUTPUT_DIR, "pumseq_test.tsv.gz")
    df = counter.run(output_path=output_file)
    
    print(f"\nDone! Results saved to: {output_file}")
    print(f"Output rows: {df['output_rows'].values[0]:,}")
