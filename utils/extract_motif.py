#!/usr/bin/env python3
"""
Add motif information to a site table.

This script uses a reference FASTA file to extract upstream and downstream sequence
windows for each row based on `chrom` and `pos` (1-based coordinates).
"""
import argparse
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam


COMP_BASE = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N',
}


def reverse_complement(seq: str) -> str:
    """Return the reverse-complement DNA sequence."""
    return ''.join(COMP_BASE.get(base, 'N') for base in reversed(seq.upper()))


def get_rna_ref_base(ref_base: str, strand: Optional[str]) -> str:
    """Convert a FASTA-space reference base into its RNA-oriented reference base."""
    ref_base = str(ref_base).upper()
    if strand == '+':
        return COMP_BASE.get(ref_base, 'N')
    return ref_base


def add_motif_column(
    df: pd.DataFrame,
    ref_fasta_path: str,
    n: int = 2,
    m: int = 2,
    motif_col_name: str = 'motif',
    chrom_col: str = 'chrom',
    pos_col: str = 'pos',
    strand_col: str = 'strand',
) -> pd.DataFrame:
    """Add a motif column to a site table based on the reference sequence.

    Notes
    -----
    - `pos` is treated as a 1-based genomic coordinate.
    - If a `strand` column exists, motifs are reported in RNA orientation:
      `+` rows use the reverse complement of the reference window, while `-` rows
      use the reference window directly.
    - If the window extends outside sequence boundaries or an error occurs, `N`
      padding is used.
    """

    print('Starting motif annotation...')
    print(f'Input rows: {len(df):,}')
    print(f'Motif length: {n + 1 + m}-mer')

    if df.empty:
        print('Input table is empty; returning immediately.')
        return df.copy()

    for required_col, alias in ((chrom_col, 'chrom'), (pos_col, 'pos')):
        if required_col not in df.columns:
            raise KeyError(f"Missing required column `{required_col}` (default alias: {alias}).")

    motif_len = n + 1 + m
    default_motif = 'N' * motif_len
    df_result = df.copy()
    has_strand = strand_col in df_result.columns

    if has_strand:
        print(f"Detected strand column `{strand_col}`; motifs will be reported in RNA orientation")
    else:
        print(f"Strand column `{strand_col}` was not found; using the reference forward-strand motif")

    fasta_path = Path(ref_fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f'Reference FASTA not found: {ref_fasta_path}')

    def extract_motif(chrom: str, pos: Optional[int], strand: Optional[str] = None) -> str:
        if pd.isna(chrom) or pd.isna(pos):
            return default_motif
        try:
            pos_int = int(pos)
            if pos_int <= 0:
                return default_motif
        except (TypeError, ValueError):
            return default_motif

        try:
            chrom_length = fasta.get_reference_length(str(chrom))
            pos_0based = pos_int - 1
            start = max(0, pos_0based - n)
            end = min(chrom_length, pos_0based + m + 1)
            seq = fasta.fetch(str(chrom), start, end).upper()
            if len(seq) != motif_len:
                seq = seq.ljust(motif_len, 'N')
            if strand == '+':
                return reverse_complement(seq)
            return seq
        except (KeyError, ValueError, OSError):
            return default_motif

    with pysam.FastaFile(str(fasta_path)) as fasta:
        motifs = [
            extract_motif(
                row[chrom_col],
                row[pos_col],
                row[strand_col] if has_strand else None,
            )
            for _, row in df_result.iterrows()
        ]

    df_result[motif_col_name] = motifs

    if has_strand and 'ref' in df_result.columns and 'rna_ref' not in df_result.columns:
        df_result['rna_ref'] = [
            get_rna_ref_base(row['ref'], row[strand_col])
            for _, row in df_result.iterrows()
        ]

    if 'ref' in df_result.columns:
        expected_len = motif_len
        df_check = df_result[df_result[motif_col_name].str.len() == expected_len].copy()
        if not df_check.empty:
            df_check['motif_center'] = df_check[motif_col_name].str[n]
            if has_strand:
                df_check['rna_ref_check'] = [
                    get_rna_ref_base(row['ref'], row[strand_col])
                    for _, row in df_check.iterrows()
                ]
                match_rate = (
                    (df_check['rna_ref_check'] == df_check['motif_center']).sum()
                    / len(df_check)
                    * 100
                )
                print(f'Motif center matches RNA ref: {match_rate:.1f}%')
            else:
                match_rate = (
                    (df_check['ref'] == df_check['motif_center']).sum()
                    / len(df_check)
                    * 100
                )
                print(f'Motif center matches ref: {match_rate:.1f}%')

    return df_result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Add a motif column to a site table')
    parser.add_argument('--input', required=True, help='Input TSV/CSV file containing chrom and pos columns (1-based)')
    parser.add_argument('--fasta', required=True, help='Reference FASTA path')
    parser.add_argument('--upstream', type=int, default=2, help='Number of upstream motif bases (default: 2)')
    parser.add_argument('--downstream', type=int, default=2, help='Number of downstream motif bases (default: 2)')
    parser.add_argument('--chrom-col', default='chrom', help='Chromosome column name (default: chrom)')
    parser.add_argument('--pos-col', default='pos', help='Position column name, 1-based (default: pos)')
    parser.add_argument('--motif-col', default='motif', help='Output motif column name (default: motif)')
    parser.add_argument('--strand-col', default='strand', help='Strand column name (default: strand)')
    parser.add_argument('--sep', default='	', help='Input/output delimiter (default: tab)')
    parser.add_argument('--inplace', action='store_true', help='Overwrite the input file')
    parser.add_argument(
        '--output',
        help='Output path; if omitted and inplace is disabled, a .motif suffix is added automatically',
    )
    return parser.parse_args()


def main():
    args = parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        raise FileNotFoundError(f'Input file not found: {args.input}')

    df = pd.read_csv(input_path, sep=args.sep)
    df_with_motif = add_motif_column(
        df=df,
        ref_fasta_path=args.fasta,
        n=args.upstream,
        m=args.downstream,
        motif_col_name=args.motif_col,
        chrom_col=args.chrom_col,
        pos_col=args.pos_col,
        strand_col=args.strand_col,
    )

    if args.inplace:
        output_path = input_path
    else:
        if args.output:
            output_path = Path(args.output)
        else:
            output_path = input_path.with_suffix(input_path.suffix + '.motif')

    df_with_motif.to_csv(output_path, sep=args.sep, index=False)
    print(f'[OK] Motif file written to: {output_path}')


if __name__ == '__main__':
    main()
