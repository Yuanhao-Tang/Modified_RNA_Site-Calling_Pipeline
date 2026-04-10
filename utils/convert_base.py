#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path


def convert_base(in_fa: str, out_fa: str, x: str, y: str) -> None:
    """
    Replace all occurrences of base `x` with base `y` in a reference FASTA file.

    Args:
        in_fa: Input FASTA path.
        out_fa: Output FASTA path.
        x: Source base (A/T/C/G).
        y: Target base (A/T/C/G).
    """
    x = x.upper()
    y = y.upper()

    if len(x) != 1 or len(y) != 1:
        raise ValueError('X and Y must each be a single base character.')

    in_path = Path(in_fa)
    if not in_path.is_file():
        raise FileNotFoundError(f'Input file not found: {in_fa}')

    with in_path.open('r') as fin, Path(out_fa).open('w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(line)
            else:
                fout.write(line.upper().replace(x, y))


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Replace a specified base with another base across a reference FASTA file.'
    )
    parser.add_argument('in_fa', help='Input FASTA path.')
    parser.add_argument('out_fa', help='Output FASTA path.')
    parser.add_argument('x', help='Source base (A/T/C/G).')
    parser.add_argument('y', help='Target base (A/T/C/G).')
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()
    convert_base(args.in_fa, args.out_fa, args.x, args.y)


if __name__ == '__main__':
    main()
# python utils/convert_base.py input.fa output.fa A G
