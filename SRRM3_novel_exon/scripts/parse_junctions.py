#!/usr/bin/env python3
import sys
import re

def parse_cigar(cigar, pos):
    """Parse CIGAR string and extract splice junctions"""
    junctions = []
    current_pos = pos

    # Parse CIGAR operations
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')
    for match in pattern.finditer(cigar):
        length = int(match.group(1))
        op = match.group(2)

        if op in ('M', 'D', '=', 'X'):  # Consume reference
            current_pos += length
        elif op == 'N':  # Splice junction
            # Junction: from (current_pos - 1) to (current_pos + length)
            # In 1-based coords: junction spans from current_pos to current_pos + length - 1
            junc_start = current_pos  # Last base of upstream exon (1-based)
            junc_end = current_pos + length  # First base of downstream exon (1-based)
            junctions.append((junc_start, junc_end))
            current_pos += length
        elif op in ('I', 'S', 'H', 'P'):  # Don't consume reference
            pass

    return junctions

# Read SAM input from stdin
for line in sys.stdin:
    if line.startswith('@'):
        continue
    fields = line.strip().split('\t')
    if len(fields) < 6:
        continue

    chrom = fields[2]
    pos = int(fields[3])  # 1-based position
    cigar = fields[5]

    if 'N' not in cigar:
        continue

    junctions = parse_cigar(cigar, pos)
    for junc_start, junc_end in junctions:
        print(f"{chrom}\t{junc_start}\t{junc_end}")
