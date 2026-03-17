#!/usr/bin/env python3

import sys
from collections import defaultdict

input_file = sys.argv[1]
output_file = sys.argv[2]

merged = defaultdict(list)

with open(input_file) as infile:
    for line in infile:
        chrom, start, end, pic_id = line.strip().split()
        merged[pic_id].append((chrom, int(start), int(end)))

with open(output_file, "w") as out:
    for pic_id, intervals in merged.items():
        chrom = intervals[0][0]
        start = min(start for _, start, _ in intervals)
        end   = max(end for _, _, end in intervals)
        out.write(f"{chrom}\t{start}\t{end}\t{pic_id}\n")
