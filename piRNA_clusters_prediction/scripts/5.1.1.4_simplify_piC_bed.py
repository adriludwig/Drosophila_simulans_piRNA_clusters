#!/usr/bin/env python3
# This script is to create a simple 4th column for piC bed files to be able to get the
# coordinates transfered. A reference file is also created to recover the original information

import os

# Process all *cleaned.bed files in the current directory
bed_files = glob.glob("*cleaned.bed")

for bed_file in bed_files:
    output_bed = bed_file.replace("cleaned.bed", "simplified.bed")
    reference_file = bed_file.replace("cleaned.bed", "reference_with_piC.txt")

    with open(bed_file, 'r') as infile, \
         open(output_bed, 'w') as outfile, \
         open(reference_file, 'w') as refout:

        counter = 1
        #outfile.write("chrom\tstart\tend\tpiC_ID\n")
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue  # Skip invalid lines
            chrom, start, end = parts[0], parts[1], parts[2]
            identifier = f"piC{counter}"
            outfile.write(f"{chrom}\t{start}\t{end}\t{identifier}\n")
            refout.write(line.strip() + f"\t{identifier}\n")
            counter += 1

    print(f"Simplified BED created: {output_bed}")
    print(f"Reference file created: {reference_file}")