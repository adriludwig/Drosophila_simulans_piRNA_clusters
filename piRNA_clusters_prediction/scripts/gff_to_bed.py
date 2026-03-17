import sys

# Function to convert GFF to BED
def gff_to_bed(gff_file, bed_file):
    with open(gff_file, 'r') as gff, open(bed_file, 'w') as bed:
        for line in gff:
            if line.startswith("#"):  # Skip header lines
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            
            chrom = fields[0]  # Chromosome
            start = int(fields[3]) - 1  # BED is 0-based, GFF is 1-based
            end = int(fields[4])  # End position
            name = fields[8].split(";")[0] if "=" in fields[8] else "."  # Feature name
            score = fields[5] if fields[5] != "." else "0"  # Score
            strand = fields[6]  # Strand
            
            bed.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python gff_to_bed.py input.gff output.bed")
        sys.exit(1)
    
    gff_file = sys.argv[1]
    bed_file = sys.argv[2]
    
    gff_to_bed(gff_file, bed_file)
    print(f"Conversion completed: {bed_file}")
