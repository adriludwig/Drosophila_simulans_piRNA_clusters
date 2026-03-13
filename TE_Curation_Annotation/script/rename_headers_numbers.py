import sys

def rename_fasta_headers(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        count = 1
        for line in infile:
            if line.startswith('>'):
                outfile.write(f'>{count}\n')
                count += 1
            else:
                outfile.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python rename_headers_numbers.py <input_fasta> <output_fasta>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    rename_fasta_headers(input_fasta, output_fasta)
