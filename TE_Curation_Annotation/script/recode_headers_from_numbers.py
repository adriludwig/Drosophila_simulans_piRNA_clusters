import sys

def recode_fasta_headers(original_fasta, masked_fasta, output_fasta):
    """Recode masked FASTA headers back to original names based on the original FASTA."""
    # Create a mapping of new IDs to original headers
    original_headers = {}
    with open(original_fasta, 'r') as infile:
        count = 1
        for line in infile:
            if line.startswith('>'):
                header = line.strip()
                original_headers[f'>{count}'] = header  # Map new ID to original header
                count += 1

    # Recode headers in the masked output
    with open(masked_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Replace new ID with original header
                new_header = line.strip()  # Get the new header (e.g., >1)
                original_header = original_headers.get(new_header, new_header)  # Find original header
                outfile.write(f'{original_header}\n')
            else:
                outfile.write(line)  # Write sequence lines as-is

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python recode_headers_from_numbers.py <original_fasta> <masked_fasta> <output_fasta>")
        sys.exit(1)

    original_fasta = sys.argv[1]
    masked_fasta = sys.argv[2]
    output_fasta = sys.argv[3]
    recode_fasta_headers(original_fasta, masked_fasta, output_fasta)
