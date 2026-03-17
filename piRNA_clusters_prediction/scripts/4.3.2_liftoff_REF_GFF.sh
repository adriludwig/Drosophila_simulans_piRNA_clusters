#!/bin/bash
# Paths correspond to the internal computing environment and should be adapted if the scripts are reused.

# Define reference genome and GFF
ref_genome="/data/adriana/TE_Invasions/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna"
ref_gff="./GCF_filtered_no_lncRNA_noTE.gff"

# Folder with target assemblies
assembly_folder="/data/adriana/TE_Invasions/Nanopore/eliminate_contamin_last/final_assemblies/renamed"

# Loop over each assembly
for assembly in "$assembly_folder"/*final_assembly.fasta; do
    # Get base name without path and extension
    base_name=$(basename "$assembly" .fasta)

    # Define output files
    output_gff="./${base_name}_final_assembly_liftoff_noTE_nolncRNA.gff"
    unmapped_txt="./${base_name}_liftoff_noTE_nolncRNA_unmapped.txt"

    # Run Liftoff
    liftoff -g "$ref_gff" -o "$output_gff" -u "$unmapped_txt" \
            "$assembly" "$ref_genome"

    echo "Liftoff completed for $base_name → $output_gff"
done
