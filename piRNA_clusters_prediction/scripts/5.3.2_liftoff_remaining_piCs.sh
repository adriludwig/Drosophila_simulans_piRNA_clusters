#!/bin/bash
# This script performs coordinate conversion (using liftover) of piC regions
# It reads non-lifted BED files for each strain
# converts them into a fake GFF3 format, and runs LiftOff with a custom feature type. 

# Paths
REF_GENOME="/data/adriana/TE_Invasions/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna"
STRAIN_PATH="/data/adriana/TE_Invasions/Nanopore/cactus/cactus_singularity/softmasked_genomes/"

# Strain list
STRAINS="AA12 AA1 AH24 AH33 AH35 Cro18 Hin17 Lps3 Lps5 LSD01 LSD07 MD221 MD233 MSC01 MSC04 PUS102 PUS412 SGA12 SGA20 SGA26 SGA27 SGA34"

# Output folder
OUTDIR="liftoff_remaining"
mkdir -p "$OUTDIR"

# Create a feature types file for LiftOff
FEATURES_FILE="${OUTDIR}/piC_feature_types.txt"
echo "piC" > "$FEATURES_FILE"

# Function to convert BED/TSV to fake GFF3
bed_to_fake_gff3() {
    local bed_file="$1"
    local gff_file="$2"
    awk 'BEGIN{OFS="\t"} {
        id = $4":" $1 "_" $2+1 "_" $3;  # GFF is 1-based, BED is 0-based
        print $1, "fake", "piC", $2+1, $3, ".", ".", ".", "ID=" id
    }' "$bed_file" > "$gff_file"
}
# Loop through each strain and process both files
for strain in $STRAINS; do
    GENOME="${STRAIN_PATH}${strain}_final_assembly.softmasked.cleaned.fasta"
    GENOME_FAI="${STRAIN_PATH}${strain}_final_assembly.softmasked.cleaned.fasta.fai"

    for bed_file in "${strain}_non_lifted.bed"; do
        if [ -f "$bed_file" ]; then
            base_name=$(basename "$bed_file" | cut -d. -f1)
            expanded_bed="${OUTDIR}/${base_name}_expanded.bed"
	        fake_gff="${OUTDIR}/${base_name}.gff3"
            lifted_gff="${OUTDIR}/${base_name}_lifted.gff3"
            unmapped_gff="${OUTDIR}/${base_name}_unmapped.gff"
            lifted_bed="${OUTDIR}/${base_name}_lifted.bed"

            echo "Processing $bed_file for strain $strain..."
            bedtools slop -i "$bed_file" -g "$GENOME_FAI" -b 0 > "$expanded_bed"

            # Convert to GFF3
            bed_to_fake_gff3 "$expanded_bed" "$fake_gff"

            # Run liftoff, specifying the feature type file
            liftoff "$REF_GENOME" "$GENOME" \
                -g "$fake_gff" \
                -f "$FEATURES_FILE" -copies \
                -o "$lifted_gff" \
                -u "$unmapped_gff" \
                -p  24

            # Convert lifted GFF3 to BED if successful
            if [ -f "$lifted_gff" ]; then
                awk '$3 == "piC" {print $1"\t"$4-1"\t"$5"\t"$9}' "$lifted_gff" > "$lifted_bed"
                echo "LiftOver completed: $lifted_bed"
            else
                echo "LiftOff failed for $bed_file – GFF not generated."
            fi
        else
            echo "File $bed_file not found for strain $strain."
        fi
    done
done

