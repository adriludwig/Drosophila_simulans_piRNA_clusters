                                                                                                             
#!/bin/bash
#pwd: /data/adriana/TE_Invasions/Nanopore/cactus/cactus_singularity/halLiftover/non_lifted/lift_remaining_flanking/concat_lifted

mkdir -p non_lifted

# Loop through all *_merged.bed files
for lifted_file in *_lifted.bed; do
    # Extract the sample name 
    sample=$(basename "$lifted_file" _lifted.bed)

    # Define corresponding original file
    original_file="/data/adriana/TE_Invasions/Nanopore/cactus/cactus_singularity/piCs_bed/${sample}.simplified.bed"

    # Check if original file exists
    if [[ ! -f "$original_file" ]]; then
        echo "Original file not found for $sample, skipping."
        continue
    fi

    # Extract piC IDs from each file
    cut -f4 "$lifted_file" | sort > lifted_ids.txt
    cut -f4 "$original_file" | sort > original_ids.txt

    # Get IDs in original but not in lifted
    comm -23 original_ids.txt lifted_ids.txt > not_lifted_ids.txt

    # Extract full lines for non-lifted IDs from the original file
    grep -Fwf not_lifted_ids.txt "$original_file" > "non_lifted/${sample}_non_lifted.bed"

    echo "Processed $sample: saved non-lifted piCs to non_lifted/${sample}_non_lifted.bed"
done

# Clean up temporary files
rm lifted_ids.txt original_ids.txt not_lifted_ids.txt
