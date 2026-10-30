#!/bin/bash

#Parsing Protrac Results. 
#It extracts cluster coordinates and piRNA-related metrics from the clusters.gtf files generated 
#by proTRAC and calculates cluster size and RPKM values for each predicted cluster
# Paths correspond to the internal computing environment and should be adapted if the script is reused.


# Protrac output directories
parent_dirs=("protrac")
# SAM directory
sam_dir="./protrac/"

# Header for output table
header="Chromosome\tSource\tFeature\tStart\tEnd\tScore\tStrand\tFrame\tpiRNA_cluster_no\tdirectionality\tRPM\tfraction_1U_reads\tfraction_10A_reads\tfraction_piRNA_sized_reads\tCluster_size\tRPKM"

# Mapped reads cache
declare -A total_mapped_reads_cache

# Loop through protrac directories
for parent_dir in "${parent_dirs[@]}"; do
    full_parent_dir="${parent_dir}"
    echo "Processing: $full_parent_dir"

    # Loop through sample outputs
    for folder in "$full_parent_dir"/proTRAC_*; do
        if [[ -d "$folder" ]]; then
            sample_name=$(basename "$folder" | cut -d'_' -f2)

            input_gtf="${folder}/clusters.gtf"
            input_sam="${sam_dir}/${sample_name}_Protrac.sam"
            output_file="${full_parent_dir}/${sample_name}_parsed_clusters_simple.txt"


            # Add header to the output file
            echo -e "$header" > "$output_file"

            # Process clusters.gtf
            while IFS=$'\t' read -r chrom source feature start end score strand frame attributes; do
                    [[ $chrom =~ ^# ]] && continue  # Skip commented lines

                    # Extract attributes
                    piRNA_cluster_no=$(echo "$attributes" | grep -oP "piRNA cluster no: \K[0-9]+")
                    directionality=$(echo "$attributes" | grep -oP "directionality: \K[^;]+")
                    mapped_reads=$(echo "$attributes" | grep -oP "mapped sequence reads: \K[0-9\.]+" || echo "0") # These are normalized mapped read counts reported by proTRAC, such as CPM/RPM
                    fraction_1U=$(echo "$attributes" | grep -oP "fraction 1U reads: \K[0-9\.]+" || echo "0")
                    fraction_10A=$(echo "$attributes" | grep -oP "fraction 10A reads: \K[0-9\.]+" || echo "0")
                    fraction_piRNA=$(echo "$attributes" | grep -oP "fraction piRNA-sized reads: \K[0-9\.]+" || echo "0")

                    # Cluster size calculations
                    cluster_size=$((end - start + 1))
                    cluster_size_kb=$(echo "scale=6; $cluster_size /1000" | bc)
                    rpkm=$(echo "scale=6; $mapped_reads/$cluster_size_kb" | bc)

                    # Write to output file
                    echo -e "${chrom}\t${source}\t${feature}\t${start}\t${end}\t${score}\t${strand}\t${frame}\t${piRNA_cluster_no}\t${directionality}\t${mapped_reads}\t${fraction_1U}\t${fraction_10A}\t${fraction_piRNA}\t${cluster_size}\t${rpkm}" >> "$output_file"

                done < "$input_gtf"

                echo "Processed $input_gtf -> $output_file"
      
        fi
    done
done

echo "All samples were processed."

