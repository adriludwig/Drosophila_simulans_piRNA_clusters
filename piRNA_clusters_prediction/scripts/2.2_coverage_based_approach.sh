#!/bin/bash
# This script is the second part to run coverage-based approach to identify piCs
# Paths correspond to the internal computing environment and should be adapted if the script is reused.

for bam in *_combined_unique_sorted2.bam; do
    sample=$(basename "$bam" _combined_unique_sorted2.bam)
    bed="${sample}_500bp_8CPM_merge10kb.bed"

    echo "Processing BAM file: $bam"
    echo "Looking for BED file: $bed"

    sorted_bam="$bam"

	# Define log file
	log_file="${sample}_processing.log"

	#Getting coverage, RPM and RPKM of potential piC
	if [[ -f "$bed" ]]; then
    	# Start logging
    	echo "Processing sample: $sample" | tee "$log_file"
    	# get total mapped reads
    	total_reads=$(samtools view -F 4 "$sorted_bam" | wc -l)

    	echo "Total reads in BAM file: $total_reads" | tee -a "$log_file"

    	output_file="${sample}_RPM_RPKM_results_correct_CPM8.bed"

    	# Add headers to the output file
    	echo -e "Chromosome\tStart\tEnd\tsum_CPM_500bp\tCoverage\tCluster_size\tRPKM\tRPM\tTotal_Mapped_Reads" > "$output_file"
		#getting coverage, RPM and RPKM
    	while read -r chrom start end name; do
        coverage=$(samtools view -c "$sorted_bam" "$chrom:$start-$end")
        region_length=$((end - start))
        rpm=$(echo "scale=6; ($coverage / $total_reads) * 1000000" | bc)
		rpkm=$(echo "scale=6; ($coverage * 1000000000) / ($region_length * $total_reads)" | bc)

        echo -e "$chrom\t$start\t$end\t$name\t$coverage\t$region_length\t$rpkm\t$rpm\t$total_reads" >> "$output_file"
done < "$bed"

    echo "Finished processing $sample" | tee -a "$log_file"
else
    echo "Missing BED file for sample: $sample" | tee -a "$log_file" >&2
fi
done

echo -e "chr\tstart\tend\tsum_CPM_500bp\tCoverage\tCluster_size\tRPKM\tRPM\tTotal_Map_Reads\tsample" > all_samples_CPM8_coverage_based.bed

for file in *_RPM_RPKM_results_correct_CPM8.bed; do
    sample=${file%_RPM_RPKM_results_correct_CPM8.bed}
    awk -v s="$sample" 'NR > 1 {print $0 "\t" s}' "$file" >> all_samples_CPM8_coverage_based.bed
done
