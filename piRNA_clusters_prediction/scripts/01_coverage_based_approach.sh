#!/bin/bash
# This script is a first part to run coverage-based approach to identify piCs
# Paths correspond to the internal computing environment and should be adapted if the script is reused.

source /data/adriana/miniconda3/etc/profile.d/conda.sh

mapped_dir="/data/adriana/TE_Invasions/piRNA/reads/TEs_restrictive_plus_antisense_genes/only_antisense/coverage/"
reads_dir="/data/adriana/TE_Invasions/piRNA/mapping_reads_genes/filtered_reads_no_sense/final_fastq/"

# set genome path for each strain and create bowtie index
for file in ${reads_dir}*_dedup.fastq; do
	sample=$(basename "$file" _dedup.fastq)
        strain=$(echo "$sample" | cut -d'-' -f1)

	conda activate bioinfo
	
	#creating bowtie index
    echo "${strain}"
    genomes_dir="/data/adriana/TE_Invasions/Nanopore/eliminate_contamin_last/final_assemblies/"
	sample_genome2="/data/adriana/TE_Invasions/Nanopore/eliminate_contamin_last/final_assemblies/renamed/${strain}_final_assembly.fasta"
	index_prefix="${genomes_dir}${strain}_bowtie_index" 
    genome_path=$(echo $sample_genome2)
	if [[ -e "${index_prefix}.1.ebwt" ]]; then
                echo "Bowtie index already exists for $genome_path"
    else
                echo "Building Bowtie index for $genome_path..."
                bowtie-build "$genome_path" "$index_prefix"
    fi      
    echo  "sample_genome: $sample_genome"
	ls $sample_genome
        

######################################### Unique Mapping #################################################################
    #To maximise the number of unique mappings, we performed two rounds of mapping using Bowtie v1.3.1. In the first round, no mismatches 
     #were allowed. In the second round, the unmapped reads from the first round were remapped, allowing for one mismatch.
    #Finally, we used samtools merge to combine the SAM files. 	
	bowtie -p 24 -n 0 -l 12 -a -m 1 -q "${file}" -x  "$index_prefix" -S "${mapped_dir}${sample}_perfect_unique2.sam"  
    samtools view -f 4 "${mapped_dir}${sample}_perfect_unique2.sam" | samtools fastq -@ 24 - > "${mapped_dir}${sample}_unmapped_from_perfect_unique2.fastq"
    bowtie -p 24 -n 1 -l 12 -a -m 1 -q  "${mapped_dir}${sample}_unmapped_from_perfect_unique2.fastq" -x  "$index_prefix" -S "${mapped_dir}${sample}_one_mismatch_unique2.sam"
    samtools merge -f "${mapped_dir}${sample}_combined_unique2.sam" "${mapped_dir}${sample}_perfect_unique2.sam" "${mapped_dir}${sample}_one_mismatch_unique2.sam"
    samtools sort -@ 24  -o "${mapped_dir}${sample}_combined_unique_sorted2.bam" "${mapped_dir}${sample}_combined_unique2.sam"

		
########################################## estimating coverage - filter 8CPM #################################################
	#Windows of 500 bp were created with bamCoverage to calculate normalised counts (CPM). 
	#Windows with an expression of at least 8 CPM were filtered and merged if they were located within 10 kb
	#of each other using samtools. 
	##Following in part Srivastan et al 2024 protocol
   	bam_file="${mapped_dir}${sample}_combined_unique_sorted2.bam"
    bedgraph_file="${mapped_dir}${sample}_500bp.bedGraph"
    filtered_bedgraph="${mapped_dir}${sample}_500bp_8CPM.bed"
    merged_bed="${mapped_dir}${sample}_500bp_8CPM_merge10kb.bed"
    samtools index "${bam_file}"
    
	#calculate coverage
	bamCoverage -bs 500 -b "$bam_file" -of bedgraph -p 2 -o "$bedgraph_file" --normalizeUsing CPM
    awk -v OFS='\t' '{ if ($4 > 8) print $0 }' "$bedgraph_file" > "$filtered_bedgraph"
    
    #merge windows within 10kb into domains
	sort -k1,1 -k2,2n "$filtered_bedgraph" | bedtools merge -d 10000 -c 4 -o sum > "$merged_bed"
	

########################################## Remapping to estimate % 1U ######################################################################
	#To analyse the 1U bias, the genomic regions were extracted using bedtools getfasta, and the reads were remapped to these regions with the same parameters.
	#The BAM files were subsequently divided into separate files and we used FastQC to produce a quality control report for each of them
	#from where the 1U proportion of each region was extracted. 
        
    #extract the small RNA-enriched domains from the respective genome files 
    piRNA_domains="${mapped_dir}${sample}_piRNA_genome.uniq.8CPM_merge10kb_domain.fa"
    bedtools getfasta -fi "$sample_genome2" -bed "$merged_bed" > "$piRNA_domains"

    #Remapping Reads with Bowtie 1
    domain_index="${mapped_dir}${sample}_piRNA_genome.uniq.8CPM_merge10kb_domain_bowtie_index"
    remap_sam1="${mapped_dir}${sample}_remap_perfect_unique_8CPM.sam"
    remap_sam2="${mapped_dir}${sample}_remap_one_mismatch_unique_8CPM.sam"
    remap_bam="${mapped_dir}${sample}_remap_unique_8CPM.bam"
        
    bowtie-build "$piRNA_domains" "$domain_index"
    bowtie -p 24 -n 0 -l 12 -a -m 1 -q -x "$domain_index" -q "${file}"  -S "$remap_sam1"
    samtools view -f 4 "$remap_sam1"| samtools fastq -@ 24 - > "${mapped_dir}${sample}_unremapped_from_perfect_unique_8CPM.fastq"
    bowtie -p 24 -n 1 -l 14 -a -m 1 -x "$domain_index" -q "${mapped_dir}${sample}_unremapped_from_perfect_unique_8CPM.fastq"  -S "$remap_sam2"
    samtools merge -f "${mapped_dir}${sample}_remapped_combined_unique_8CPM.sam" "$remap_sam1"  "$remap_sam2"
    samtools view -bS -@ 24  "${mapped_dir}${sample}_remapped_combined_unique_8CPM.sam" > "$remap_bam"

	#split this bam file by each merged domain and obtain the 1U %
    bamtools split -in $remap_bam -reference 
	mkdir -p "${mapped_dir}bam_split"
	mv ${mapped_dir}$sample*REF*.bam ${mapped_dir}bam_split
done

#Obtaining FastQC for the potential piRNA clusters
split_dir="${mapped_dir}bam_split/"
output_dir="${split_dir}fastqc_files2/"
mkdir -p "${output_dir}"
find "${split_dir}" -maxdepth 1 -name "*.bam" | parallel -j 46 fastqc {} -t 3 -q --extract -f bam_mapped -o "${output_dir}"

fastqc_dir="${split_dir}/fastqc_files2/"

# Moving fastqc files
for bam_file in "$split_dir"/*.bam; do
    base_name=$(basename "$bam_file" .bam)
    fastqc_data_file="${fastqc_dir}/${base_name}_fastqc/fastqc_data.txt"
	fastqc_data_file_mv="${split_dir}/${base_name}_fastqc_data.txt"

 	echo "Checking: $fastqc_data_file"

    if [ -f "$fastqc_data_file" ]; then
    	mv "$fastqc_data_file" "$fastqc_data_file_mv"
     	echo "Moved: $fastqc_data_file → $fastqc_data_file_mv"
    else
   		echo "FastQC data file NOT found for: $bam_file"
 fi
done

#Extracting 1U% info
# Find bam split files
find "${split_dir}" -maxdepth 1  -name "*REF*.bam" | while read bam_file; do
base_name=$(basename "$bam_file" .bam)

# Find the corresponding fastqc_data.txt file using the base name
fastqc_file="${bam_file%/*}/${base_name}_fastqc_data.txt"

# Check if the corresponding fastqc_data.txt file exists
if [ -f "$fastqc_file" ] && [ ! -f "${fastqc_file%.fastqc_data.txt}.1U.txt" ]; then
	echo "processing: ${bam_file}"
	# Process the fastqc_data.txt 
	grep -A 1 "^#Base[[:blank:]]*G" "$fastqc_file" | \
	tail -n 1 | \
	awk -v fname=$(basename "$fastqc_file") '{print $0 "\t" fname}' | \
	cut -f4,6 | \
	sed 's/_remap_unique_8CPM.REF_/\t/' | \
	sed 's/:/\t/' | \
	sed 's/_fastqc_data.txt//' | \
	awk -F '\t' -v OFS='\t' '{split($4, range, "-"); print $1, $2, $3, range[1], range[2]}' > "${fastqc_file%.fastqc_data.txt}.1U.txt"

else
    if [ ! -f "$fastqc_file" ]; then
        echo "No matching fastqc_data.txt file for ${bam_file}"
  	elif [ -f "${fastqc_file%.fastqc_data.txt}.1U.txt" ]; then
      echo "Output file already exists. Skipping processing."
 	fi
fi
done

#concatenate files from all strains into a single file
find "${split_dir}" -maxdepth 1 -type f -name '*1U.txt' -exec cat {} + > all_samples_1U_coverage_based.txt




