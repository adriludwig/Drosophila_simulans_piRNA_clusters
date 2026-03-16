#!/bin/bash
# Paths correspond to the internal computing environment and should be adapted if the script is reused.

source /data/adriana/miniconda3/etc/profile.d/conda.sh
mapped_dir="/data/adriana/TE_Invasions/piRNA/reads/TEs_restrictive_plus_antisense_genes/only_antisense/protrac/"
reads_dir="/data/adriana/TE_Invasions/piRNA/mapping_reads_genes/filtered_reads_no_sense/final_fastq/"

##############################Mapping - allowing multi-mapping################################################
for file in ${reads_dir}*_dedup.fastq; do
        sample=$(basename "$file" _dedup.fastq)
        strain=$(echo "$sample" | cut -d'-' -f1)

        #Mapping reads against the genomes
        echo "${strain}"
        sample_genome="/data/adriana/TE_Invasions/Nanopore/eliminate_contamin_last/final_assemblies/renamed/${strain}_final_assembly.fasta"
        genome_path=$(echo $sample_genome)
        index_prefix="${genome_path}${strain}_bowtie_index" 
        conda activate bioinfo
        if [[ -e "${index_prefix}.1.ebwt" ]]; then
                echo "Bowtie index already exists for $genome_path"
        else
                echo "Building Bowtie index for $genome_path..."
                bowtie-build "$genome_path" "$index_prefix"
        fi      
        echo  "sample_genome: $sample_genome"
        ls $sample_genome
        


##############################PROTRAC#########################################
	    te_annot_dir="/data/adriana/TE_Invasions/piRNA/reads/data_nov24/TE_libraries/"
        RM="${te_annot_dir}RM_${strain}"
        mapped_genome="${mapped_dir}${sample}_Protrac.sam"
        bowtie -n 1 -l 12 -a --best --strata -y --chunkmbs 1024 -S "$index_prefix" "${file}"  >"$mapped_genome"
                
        perl ~/bin/proTRAC_2.4.3.pl -swsize 500 -swincr 100 -clsize 500 -1Tor10A 0.6 -pimin 23 -clhitsn 10 -pdens 0.05 -pti -map "${mapped_genome}"  -genome $sample_genome -repeatmasker $RM &  # Run in background

    	# Limit concurrent jobs to avoid CPU overload
    	if [[ $(jobs -r -p | wc -l) -ge 10 ]]; then
        	wait -n  # Wait for one job to finish before starting another
   		fi
done

