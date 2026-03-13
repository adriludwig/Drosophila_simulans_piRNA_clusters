#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

# TE annotation pipeline
#
# Performs a two-step annotation of TEs in genome assemblies.
# First, conserved TE families are annotated using a curated D. simulans TE library.
# Confirmed TE fragments are then masked and the genome is re-annotated using a
# broader Drosophila TE library to detect more divergent TE fragments.
#
# Dependencies (All tools must be available in the system PATH):
#   RepeatMasker
#   BLAST+
#   BEDtools
#   Python3
#
#Additional Python helper scripts required by this workflow:
#- rename_headers_numbers.py
#- recode_headers_from_numbers.py

initial_dir=$(pwd)

# Parsing command-line arguments
while getopts "f:p:g:t:c:h" opt; do
    case $opt in
        f) SAMPLE_LIST_FILE="$OPTARG" ;;
        p) path_genome="$OPTARG" ;;
        g) GENOME_SUFFIX="$OPTARG" ;;
        t) specific_TE_library="$OPTARG" ;;
        c) general_TE_lib="$OPTARG" ;;
        h) echo "Usage: $0 -f SAMPLE_LIST_FILE -p PATH_GENOME -g GENOME_SUFFIX -t SPECIFIC_TE_LIBRARY -c GENERAL_TE_LIBRARY"
           exit 0 ;;
        *) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
    esac
done

# Check that required arguments are provided
if [[ -z "$SAMPLE_LIST_FILE" || -z "$path_genome" || -z "$GENOME_SUFFIX"|| -z "$specific_TE_library" || -z "$general_TE_lib" ]]; then
    echo "Missing required argument(s)."
    echo "Usage: $0 -f SAMPLE_LIST_FILE -p PATH_GENOME -g GENOME_SUFFIX -t SPECIFIC_TE_LIBRARY -c GENERAL_TE_LIBRARY"
    exit 1
fi

#Check if the BLAST database exists for specific_TE_library
if [ ! -f "${specific_TE_library}.nin" ]; then
   echo "Creating BLAST database for specific TE library: $specific_TE_library"
  makeblastdb -in "$specific_TE_library" -dbtype nucl -out "$specific_TE_library" -title "$specific_TE_library"
else
   echo "BLAST database already exists for specific TE library: $specific_TE_library"
fi

# Check if the BLAST database exists for general_TE_lib
if [ ! -f "${general_TE_lib}.nin" ]; then
    echo "Creating BLAST database for general TE library: $general_TE_lib"
    makeblastdb -in "$general_TE_lib" -dbtype nucl -out "$general_TE_lib" -title "$general_TE_lib"
else
    echo "BLAST database already exists for general TE library: $general_TE_lib"
fi

# Read samples from the list file and process each one
while IFS= read -r sample || [[ -n "$sample" ]]; do
    cd "$initial_dir" 
    echo "Processing sample: $sample"
    
    # Set paths based on existing variables
    genome="${path_genome}/${sample}${GENOME_SUFFIX}"
    
    # Ensure the genome file exists
    if [[ ! -f "$genome" ]]; then
        echo "Genome file ${genome} does not exist. Skipping RepeatMasker for ${sample}."
        continue
    fi
###############################################################################
# Module 1: Annotation of conserved TE families using the curated D. simulans TE library
###############################################################################

    RM_dir="RM_${sample}_nolow/"
    RM_file="${initial_dir}/${RM_dir}${sample}${GENOME_SUFFIX}.out"
    BED_file="${sample}_genomic.RM.bed"
    RM_fasta="${sample}_genomic.RM.fasta"
    if [ -f "$genome" ]; then
	RepeatMasker -s -norna -div 20 -nolow -gff -no_is -pa 48 -a -lib "${specific_TE_library}" "${genome}" -dir "${RM_dir}"
	else
        echo "Genome file ${genome} does not exist. Skipping RepeatMasker for ${sample}."
        continue
	fi

     #accessing results fromRepeatMasker 
    cd "$RM_dir" || { echo "Failed to cd into $RM_dir"; exit 1; }

    #rmsk to bed; rmsk 1-based; bed 0-based
	if [ -f "${RM_file}" ]; then	
		awk 'OFS= "\t" {print $5,$6-1,$7,$10}' "${RM_file}" |tail -n +4 > "${BED_file}"
 	else
       echo "RepeatMasker output file ${RM_file} not found. Skipping bed creation for ${sample}."
      continue
    fi

	#retrieving all annotated TE sequences
    if [ -f "${BED_file}" ]; then
       bedtools getfasta -fi "$genome" -bed "${BED_file}" -fo "${RM_fasta}" -name
    else
       echo "BED file ${BED_file} does not exist. Skipping bedtools getfasta for ${sample}."
      continue
    fi
	
	#renaming headers with numbers
	if [ -f "${RM_fasta}" ]; then
		python rename_headers_numbers.py "${RM_fasta}" "${RM_fasta}_renamed"
	else
       echo "RM fasta file ${RM_fasta} not found. Skipping renaming for ${sample}."
      continue
 fi
	#running repeatmasker with -noint to mask only simple repeats and low complexity sequences
	if [ -f "${RM_fasta}_renamed" ]; then
		RepeatMasker -s -norna -gff -no_is -noint -pa 48 "${RM_fasta}_renamed"
	else
       echo "Renamed fasta file ${RM_fasta}_renamed not found. Skipping RepeatMasker for ${sample}."
      continue
 fi

	#renaming the headers back to original names
	if [ -f "${RM_fasta}_renamed.masked" ]; then
	python recode_headers_from_numbers.py "${RM_fasta}" "${RM_fasta}_renamed.masked" "${RM_fasta}_correct_name"
	else
       echo "Masked file ${RM_fasta}_renamed.masked not found. Skipping header recoding for ${sample}."
      continue
  fi

	#blast of annotated TEs against the Manual curated Dsim library to eliminate non-TE sequences eventually identified by RepeatMasker
	if [ -f "${RM_fasta}_correct_name" ]; then
	blastn -query  "${RM_fasta}_correct_name" -db "${specific_TE_library}" -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -num_threads 48 -word_size 7 -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"  -reward 2  -penalty -3  -gapopen 5  -gapextend 2 -out ${sample}_RM_blast
	else
        echo "Correct name fasta file ${RM_fasta}_correct_name not found. Skipping blastn for ${sample}."
        continue
   	fi

	#processing blast output and saving as bed file; first position from previous bed
	if [ -f "${sample}_RM_blast" ]; then
	awk -F "::" 'OFS="\t"{print $1, $2}' "${sample}_RM_blast" |awk -F ":" 'OFS="\t" {print $1, $2}'|awk 'OFS="\t" {split($3, a, "-"); $3 = a[1] "\t" a[2]; print $0}' >"${sample}_RM_processed_blast"
	awk -F'\t' '{split($6, a, "#"); print $2, $3, $4, a[1]}' OFS='\t' "${sample}_RM_processed_blast" >"${sample}_RM_processed_blast.bed"
	else
       echo "Blast output ${sample}_RM_blast not found. Skipping further processing for ${sample}."
      continue
	fi
    	
	#separate TEs, sort, merge, concatenate back to a single file
       if [[ -f "${sample}_RM_processed_blast.bed" ]]; then
       sort -k1,1 -k2,2n -k4,4 "${sample}_RM_processed_blast.bed" | awk '{print > $4".tmp.bed"}' 
      for file in *.tmp.bed; do bedtools merge -d 50  -i "$file" > "merged_$file"; done
        cat merged_*.tmp.bed|bedtools sort -i - > "${sample}_merged_output.bed"
       rm *.tmp.bed
      else
  	echo "File ${sample}_RM_processed_blast.bed does not exist. Skipping sort and split step for ${sample}."
   	continue
	fi
	
	#filtering fragments with minimum of 80 bp
       if [[ -f "${sample}_merged_output.bed" ]]; then
	awk '$3-$2>79' "${sample}_merged_output.bed" > "${sample}_merged_output.min80.bed"
       else
	echo "File ${sample}_merged_output.bed does not exist. Skipping fragment filtering for ${sample}."
   	continue
	fi

##	#recovering merged fragments
      if [[ -f "${sample}_merged_output.min80.bed" && -f "${genome}" ]]; then
	bedtools getfasta -fi "$genome" -bed "${sample}_merged_output.min80.bed"  -fo "${sample}_RM_BLAST_MERGED.fasta" 
       else
	echo "File ${sample}_merged_output.min80.bed or ${genome} does not exist. Skipping fasta recovery for ${sample}."
 	continue
	fi
	
	#blast to recover the classification of the confirmed TE fragments
	if [[ -f "${sample}_RM_BLAST_MERGED.fasta" && -f "${specific_TE_library}" ]]; then
       blastn -query  "${sample}_RM_BLAST_MERGED.fasta" -db "${specific_TE_library}" -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -num_threads 48 -word_size 7 -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart  send evalue bitscore" -reward 2  -penalty -3  -gapopen 5  -gapextend 2 -out ${sample}_RM_Dsim_blast_filtered  
      else
	echo "File ${sample}_RM_BLAST_MERGED.fasta or ${specific_TE_library} does not exist. Skipping BLAST step for ${sample}."
 	continue
	fi
	
	#blast to bed; first position is already 0-based
       if [[ -f "${sample}_RM_Dsim_blast_filtered" ]]; then
	awk 'BEGIN {OFS="\t"} {split($1, a, "[:\\-]"); print a[1], a[2], a[3], $2, $3, $4, $5, $6, $7}' "${sample}_RM_Dsim_blast_filtered"> "${sample}_RM_Dsim_blast_filtered.bed_complete"
       awk 'OFS="\t" {print $1,$2,$3,$4}' "${sample}_RM_Dsim_blast_filtered.bed_complete" >"${sample}_RM_Dsim_blast_filtered.bed"
      else
	echo "File ${sample}_RM_Dsim_blast_filtered does not exist. Skipping BLAST-to-BED conversion for ${sample}."
 	continue
	fi
	
###############################################################################
# Module 2: Annotation of more divergent TE fragments using a broader
# Drosophila TE library after masking TE fragments identified in Module 1.
###############################################################################	

	#mask genome with the confirmed Dsim TE fragments
	if [[ -f "${sample}_RM_Dsim_blast_filtered.bed" && -f "${genome}" ]]; then
	bedtools maskfasta -fi "${genome}" -bed "${sample}_RM_Dsim_blast_filtered.bed" -fo "${sample}_maskedTEDsim.fasta"
       else
		echo "File ${sample}_RM_Dsim_blast_filtered.bed or ${genome} does not exist. Skipping genome masking for ${sample}."
 	continue
	fi

	#run repeatmasker with Dsim + Dmel + other species library in the masked genome
       if [[ -f "${sample}_maskedTEDsim.fasta" ]]; then
       RepeatMasker -s -norna -gff -no_is -nolow -pa 48 -lib "${general_TE_lib}"  "${sample}_maskedTEDsim.fasta" -dir ./"RM_${sample}_lib_others"
	else
 	echo "File ${sample}_maskedTEDsim.fasta does not exist. Skipping RepeatMasker step for ${sample}."
  	continue
	fi
	
#	#rmsk to bed; rmsk 1-based; bed 0-based
	if [[ -f "RM_${sample}_lib_others/${sample}_maskedTEDsim.fasta.out" ]]; then
		awk 'OFS= "\t" {print $5,$6-1,$7,$10}' "RM_${sample}_lib_others"/"${sample}_maskedTEDsim.fasta.out"  |tail -n +4 > "${sample}_maskedTEDsim.bed" 
	#retrieving all annotated TE sequences
	bedtools getfasta -fi "$genome" -bed "${sample}_maskedTEDsim.bed"  -fo "${sample}_RM_others.fasta"  -name
	else
   	echo "RepeatMasker output file RM_${sample}_lib_others/${sample}_maskedTEDsim.fasta.out does not exist. Skipping TE conversion to BED."
   	continue
	fi

	#renaming headers with numbers
	if [[ -f "${sample}_RM_others.fasta" ]]; then
		python rename_headers_numbers.py "${sample}_RM_others.fasta" "${sample}_RM_others.fasta_renamed"
	else
 	 	echo "Fasta file ${sample}_RM_others.fasta does not exist. Skipping header renaming."
   	continue
	fi

	#running repeatmasker with -noint to mask only simple repeats and low complexity sequences
	if [[ -f "${sample}_RM_others.fasta_renamed" ]]; then
		RepeatMasker -s -norna -gff -no_is -noint -pa 48 "${sample}_RM_others.fasta_renamed"
	else
 		echo "Renamed fasta file ${sample}_RM_others.fasta_renamed does not exist. Skipping RepeatMasker."
   	continue
	fi
	
	#renaming the headers back to original names
	if [[ -f "${sample}_RM_others.fasta" && -f "${sample}_RM_others.fasta_renamed.masked" ]]; then
		python recode_headers_from_numbers.py "${sample}_RM_others.fasta" "${sample}_RM_others.fasta_renamed.masked" "${sample}_RM_others_correct_name.fasta"
	else
	 	 echo "One or both of the files ${sample}_RM_others.fasta or ${sample}_RM_others.fasta_renamed.masked do not exist. Skipping header recoding."
  		 continue
	fi
	
	#blast of annotated TEs against the manual curated library to eliminate non-TE sequences eventually identified by RepeatMasker
	if [[ -f "${sample}_RM_others_correct_name.fasta" && -f "${general_TE_lib}" ]]; then
		blastn -query  "${sample}_RM_others_correct_name.fasta" -db "${general_TE_lib}" -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -num_threads 48 -word_size 7 -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"  -reward 2  -penalty -3  -gapopen 5  -gapextend 2 -out "${sample}_RM_blast_otherTEs"
	else
  	echo "Required files ${sample}_RM_others_correct_name.fasta or ${general_TE_lib} do not exist. Skipping BLAST step."
   	continue
	fi
	
	#processing blast output and saving as bed file; first position from previous bed
	if [[ -f "${sample}_RM_blast_otherTEs" ]]; then
		awk -F "::" 'OFS="\t"{print $1, $2}' "${sample}_RM_blast_otherTEs" |awk -F ":" 'OFS="\t" {print $1, $2}'|awk 'OFS="\t" {split($3, a, "-"); $3 = a[1] "\t" a[2]; print $0}' >"${sample}_RM_processed_blast_otherTEs"
		awk 'OFS="\t" {print $2, $3, $4, $1}' "${sample}_RM_processed_blast_otherTEs" >"${sample}_RM_processed_blast_otherTEs.bed"
		
#		#separate TEs, sort, merge, concatenate back to a single file
		sort -k1,1 -k2,2n -k4,4 "${sample}_RM_processed_blast_otherTEs.bed" | awk '{print > $4".tmp.bed"}' 
		for file in *.tmp.bed; do bedtools merge -d 50  -i "$file" > "merged_$file"; done
		cat merged_*.tmp.bed|bedtools sort -i - > "${sample}_merged_output_otherTEs.bed"
		rm *.tmp.bed
	else
  	echo "BLAST output file ${sample}_RM_blast_otherTEs does not exist. Skipping processing."
   	continue
	fi

	#filtering fragments with minimum of 80 bp of sequence
	if [[ -f "${sample}_merged_output_otherTEs.bed" ]]; then
		awk '$3-$2>79' "${sample}_merged_output_otherTEs.bed" > "${sample}_merged_output_otherTEs.min80.bed"
#		#recovering merged fragments
		bedtools getfasta -fi "$genome" -bed "${sample}_merged_output_otherTEs.min80.bed"  -fo "${sample}_RM_BLAST_MERGED_otherTEs.fasta" 
	else
  	echo "File ${sample}_merged_output_otherTEs.bed does not exist. Skipping filtering."
  	continue
	fi
	
	#blast to recover the classification of the confirmed TE fragments
	if [[ -f "${sample}_RM_BLAST_MERGED_otherTEs.fasta" && -f "${general_TE_lib}" ]]; then
		blastn -query  "${sample}_RM_BLAST_MERGED_otherTEs.fasta" -db "${general_TE_lib}" -evalue 1e-10 -max_target_seqs 1 -max_hsps 1 -num_threads 48 -word_size 7 -outfmt "6 qseqid qlen sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"  -reward 2  -penalty -3  -gapopen 5  -gapextend 2 -out "${sample}_RM_others_blast_filtered"
	else
 		 echo "Required files ${sample}_RM_BLAST_MERGED_otherTEs.fasta or ${general_TE_lib} not found. Skipping BLAST."
  		 continue
	fi
	
	#processing blast output and saving as bed file; first position from previous bed
	if [[ -f "${sample}_RM_others_blast_filtered" ]]; then
		awk 'BEGIN {OFS="\t"} {split($1, a, "[:\\-]"); print a[1], a[2], a[3], $2, $3, $4, $5, $6, $7}' "${sample}_RM_others_blast_filtered" > "${sample}_RM_others_blast_filtered.bed"
	else
	echo "File ${sample}_RM_others_blast_filtered does not exist. Skipping processing."
	continue
fi
	
	#getting complete annotation
	#Getting strand info
	awk 'BEGIN {OFS="\t"} {split($1, a, "[:\\-]"); strand = ($10 > $11) ? "-" : "+"; print a[1], a[2], a[3], $3, ".", strand}' "${sample}_RM_Dsim_blast_filtered" > "${sample}_RM_Dsim_blast_first_step_final.bed"
        awk 'BEGIN {OFS="\t"} {split($1, a, "[:\\-]"); strand = ($10 > $11) ? "-" : "+"; print a[1], a[2], a[3], $3, ".", strand}' "${sample}_RM_others_blast_filtered" > "${sample}_RM_others_blast_second_step_final.bed"
        cat "${sample}_RM_Dsim_blast_first_step_final.bed" "${sample}_RM_others_blast_second_step_final.bed"  |bedtools sort -i - > "${sample}_final_TE_annotation.bed"
	awk '$3-$2>79' "${sample}_final_TE_annotation.bed" >"${sample}_final_TE_annotation_morestringent.bed"
	echo "Finished ${sample}"
	
	
	 cd "$initial_dir" 
	
done  < "$SAMPLE_LIST_FILE"
