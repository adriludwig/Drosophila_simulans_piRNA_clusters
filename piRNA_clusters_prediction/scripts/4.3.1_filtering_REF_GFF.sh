#!/bin/bash
# This script removes lncRNA genes and genes whose exons show TE matches.
# Paths correspond to the internal computing environment and should be adapted if the scripts are reused.

###########eliminating lncRNAS##############################################################
#selecting IDs corresponding to lnRNA genes
awk -F'\t' '$3 == "gene" && /gene_biotype=lncRNA/' GCF_016746395.2_Prin_Dsim_3.1_genomic.gff | \
sed -n 's/.*ID=\([^;]*\).*/\1/p' > gene_ids_from_lncRNAs.txt

#eliminating lncRNAs from GFF
agat_sp_filter_feature_from_kill_list.pl \
  --gff GCF_016746395.2_Prin_Dsim_3.1_genomic.gff \
  --kill_list gene_ids_from_lncRNAs.txt \
  --output GCF_filtered_no_lncRNA.gff

############# Removing genes with TE matches in their exons ##########################∑
#  INPUT FILES 
ref_gff="GCF_filtered_no_lncRNA.gff"
ref_genome="/data/adriana/TE_Invasions/references/GCF_016746395.2_Prin_Dsim_3.1_genomic.fna"
TEs="/data/adriana/TE_Invasions/piRNA/TE_library/Manual_Curated_Dsim_TE_lib.fasta"

#  OUTPUT FILES 
exons_gff="GCF_016746395.2_exons_only.gff"
exons_fasta="GCF_016746395.2_exons_only.fasta"
blast_db="exons_blast_db"
blast_output="TEs_vs_exons_initial_blast.tsv"

#  STEP 1: Extract only exon features from GFF 
echo "Extracting exons from GFF..."
grep -P "\texon\t" "$ref_gff" > "$exons_gff"

#  STEP 2: Extract exon sequences from the genome 
echo "Extracting exon sequences from genome..."
gffread "$exons_gff" -g "$ref_genome" -w "$exons_fasta"

#  STEP 3: Create a BLAST database from exon sequences 
echo "Creating BLAST database..."
makeblastdb -in "$exons_fasta" -dbtype nucl -out "$blast_db"

#  STEP 4: Run BLAST of TEs against exon sequences 
echo "Running BLAST of TEs vs exons..."
blastn -query "$TEs" -db "$blast_db" -out "$blast_output" \
  -evalue 1e-10 -word_size 11 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

#  STEP 5: Report 
hit_count=$(wc -l < "$blast_output")
echo " Done. Found $hit_count TE hits to exons in the original GFF."
echo " BLAST results saved to: $blast_output"

#  STEP 6: Selecting gene IDs
cut -f2 TEs_vs_exons_initial_blast.tsv | sort -u > exon_ids_to_remove_2.txt

#  STEP 7: Eliminating genes matching TEs
agat_sp_filter_feature_from_kill_list.pl \
  --gff $ref_gff \
  --kill_list exon_ids_to_remove_2.txt\
  --output GCF_filtered_no_lncRNA_noTE.gff