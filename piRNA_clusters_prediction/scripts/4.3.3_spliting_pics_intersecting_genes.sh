# Process each piC file
# Script to split piCs if they contain genes
# Paths correspond to the internal computing environment and should be adapted if the script is reused.

for pic_file in *_coverage_protrac_merged.bed; do
    sample="${pic_file%_coverage_protrac_merged.bed}"

    # Match the corresponding gene file
    gene_bed="../Gene_annot_strains/${sample}_final_assembly_liftoff_noTE_nolncRNA.gff.bed"

    # Output files
    outside="${sample}_pic_outside_gene.bed"
    outside_filtered="${sample}_pic_outside_gene_2kb.bed"
    outside_final="${sample}_pic_outside_gene_2kb_final.bed"
    te_annot="../../../TE_Curation_Annotation_Final/TE_annot/${sample}_final_TE_annotation_morestringent.bed"
    echo "Processing $sample"

    # Sort files for bedops
    sort-bed "$pic_file" > tmp_pic_sorted.bed
    sort-bed "$gene_bed" > tmp_gene_sorted.bed

    # Step 1: Remove piCs that overlap any gene
     bedtools subtract -a tmp_pic_sorted.bed -b tmp_gene_sorted.bed > "$outside"	

    #Step 2: Filter by size
    awk '$3-$2>1999' "$outside" >"$outside_filtered"

    #Step 4: Remove piCs that do no overlap with TEs
    bedtools intersect -a "${outside_filtered}" -b "${te_annot}" -wa |sort|uniq >"${outside_final}"
done

# Clean up temporary files
rm -f tmp_pic_sorted.bed tmp_gene_sorted.bed

