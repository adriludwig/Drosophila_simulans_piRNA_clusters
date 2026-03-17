#!/bin/bash
#This script soft-masks genome assemblies using TE annotation coordinates.  
# Paths correspond to the internal computing environment and should be adapted if the script is reused.

for fasta in /data/adriana/TE_Invasions/Nanopore/final_assemblies/*_final_assembly.fasta; do
    strain=$(basename "$fasta" _final_assembly.fasta)

    # TE annotation
    bed="/data/adriana/TE_Invasions/piRNA/repeatmasker/lib1.6/${strain}_final_TE_annotation_morestringent_new.bed"

    # output
    out="./softmasked_genomes/${strain}_final_assembly.softmasked.fasta"

    bedtools maskfasta -fi "$fasta" -bed "$bed" -soft -fo "$out"
done
