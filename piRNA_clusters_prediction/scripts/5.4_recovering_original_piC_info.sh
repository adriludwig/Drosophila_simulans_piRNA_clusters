#!/bin/bash

orig_dir="/data/adriana/TE_Invasions/Nanopore/cactus/cactus_singularity/piCs_bed"
outdir="/data/adriana/TE_Invasions/Nanopore/cactus/cactus_singularity/halLiftover/final_piCs_lifted"

mkdir -p "$outdir"

for hal_file in ./*_final_lifted.bed; do
    sample=$(basename "$hal_file" _final_lifted.bed)
    orig_file="$orig_dir/${sample}.reference_with_piC.txt"
    out_file="$outdir/${sample}_final_lifted_ed.bed"

    awk -F'\t' '
        NR==FNR { ids[$4] = $0; next }
        {
            for (id in ids) {
                pattern = "\\<" id "\\>"
                if ($0 ~ pattern) {
                    print ids[id] "\t" $0
                }
            }
        }
    ' "$hal_file" "$orig_file"|sed 's/;\t/;/g' > "$out_file"
done
