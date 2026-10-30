#!/bin/bash
set -euo pipefail

# ----------------------------
# CONFIG
# ----------------------------
MIRNA_TSV="/data/adriana/TE_Invasions/piRNA/miRNA/miRNA_total_per_replicate.tsv"
BAM_DIR="/data/adriana/TE_Invasions/test_ping_pong/results_multimappers_corrected_lib"
TE_FAI="/data/adriana/TE_Invasions/test_ping_pong/TE_lib_1.6_correct.fasta.fai"

OUTDIR="/data/adriana/TE_Invasions/piRNA/miRNA/bw_TE_miRNA_norm"
BEDDIR="$OUTDIR/bedgraph"
CHROMSIZES="$OUTDIR/TE_lib_1.6_correct.chrom.sizes"

mkdir -p "$OUTDIR" "$BEDDIR"

# ----------------------------
# CHROM SIZES (TEs as "chromosomes")
# ----------------------------
cut -f1,2 "$TE_FAI" > "$CHROMSIZES"

echo "[INFO] Output dir: $OUTDIR"
echo "[INFO] Using miRNA normalizers from: $MIRNA_TSV"
echo "[INFO] Using chrom sizes: $CHROMSIZES"

# ----------------------------
# LOOP OVER BAMs 
# ----------------------------
for bam in "$BAM_DIR"/*.sorted.bam; do
  rep=$(basename "$bam" .sorted.bam)

  miRNA_total=$(awk -v r="$rep" '$1==r {print $2}' "$MIRNA_TSV" | head -n 1)

  if [[ -z "${miRNA_total}" ]]; then
    echo "[SKIP] $rep (no miRNA normalizer)"
    continue
  fi

  scale=$(awk -v m="$miRNA_total" 'BEGIN{printf("%.12f", 1000000.0/m)}')

  echo "[PROCESS] $rep miRNA_total=$miRNA_total scale_to_PPM=$scale"

  plus_bg="$BEDDIR/${rep}.plus.PPMmiRNA.bedGraph"
  minus_bg="$BEDDIR/${rep}.minus.PPMmiRNA.bedGraph"
  plus_bw="$OUTDIR/${rep}.plus.PPMmiRNA.bw"
  minus_bw="$OUTDIR/${rep}.minus.PPMmiRNA.bw"

  # PLUS strand bedGraph (scaled) 
  bedtools genomecov -bg -strand + -ibam "$bam" \
    | awk -v s="$scale" 'BEGIN{OFS="\t"} {$4=$4*s; print}' \
    | LC_ALL=C sort -k1,1 -k2,2n > "$plus_bg"

  # MINUS strand bedGraph (scaled)
  bedtools genomecov -bg -strand - -ibam "$bam" \
    | awk -v s="$scale" 'BEGIN{OFS="\t"} {$4=$4*s; print}' \
    | LC_ALL=C sort -k1,1 -k2,2n > "$minus_bg"

  #bedGraphToBigWig "$plus_bg" "$CHROMSIZES" "$plus_bw"
  #bedGraphToBigWig "$minus_bg" "$CHROMSIZES" "$minus_bw"
  
echo "[DONE] $rep"
done

echo "[INFO] All replicates processed. bedGraphs in: $BEDDIR"


  

