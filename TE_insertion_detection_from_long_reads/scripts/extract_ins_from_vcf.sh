#!/usr/bin/env bash
set -euo pipefail

VCF="$1"
OUT_PREFIX="$2"

# 1) Extract insertion sequences (only if ALT is an actual DNA sequence, not <INS>)
# Output FASTA headers encode a stable key to join later:
# >CHROM|POS|ID
bcftools view -i 'INFO/SVTYPE="INS" && ALT!="<INS>"' "${VCF}" \
  | bcftools query -f '>%CHROM|%POS|%ID\n%ALT\n' \
  > "${OUT_PREFIX}.ins.fa"

# 2) Export a metadata table for ALL INS 
# Fields: CHROM POS ID SVLEN SUPPORT SUPPORT_LONG AF GT GQ DR DV
bcftools query -i 'INFO/SVTYPE="INS"' \
  -f '%CHROM\t%POS\t%ID\t%INFO/SVLEN\t%INFO/SUPPORT\t%INFO/SUPPORT_LONG\t%INFO/AF\t[%GT]\t[%GQ]\t[%DR]\t[%DV]\n' \
  "${VCF}" \
  > "${OUT_PREFIX}.ins.meta.tsv"

# Add a header
sed -i '1iCHROM\tPOS\tID\tSVLEN\tSUPPORT\tSUPPORT_LONG\tAF\tGT\tGQ\tDR\tDV' "${OUT_PREFIX}.ins.meta.tsv"

echo "[INFO] Wrote:"
echo "  ${OUT_PREFIX}.ins.fa"
echo "  ${OUT_PREFIX}.ins.meta.tsv"

