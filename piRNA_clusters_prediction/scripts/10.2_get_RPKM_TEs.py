import os
import pysam
import pandas as pd
from collections import defaultdict

# ----------------------------
# PATHS
# ----------------------------
miRNA_table = os.path.expanduser("/data/adriana/TE_Invasions/piRNA/miRNA/miRNA_total_per_replicate.tsv")
bam_dir = "/data/adriana/TE_Invasions/test_ping_pong/results_multimappers_corrected_lib"
te_fai = "/data/adriana/TE_Invasions/test_ping_pong/TE_lib_1.6_correct.fasta.fai"
out_dir = os.path.expanduser("~/TE_Invasions/piRNA/miRNA")

# ----------------------------
# LOAD miRNA NORMALIZERS
# ----------------------------
mirna_df = pd.read_csv(miRNA_table, sep="\t", header=None, names=["replicate", "total_miRNA"])
mirna_dict = dict(zip(mirna_df.replicate, mirna_df.total_miRNA))

# ----------------------------
# LOAD TE LENGTHS
# ----------------------------
te_len = {}
with open(te_fai) as fh:
    for line in fh:
        f = line.strip().split("\t")
        te_len[f[0]] = int(f[1])

# ----------------------------
# PROCESS EACH BAM
# ----------------------------
for bam_file in sorted(f for f in os.listdir(bam_dir) if f.endswith(".sorted.bam")):

    replicate = bam_file.replace(".sorted.bam", "")

    # skip contaminated / missing miRNA-normalizer replicates
    if replicate not in mirna_dict:
        print(f"[SKIP] {replicate} (no miRNA normalizer)")
        continue

    print(f"[PROCESS] {replicate}")

    bam_path = os.path.join(bam_dir, bam_file)
    bam = pysam.AlignmentFile(bam_path, "rb")

    sense = defaultdict(int)
    antisense = defaultdict(int)
    total_TE_reads = 0

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        te = bam.get_reference_name(read.reference_id)
        total_TE_reads += 1

        if read.is_reverse:
            antisense[te] += 1
        else:
            sense[te] += 1

    bam.close()

    rows = []
    miRNA_total = mirna_dict[replicate]

    for te in te_len:
        s = sense.get(te, 0)
        a = antisense.get(te, 0)
        length_kb = te_len[te] / 1000

        rpkm_s = (s * 1e6) / (length_kb * miRNA_total) if miRNA_total > 0 else 0
        rpkm_a = (a * 1e6) / (length_kb * miRNA_total) if miRNA_total > 0 else 0

        rows.append([
            replicate,
            te,
            te_len[te],
            s,
            a,
            s + a,
            rpkm_s,
            rpkm_a,
            miRNA_total,
            total_TE_reads
        ])

    out_df = pd.DataFrame(rows, columns=[
        "replicate",
        "TE_ID",
        "TE_length",
        "sense_reads",
        "antisense_reads",
        "total_TE_reads_per_TE",
        "RPKM_sense_miRNA_norm",
        "RPKM_antisense_miRNA_norm",
        "total_mapped_miRNA",
        "total_mapped_TE"
    ])

    out_file = os.path.join(out_dir, f"{replicate}.TE_RPKM_miRNA_norm.tsv")
    out_df.to_csv(out_file, sep="\t", index=False)

    print(f"[DONE] {out_file}")

print("All replicates processed.")

