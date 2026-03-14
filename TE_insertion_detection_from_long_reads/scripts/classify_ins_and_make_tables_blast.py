#!/usr/bin/env python3
import argparse
import csv
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


@dataclass
class BlastHit:
    te_id: str
    pident: float
    aln_len: int
    bitscore: float
    evalue: float


def load_meta(meta_tsv: str) -> Dict[str, Dict[str, str]]:
    data: Dict[str, Dict[str, str]] = {}
    with open(meta_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = f"{row['CHROM']}|{row['POS']}|{row['ID']}"
            data[key] = row
    return data


def parse_blast_tsv(blast_tsv: str) -> Dict[str, BlastHit]:
    """
    Expects outfmt 6 fields:
      qseqid sseqid pident length evalue bitscore
    Returns best hit per qseqid based on:
      1) higher bitscore
      2) higher aligned length
      3) higher percent identity
    """
    best: Dict[str, BlastHit] = {}

    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        return best

    with open(blast_tsv, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            qseqid, sseqid, pident, length, evalue, bitscore = line.rstrip("\n").split("\t")
            hit = BlastHit(
                te_id=sseqid,
                pident=float(pident),
                aln_len=int(length),
                evalue=float(evalue),
                bitscore=float(bitscore),
            )

            if qseqid not in best:
                best[qseqid] = hit
            else:
                prev = best[qseqid]
                if (hit.bitscore, hit.aln_len, hit.pident) > (prev.bitscore, prev.aln_len, prev.pident):
                    best[qseqid] = hit

    return best


def write_tables(
    out_dir: str,
    sample: str,
    meta: Dict[str, Dict[str, str]],
    hits: Dict[str, BlastHit],
    min_ident: float,
    min_aln_len: int,
) -> None:
    os.makedirs(out_dir, exist_ok=True)

    te_rows: Dict[str, List[Dict[str, str]]] = {}
    unclassified: List[Dict[str, str]] = []

    cols = [
        "CHROM", "POS", "ID", "SVLEN", "SUPPORT", "SUPPORT_LONG", "AF",
        "GT", "GQ", "DR", "DV",
        "TE", "TE_pident", "TE_aln_len", "TE_evalue", "TE_bitscore",
    ]

    for key, row in meta.items():
        out_row = dict(row)

        if key not in hits:
            out_row["TE"] = "NO_HIT"
            out_row["TE_pident"] = ""
            out_row["TE_aln_len"] = ""
            out_row["TE_evalue"] = ""
            out_row["TE_bitscore"] = ""
            unclassified.append(out_row)
            continue

        h = hits[key]
        out_row["TE"] = h.te_id
        out_row["TE_pident"] = f"{h.pident:.2f}"
        out_row["TE_aln_len"] = str(h.aln_len)
        out_row["TE_evalue"] = f"{h.evalue:.2e}"
        out_row["TE_bitscore"] = f"{h.bitscore:.1f}"

        if (h.pident >= min_ident) and (h.aln_len >= min_aln_len):
            te_rows.setdefault(h.te_id, []).append(out_row)
        else:
            out_row["TE"] = "UNCLASSIFIED"
            unclassified.append(out_row)

    def write_tsv(path: str, rows: List[Dict[str, str]]) -> None:
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c, "") for c in cols})

    for te_id, rows in te_rows.items():
        te_safe = te_id.replace("/", "_").replace("#", "_")
        write_tsv(os.path.join(out_dir, f"{sample}.{te_safe}.tsv"), rows)

    write_tsv(os.path.join(out_dir, f"{sample}.UNCLASSIFIED.tsv"), unclassified)


def main() -> None:
    ap = argparse.ArgumentParser(description="Classify INS sequences by BLAST to TE DB and output per-TE TSV tables.")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--blast", required=True, help="BLAST outfmt 6 TSV")
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_ident", type=float, default=80.0, help="Minimum percent identity (default 80).")
    ap.add_argument("--min_aln_len", type=int, default=80, help="Minimum aligned length in bp (default 80).")
    args = ap.parse_args()

    meta = load_meta(args.meta)
    hits = parse_blast_tsv(args.blast)

    write_tables(
        out_dir=args.out_dir,
        sample=args.sample,
        meta=meta,
        hits=hits,
        min_ident=args.min_ident,
        min_aln_len=args.min_aln_len,
    )


if __name__ == "__main__":
    main()

