#!/usr/bin/env python3
import argparse
import csv
import subprocess
from typing import Dict, List, Optional, Tuple


def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    parts: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            seqs[name] = "".join(parts).upper()
    return seqs


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.upper().translate(comp)[::-1]


def tsd_len_from_te(te_name: str) -> Optional[int]:
    t = te_name.lower()
    if "shellder" in t:
        return 4
    if "spoink" in t:
        return 4
    # P element / P-derived naming
    if t == "p_ds" in t:
        return 8
    return None


def load_te_map(te_map_tsv: str) -> Dict[str, str]:
    """
    Requires columns: CHROM POS ID TE
    Returns key CHROM|POS|ID -> TE
    """
    te_for_key: Dict[str, str] = {}
    with open(te_map_tsv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"CHROM", "POS", "ID", "TE"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise SystemExit(f"TE map TSV must contain columns: {sorted(required)}")
        for row in reader:
            key = f"{row['CHROM']}|{row['POS']}|{row['ID']}"
            te_for_key[key] = row["TE"]
    return te_for_key


def fetch_ref_window(ref_fa: str, chrom: str, start: int, end: int) -> str:
    region = f"{chrom}:{start}-{end}"
    res = subprocess.run(
        ["samtools", "faidx", ref_fa, region],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    seq = "".join([l.strip() for l in res.stdout.splitlines() if l.strip() and not l.startswith(">")]).upper()
    return seq


def ref_kmers_around_pos(ref_fa: str, chrom: str, pos_1based: int, k: int) -> Tuple[str, str]:
    """
    left_k  = ref[POS-k .. POS-1]
    right_k = ref[POS .. POS+k-1]
    """
    start = max(1, pos_1based - k)
    end = pos_1based + k - 1
    seq = fetch_ref_window(ref_fa, chrom, start, end)

    idx = pos_1based - start  # 0-based index corresponding to POS
    left_k = seq[max(0, idx - k):idx]
    right_k = seq[idx:idx + k]

    if len(left_k) < k:
        left_k = ("N" * (k - len(left_k))) + left_k
    if len(right_k) < k:
        right_k = right_k + ("N" * (k - len(right_k)))

    return left_k, right_k


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Extract TSD candidates from ALT ends and match to REF around POS (with POS±W)."
    )
    ap.add_argument("--ins_fa", required=True, help="FASTA headers CHROM|POS|ID, seq=ALT.")
    ap.add_argument("--te_map", required=True, help="TSV with CHROM POS ID TE (e.g. sample.P.tsv).")
    ap.add_argument("--ref_fa", required=True, help="Reference FASTA indexed with samtools faidx.")
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--pos_window", type=int, default=2, help="Test POS offsets in [-W..+W], default 2.")
    args = ap.parse_args()

    ins = read_fasta(args.ins_fa)
    te_map = load_te_map(args.te_map)

    # Output includes raw kmers + inferred final TSD match
    out_cols = [
        "key", "CHROM", "POS", "ID", "TE", "TSD_len",
        "ALT_first_k", "ALT_last_k",
        "REF_left_k_at_POS", "REF_right_k_at_POS",
        "MATCH_any",
        "FINAL_TSD_ALT", "FINAL_TSD_REF",
        "FINAL_ALT_END", "FINAL_REF_SIDE", "FINAL_POS_OFFSET",
        "ALL_MATCHES",
    ]

    with open(args.out_tsv, "w", encoding="utf-8", newline="") as out:
        out.write("\t".join(out_cols) + "\n")

        # Only evaluate calls present in te_map (e.g. only P calls if you pass sample.P.tsv)
        for key, te in te_map.items():
            chrom, pos_str, vid = key.split("|", 2)
            pos = int(pos_str)

            k = tsd_len_from_te(te)
            if k is None:
                out.write("\t".join([
                    key, chrom, pos_str, vid, te, "",
                    "", "", "", "",
                    "NA",
                    "", "", "", "", "",
                    "No_TSD_length_for_TE"
                ]) + "\n")
                continue

            if key not in ins:
                # Usually ALT=<INS> (no inserted sequence in ins.fa)
                out.write("\t".join([
                    key, chrom, pos_str, vid, te, str(k),
                    "", "", "", "",
                    "NA",
                    "", "",
                    "", "", "",
                    "No_inserted_sequence_in_ins_fa_(ALT_<INS>_or_filtered)"
                ]) + "\n")
                continue

            alt_seq = ins[key]
            alt_first = alt_seq[:k]
            alt_last = alt_seq[-k:]

            # Raw ref kmers at POS (offset 0) for transparency
            ref_left_0, ref_right_0 = ref_kmers_around_pos(args.ref_fa, chrom, pos, k)

            alt_candidates = {
                "ALT_first": alt_first,
                "ALT_last": alt_last,
                "ALT_first_rc": reverse_complement(alt_first),
                "ALT_last_rc": reverse_complement(alt_last),
            }

            best = None  # tuple: (off_abs, off, bits, alt_name, ref_side, tsd_alt, tsd_ref)
            matches: List[str] = []

            for off in range(-args.pos_window, args.pos_window + 1):
                test_pos = pos + off
                if test_pos < 1:
                    continue
                ref_left, ref_right = ref_kmers_around_pos(args.ref_fa, chrom, test_pos, k)
                ref_candidates = {"REF_left": ref_left, "REF_right": ref_right}

                for alt_name, tsd_alt in alt_candidates.items():
                    if "N" in tsd_alt:
                        continue
                    for ref_side, tsd_ref in ref_candidates.items():
                        if "N" in tsd_ref:
                            continue
                        if tsd_alt == tsd_ref:
                            matches.append(f"{alt_name}=={ref_side}@offset{off}")

                            # Choose "best" match:
                            # Primary: smallest absolute offset (prefer exact POS)
                            # Secondary: prefer non-RC over RC (arbitrary but stable)
                            # Tertiary: prefer ALT_last over ALT_first (arbitrary but stable)
                            off_abs = abs(off)
                            pref_rc = 1 if alt_name.endswith("_rc") else 0
                            pref_end = 0 if "ALT_last" in alt_name else 1  # ALT_last preferred
                            score = (off_abs, pref_rc, pref_end)

                            if best is None or score < best[0]:
                                best = (score, off, alt_name, ref_side, tsd_alt, tsd_ref, ref_left, ref_right)

            match_any = "YES" if matches else "NO"
            all_matches = ",".join(matches) if matches else "no_match"

            if best is None:
                out.write("\t".join([
                    key, chrom, pos_str, vid, te, str(k),
                    alt_first, alt_last,
                    ref_left_0, ref_right_0,
                    "NO",
                    "", "",
                    "", "", "",
                    all_matches
                ]) + "\n")
            else:
                _, off, alt_name, ref_side, tsd_alt, tsd_ref, ref_left_best, ref_right_best = best
                out.write("\t".join([
                    key, chrom, pos_str, vid, te, str(k),
                    alt_first, alt_last,
                    ref_left_0, ref_right_0,
                    "YES",
                    tsd_alt, tsd_ref,
                    alt_name, ref_side, str(off),
                    all_matches
                ]) + "\n")


if __name__ == "__main__":
    main()

