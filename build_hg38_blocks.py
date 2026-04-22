"""Derive hg38 LD-block boundaries from SBayesRC's per-SNP Block assignments.

The sbayesrc-liftover release (v1.0) provides a per-SNP table with both
hg19 and hg38 positions and the original SBayesRC `Block` id. That already
tells us, for every SBayesRC SNP, which 4 cM block it belongs to. We can
therefore skip interval liftover entirely and define each block's hg38
bounds directly from its SNPs' hg38 positions.

For each block we compute min_i = min(pos_hg38) and max_i = max(pos_hg38)
across its SNPs, then choose StartBP and EndBP that strictly bracket every
SNP in the block:

    StartBP_i = min_i - 1, decremented further while StartBP_i coincides
                with any SBayesRC panel SNP position (for defensive
                strictness against dense SNP regions).
    EndBP_i   = max_i + 1, incremented further under the same rule.

That yields `StartBP_i < pos < EndBP_i` *strictly* for every SNP in block
i, with no SNP ever equal to any boundary — independent of which
half-open convention SBayesRC's `LDstep1` uses internally.

Rows with `pos_hg38 == -1` (unlifted) or empty `Block` are ignored.

Output: one line per block in SBayesRC's space-delimited format,
renumbered 1..N in `(chrom, StartBP)` order:

    Block Chrom StartBP EndBP
    1 1 866280 3665398
    ...
"""

import argparse
import bisect
import csv
from collections import defaultdict


def chrom_sort_key(c):
    """Sort autosomes numerically; anything non-numeric sorts after."""
    return (0, int(c)) if c.isdigit() else (1, c)


def verify_round_trip(liftover_csv, alignment_csv, pos_file):
    """Round-trip QC on the block file we just wrote, using the canonical
    alignment file (data/support/sbayesrc_hg38.csv) as the hg38 position
    source — NOT the liftover CSV we derived the blocks from.

    Two independent inputs:
      * `liftover_csv`    — supplies only the original SBayesRC `Block` id
                            per rsid (the hg19 4 cM partition).
      * `alignment_csv`   — supplies the canonical hg38 (chrom, pos) per
                            rsid. This is the same file the rest of the
                            pipeline uses for hg38 coords.

    For every SNP in the alignment file:
      1. Bisect-assign it into a new block using its alignment (chrom, pos).
      2. Look up its original Block from the liftover CSV.
      3. Group new-ids by original Block id.

    Passes iff:
      - 0 SNPs fall outside any block.
      - 0 SNPs land exactly on a boundary.
      - Every original Block's SNPs all bisect-assign to exactly one new id.
      - The resulting original->new mapping is a 1-to-1 bijection over the
        full set of new ids in the .pos file.
      - The two source files agree on (chrom, pos_hg38) for every rsid
        they share.

    Any failure raises SystemExit.
    """

    # ── Re-read the .pos file into per-chrom bisect arrays. ────────────────
    starts_by_chrom = defaultdict(list)
    ends_by_chrom   = defaultdict(list)
    ids_by_chrom    = defaultdict(list)
    with open(pos_file) as f:
        next(f)  # skip header
        for line in f:
            blk, chrom, s, e = line.split()
            starts_by_chrom[chrom].append(int(s))
            ends_by_chrom[chrom].append(int(e))
            ids_by_chrom[chrom].append(int(blk))
    for chrom, starts in starts_by_chrom.items():
        assert starts == sorted(starts), f"chrom {chrom} starts not sorted"

    # ── rsid -> (orig Block, chrom, pos_hg38) from the liftover CSV. ──────
    # We only need Block for the verification itself; chrom/pos are kept
    # solely to cross-check agreement with the alignment CSV below.
    orig_block = {}          # rsid -> int
    liftover_coord = {}      # rsid -> (chrom, pos_hg38)
    with open(liftover_csv) as f:
        for row in csv.DictReader(f):
            if not row["Block"]:
                continue
            try:
                p = int(row["pos_hg38"])
            except (ValueError, TypeError):
                continue
            if p < 0:
                continue
            rsid = row["ID"]
            orig_block[rsid] = int(row["Block"])
            liftover_coord[rsid] = (row["chrom"], p)

    # ── Scan the alignment CSV, cross-check, and bisect. ──────────────────
    n_total = 0
    n_missing_block = 0    # rsid not in liftover CSV (unlifted / dropped)
    n_pos_mismatch = 0     # alignment (chrom, pos) disagrees with liftover
    n_outside = 0          # bisect places SNP outside any block
    n_on_boundary = 0      # SNP lands exactly on a boundary
    orig_to_new_ids = defaultdict(set)   # orig Block -> set of new ids
    mismatch_examples = []

    with open(alignment_csv) as f:
        for row in csv.DictReader(f):
            rsid  = row["rsid"]
            chrom = row["chrom"]
            pos   = int(row["pos"])
            n_total += 1

            if rsid not in orig_block:
                n_missing_block += 1
                continue

            # Cross-check source-file agreement.
            lc, lp = liftover_coord[rsid]
            if lc != chrom or lp != pos:
                n_pos_mismatch += 1
                if len(mismatch_examples) < 5:
                    mismatch_examples.append(
                        f"{rsid}: alignment=({chrom},{pos}) "
                        f"liftover=({lc},{lp})")

            starts = starts_by_chrom.get(chrom)
            ends   = ends_by_chrom.get(chrom)
            ids    = ids_by_chrom.get(chrom)
            if not starts:
                n_outside += 1
                continue
            i = bisect.bisect_right(starts, pos) - 1
            if i < 0 or pos >= ends[i]:
                n_outside += 1
                continue
            if pos == starts[i] or pos == ends[i]:
                n_on_boundary += 1
                continue
            orig_to_new_ids[orig_block[rsid]].add(ids[i])

    # ── Derive the orig->new mapping and bijection metrics. ────────────────
    split_origs  = {o: s for o, s in orig_to_new_ids.items() if len(s) > 1}
    orig_to_new  = {o: next(iter(s)) for o, s in orig_to_new_ids.items()
                    if len(s) == 1}
    new_ids_used = set(orig_to_new.values())
    identity     = sum(1 for o, n in orig_to_new.items() if o == n)
    expected_ids = set(b for ids in ids_by_chrom.values() for b in ids)

    # ── Report. ───────────────────────────────────────────────────────────
    print("\n=== Round-trip QC: alignment-file positions vs derived blocks ===")
    print(f"  Position source:   {alignment_csv}")
    print(f"  Block-id source:   {liftover_csv}")
    print(f"  Boundary file:     {pos_file}")
    print()
    print(f"  SNPs in alignment CSV:                         {n_total:,}")
    print(f"  SNPs with no Block id in the liftover CSV:     {n_missing_block}")
    print(f"  SNPs with chrom/pos disagreement across files: {n_pos_mismatch}")
    print(f"  SNPs outside any block:                        {n_outside}")
    print(f"  SNPs landing on a block boundary:              {n_on_boundary}")
    print(f"  Original blocks split across >1 new block:     {len(split_origs)}")
    print(f"  Original blocks mapping cleanly to one new id: {len(orig_to_new)}")
    print(f"  Distinct new ids reached by the mapping:       {len(new_ids_used)}")
    print(f"  Blocks where new_id == original_id:            "
          f"{identity} / {len(orig_to_new)}")
    if mismatch_examples:
        print("  First chrom/pos mismatches:")
        for m in mismatch_examples:
            print(f"    {m}")

    # ── Gates. ────────────────────────────────────────────────────────────
    if n_pos_mismatch:
        raise SystemExit(
            f"FAIL: {n_pos_mismatch} rsids have disagreeing (chrom, pos) "
            f"between {alignment_csv} and {liftover_csv}")
    if n_outside:
        raise SystemExit(
            f"FAIL: {n_outside} alignment-file SNPs fell outside all blocks")
    if n_on_boundary:
        raise SystemExit(
            f"FAIL: {n_on_boundary} alignment-file SNPs landed on a block "
            f"boundary; strict-bounds guarantee violated")
    if split_origs:
        raise SystemExit(
            f"FAIL: {len(split_origs)} original blocks split across multiple "
            f"new blocks. Examples: {list(split_origs.items())[:3]}")
    if len(orig_to_new) != len(new_ids_used):
        raise SystemExit(
            "FAIL: original->new block mapping is not 1-to-1 (two originals "
            "collide on the same new id)")
    if new_ids_used != expected_ids:
        missing = expected_ids - new_ids_used
        raise SystemExit(
            f"FAIL: mapping doesn't cover every new block; missing {missing}")

    print("  PASS: alignment-file positions preserve SBayesRC's block "
          "partition under the derived hg38 boundaries")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--liftover", required=True,
                    help="sbayesrc_liftover_results.csv")
    ap.add_argument("--alignment", required=True,
                    help="data/support/sbayesrc_hg38.csv "
                         "(canonical hg38 positions used by the pipeline)")
    ap.add_argument("--output", required=True,
                    help="Output .pos file")
    args = ap.parse_args()

    # Per-block state keyed by the original SBayesRC Block id.
    # Values: {"chrom": str, "min": int, "max": int, "n": int}
    state = {}
    # Per-chrom set of occupied SBayesRC panel positions (used below for the
    # shift-until-free boundary expansion).
    panel_positions = defaultdict(set)
    n_rows = 0
    n_no_block = 0
    n_unlifted = 0
    n_used = 0

    with open(args.liftover) as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_rows += 1
            blk_str = row["Block"]
            if not blk_str:
                n_no_block += 1
                continue
            try:
                pos = int(row["pos_hg38"])
            except (ValueError, TypeError):
                n_unlifted += 1
                continue
            if pos < 0:
                # sbayesrc-liftover uses -1 as the unlifted sentinel
                n_unlifted += 1
                continue
            blk = int(blk_str)
            chrom = row["chrom"]
            n_used += 1
            panel_positions[chrom].add(pos)
            s = state.get(blk)
            if s is None:
                state[blk] = {"chrom": chrom, "min": pos, "max": pos, "n": 1}
                continue
            if s["chrom"] != chrom:
                raise SystemExit(
                    f"Block {blk} has SNPs on multiple chroms: "
                    f"{s['chrom']} and {chrom}")
            if pos < s["min"]:
                s["min"] = pos
            if pos > s["max"]:
                s["max"] = pos
            s["n"] += 1

    print(f"Rows read:                {n_rows:,}")
    print(f"  with empty Block:       {n_no_block:,}")
    print(f"  with unlifted pos_hg38: {n_unlifted:,}")
    print(f"  used:                   {n_used:,}")
    print(f"Distinct blocks kept:     {len(state)}")

    # Expand each block's StartBP and EndBP away from any panel-SNP position
    # so every SNP strictly satisfies StartBP < pos < EndBP.
    start_shifts = 0
    end_shifts = 0
    for s in state.values():
        occupied = panel_positions[s["chrom"]]
        start = s["min"] - 1
        while start in occupied:
            start -= 1
            start_shifts += 1
        end = s["max"] + 1
        while end in occupied:
            end += 1
            end_shifts += 1
        s["startbp"] = start
        s["endbp"] = end
    print(f"Boundary shifts to clear panel-SNP collisions: "
          f"{start_shifts} at StartBP, {end_shifts} at EndBP")

    # Sort blocks by (chrom, StartBP) and renumber 1..N.
    records = []
    for _orig_blk, s in state.items():
        records.append((chrom_sort_key(s["chrom"]), s["startbp"], s["endbp"],
                        s["chrom"], s["n"]))
    records.sort(key=lambda r: (r[0], r[1]))

    # Safety check: adjacent blocks on the same chrom must not overlap.
    for i in range(len(records) - 1):
        _cka, s_a, e_a, ch_a, _ = records[i]
        _ckb, s_b, e_b, ch_b, _ = records[i + 1]
        if ch_a == ch_b and e_a >= s_b:
            raise SystemExit(
                f"Adjacent blocks on chr{ch_a} overlap after boundary "
                f"expansion: block ending at {e_a} meets block starting at "
                f"{s_b} (blocks must satisfy EndBP_i < StartBP_{{i+1}}).")

    with open(args.output, "w") as out:
        out.write("Block Chrom StartBP EndBP\n")
        for i, (_ck, start, end, chrom, _n) in enumerate(records, 1):
            out.write(f"{i} {chrom} {start} {end}\n")
    print(f"Wrote {args.output} ({len(records)} blocks)")

    # Per-chromosome summary.
    per_chrom = defaultdict(list)
    for _ck, start, end, chrom, n in records:
        per_chrom[chrom].append((start, end, n))

    print("\nPer-chromosome summary:")
    print(f"  {'chrom':>5}  {'blocks':>6}  {'span_low':>12}  {'span_high':>12}  {'snps':>10}")
    total_snps = 0
    for chrom in sorted(per_chrom, key=chrom_sort_key):
        spans = sorted(per_chrom[chrom])
        n_snps = sum(s[2] for s in spans)
        total_snps += n_snps
        print(f"  {chrom:>5}  {len(spans):>6}  {spans[0][0]:>12,}  "
              f"{spans[-1][1]:>12,}  {n_snps:>10,}")
    print(f"  {'total':>5}  {sum(len(v) for v in per_chrom.values()):>6}  "
          f"{'':>12}  {'':>12}  {total_snps:>10,}")

    # Final QC: verify the written .pos file against the canonical alignment
    # file (hg38 positions used by the rest of the pipeline), using the
    # liftover CSV only as the source of original Block ids. Fails loudly
    # on any inconsistency.
    verify_round_trip(args.liftover, args.alignment, args.output)


if __name__ == "__main__":
    main()
