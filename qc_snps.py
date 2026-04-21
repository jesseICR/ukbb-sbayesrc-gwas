"""QC-filter the WGS vs TopMed allele-frequency comparison to a high-quality SNP set.

Three filters are applied (SBayesRC paper, Zheng et al. 2024, Nat Genet):
  1. Low MAF in WGS — drops rare-variant noise.
  2. WGS vs TopMed-imputed freq diff — drops variants where TopMed imputation in
     UKB disagrees with the WGS truth on our own EUR samples.
  3. WGS vs HRC-imputed (SBayesRC) freq diff — drops variants where our WGS
     frequencies disagree with the SBayesRC paper's own HRC-imputed reference
     panel (white-British EUR, hg19; lifted to hg38 and allele-aligned here).

Input:
  - freq-compare CSV produced by compute_freq_compare.sh, columns:
      chrom, variant_id, ref, alt, wgs_alt_cts, wgs_obs_ct, imp_alt_cts, imp_obs_ct
  - SBayesRC liftover CSV (jesseICR/sbayesrc-liftover release v1.0), columns
    used here: ID, A1, A2, dbsnp_ref, A1Freq.
    Note: A1/A2 are the original hg19 SBayesRC alleles; dbsnp_ref is the hg38
    reference allele; A1Freq is the SBayesRC HRC-imputed frequency of A1.

Output:
  CSV copy of the input freq-compare rows that pass all three filters,
  annotated with the derived columns alt_freq_wgs, alt_freq_topmed, abs_diff_wgs_vs_topmed,
  and alt_freq_hrc (SBayesRC's HRC-imputed frequency of the hg38 ALT allele).
"""

import argparse
import sys

import pandas as pd


# Watson–Crick complement — only used for palindromic / strand-flipped SNPs
# where neither hg19 A1 nor A2 matches the hg38 dbsnp_ref directly.
COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--freq-compare", required=True,
                   help="Per-variant WGS vs TopMed-imputed freq-compare CSV")
    p.add_argument("--sbayesrc-liftover", required=True,
                   help="SBayesRC liftover CSV with hg19→hg38 allele mapping and A1Freq")
    p.add_argument("--maf-threshold", type=float, required=True,
                   help="Drop variants with WGS MAF strictly below this threshold")
    p.add_argument("--freq-diff-threshold", type=float, required=True,
                   help="Drop variants whose |WGS − TopMed-imputed| alt-freq diff exceeds this")
    p.add_argument("--sbayesrc-freq-diff-threshold", type=float, required=True,
                   help="Drop variants whose |WGS − HRC-imputed (SBayesRC)| alt-freq diff exceeds this")
    p.add_argument("--output", required=True, help="Output QC-passed CSV path")
    return p.parse_args()


def compute_wgs_and_topmed_frequencies(variants: pd.DataFrame) -> pd.DataFrame:
    """Add alt_freq_wgs, alt_freq_topmed, and abs_diff_wgs_vs_topmed columns."""
    variants["alt_freq_wgs"] = variants.wgs_alt_cts / variants.wgs_obs_ct
    variants["alt_freq_topmed"] = variants.imp_alt_cts / variants.imp_obs_ct
    variants["abs_diff_wgs_vs_topmed"] = (
        variants.alt_freq_topmed - variants.alt_freq_wgs
    ).abs()
    return variants


def compute_hrc_alt_frequency(sbayesrc: pd.DataFrame) -> pd.DataFrame:
    """Align SBayesRC's hg19 A1Freq onto the hg38 ALT allele defined by dbsnp_ref.

    For every variant, exactly one of the following holds:
      - A2 == dbsnp_ref  → A1 is the ALT allele, so alt_freq_hrc = A1Freq
      - A1 == dbsnp_ref  → A2 is the ALT allele, so alt_freq_hrc = 1 − A1Freq
      - Neither (palindromic A/T or C/G with strand flip, or true mismatch):
          complement A1; if complement(A1) == dbsnp_ref, then A1 (on the hg38
          strand) represents the REF, so alt_freq_hrc = 1 − A1Freq. Otherwise
          A1 represents the ALT on the hg38 strand and alt_freq_hrc = A1Freq.
    """
    match_forward = sbayesrc.dbsnp_ref == sbayesrc.A2
    match_reverse = sbayesrc.dbsnp_ref == sbayesrc.A1
    match_neither = ~(match_forward | match_reverse)

    # Default: A2 == dbsnp_ref (most common) → A1 is ALT → alt_freq_hrc = A1Freq.
    sbayesrc["alt_freq_hrc"] = sbayesrc.A1Freq
    sbayesrc.loc[match_reverse, "alt_freq_hrc"] = 1 - sbayesrc.loc[match_reverse, "A1Freq"]

    # Strand-flip / palindromic handling: compare complement(A1) against dbsnp_ref.
    a1_complement = sbayesrc.A1.map(COMPLEMENT)
    needs_complement_flip = match_neither & (a1_complement == sbayesrc.dbsnp_ref)
    sbayesrc.loc[needs_complement_flip, "alt_freq_hrc"] = (
        1 - sbayesrc.loc[needs_complement_flip, "A1Freq"]
    )

    print("  Allele-alignment counts in liftover table:")
    print(f"    A2 == dbsnp_ref (alt_freq_hrc = A1Freq):         {int(match_forward.sum()):>10,}")
    print(f"    A1 == dbsnp_ref (alt_freq_hrc = 1 − A1Freq):     {int(match_reverse.sum()):>10,}")
    print(f"    Neither matches (complement/strand-flip handled): {int(match_neither.sum()):>10,}")
    return sbayesrc


def main():
    args = parse_args()

    print(f"Loading freq-compare CSV: {args.freq_compare}")
    variants = pd.read_csv(args.freq_compare)
    print(f"  {len(variants):,} variants loaded")

    variants = compute_wgs_and_topmed_frequencies(variants)

    # --- Filter 1: low MAF in WGS ---
    mask_low_maf = variants.alt_freq_wgs.map(
        lambda f: min(f, 1 - f) < args.maf_threshold
    )

    # --- Filter 2: WGS vs TopMed-imputed freq diff ---
    mask_fail_topmed = variants.abs_diff_wgs_vs_topmed > args.freq_diff_threshold

    # --- Compute alt_freq_hrc from SBayesRC liftover, then Filter 3 ---
    print(f"\nLoading SBayesRC liftover CSV: {args.sbayesrc_liftover}")
    sbayesrc = pd.read_csv(args.sbayesrc_liftover)
    print(f"  {len(sbayesrc):,} rows loaded")

    sbayesrc = sbayesrc[sbayesrc.ID.isin(variants.variant_id)].copy()
    print(f"  {len(sbayesrc):,} rows after intersecting with freq-compare variants")

    sbayesrc = compute_hrc_alt_frequency(sbayesrc)

    variants["alt_freq_hrc"] = variants.variant_id.map(
        sbayesrc.set_index("ID").alt_freq_hrc
    )
    mask_fail_hrc = (
        variants.alt_freq_wgs - variants.alt_freq_hrc
    ).abs() > args.sbayesrc_freq_diff_threshold

    mask_fail_any = mask_low_maf | mask_fail_topmed | mask_fail_hrc

    # --- Report ---
    filter_1_label = f"Filter 1 — WGS MAF < {args.maf_threshold}"
    filter_2_label = f"Filter 2 — |WGS − TopMed-imputed| alt-freq > {args.freq_diff_threshold}"
    filter_3_label = f"Filter 3 — |WGS − HRC-imputed (SBayesRC)| alt-freq > {args.sbayesrc_freq_diff_threshold}"
    label_width = max(len(filter_1_label), len(filter_2_label), len(filter_3_label)) + 2

    print(f"\nQC filter results:")
    print(f"  {filter_1_label:<{label_width}} {int(mask_low_maf.sum()):>10,} flagged")
    print(f"  {filter_2_label:<{label_width}} {int(mask_fail_topmed.sum()):>10,} flagged")
    print(f"  {filter_3_label:<{label_width}} {int(mask_fail_hrc.sum()):>10,} flagged")
    print(f"  {'Overlap: flagged by BOTH TopMed AND HRC freq-diff':<{label_width}} "
          f"{int((mask_fail_topmed & mask_fail_hrc).sum()):>10,}")
    print(f"  {'Overlap: low MAF AND WGS-vs-TopMed freq diff':<{label_width}} "
          f"{int((mask_low_maf & mask_fail_topmed).sum()):>10,}")
    print(f"  {'Total flagged by ANY of the 3 filters':<{label_width}} "
          f"{int(mask_fail_any.sum()):>10,}")

    kept = variants.loc[~mask_fail_any]
    print(f"\nFinal: {len(kept):,} variants pass QC → {args.output}")

    if kept["alt_freq_hrc"].isna().any():
        n_missing = int(kept["alt_freq_hrc"].isna().sum())
        print(
            f"WARNING: {n_missing:,} kept variants have no SBayesRC alt_freq_hrc "
            "(not present in liftover CSV). These passed the HRC freq-diff filter "
            "trivially because |NaN| is not > threshold.",
            file=sys.stderr,
        )

    kept.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
