"""Standardize a single pvar file: remap IDs to rsids, keep only 5 core columns.

Reads an alignment CSV to build a chrom+pos+ref+alt -> rsid mapping, then
rewrites the input pvar to contain only #CHROM, POS, ID, REF, ALT columns
with rsid values in the ID column.
"""

import argparse
from pathlib import Path

import pandas as pd


def get_chrom_pos_ref_alt(df):
    return (
        df["chrom"].astype(str)
        + "_"
        + df["pos"].astype(str)
        + "_"
        + df["ref"].astype(str)
        + "_"
        + df["alt"].astype(str)
    )


def count_comment_lines(path):
    """Count leading lines that start with '##' (plink2 pvar metadata)."""
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                n += 1
            else:
                break
    return n


def main():
    parser = argparse.ArgumentParser(
        description="Standardize a pvar file with rsid mapping."
    )
    parser.add_argument(
        "--pvar", type=Path, required=True, help="Input pvar file path"
    )
    parser.add_argument(
        "--alignment", type=Path, required=True, help="Alignment CSV file path"
    )
    parser.add_argument(
        "--chrom", type=int, default=None,
        help="Chromosome number to filter alignment CSV",
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="Output pvar file path"
    )
    args = parser.parse_args()

    # Build cpra -> rsid mapping from alignment CSV
    df = pd.read_csv(args.alignment)
    if args.chrom is not None:
        df = df[df["chrom"] == args.chrom]
    df["cpra"] = get_chrom_pos_ref_alt(df)
    assert df["cpra"].is_unique, "Duplicate cpra values in alignment file"
    assert df["rsid"].is_unique, "Duplicate rsid values in alignment file"
    cpra_to_rsid = df.set_index("cpra")["rsid"].to_dict()

    # Read pvar, skipping ## comment lines
    skip = count_comment_lines(args.pvar)
    dv = pd.read_csv(args.pvar, sep="\t", skiprows=skip)

    dv.rename(
        columns={
            "#CHROM": "chrom",
            "POS": "pos",
            "REF": "ref",
            "ALT": "alt",
            "ID": "var_id",
        },
        inplace=True,
    )

    dv["cpra"] = get_chrom_pos_ref_alt(dv)

    missing = set(dv["cpra"]) - set(cpra_to_rsid.keys())
    assert not missing, f"{len(missing)} variants not found in alignment file"

    dv["rsid"] = dv["cpra"].map(cpra_to_rsid)

    dv = dv[["chrom", "pos", "rsid", "ref", "alt"]]
    dv = dv.rename(
        columns={
            "chrom": "#CHROM",
            "pos": "POS",
            "rsid": "ID",
            "ref": "REF",
            "alt": "ALT",
        }
    )

    dv.to_csv(args.output, sep="\t", index=False)
    print(f"  Mapped {len(dv)} variants")


if __name__ == "__main__":
    main()
