"""Generate per-chromosome DRAGEN variant ID files from SBayesRC alignment data.

Reads the combined alignment CSV (sbayesrc_hg38.csv) containing all chromosomes
and writes one text file per chromosome (chr1.txt … chr22.txt) into
data/dragen_ids/, each containing DRAGEN-format variant IDs
(e.g. DRAGEN:chr1:12345:A:G).

Idempotent: skips chromosomes whose output file already exists.
"""

import argparse
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILE = SCRIPT_DIR / "data" / "support" / "sbayesrc_hg38.csv"
OUTPUT_DIR = SCRIPT_DIR / "data" / "dragen_ids"


def make_dragen_id(row: pd.Series) -> str:
    return f"DRAGEN:chr{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-file", type=Path, default=INPUT_FILE,
        help="Alignment CSV file (all chromosomes, columns: chrom,pos,ref,alt,rsid)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=OUTPUT_DIR,
        help="Directory to write per-chromosome DRAGEN ID files",
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    full_df = pd.read_csv(args.input_file)

    for chrom in range(1, 23):
        output_path = args.output_dir / f"chr{chrom}.txt"
        if output_path.exists() and output_path.stat().st_size > 0:
            print(f"chr{chrom}: skipping — {output_path} already exists")
            continue

        df = full_df[full_df["chrom"] == chrom]
        df["dragen_id"] = df.apply(make_dragen_id, axis=1)

        ids = df["dragen_id"].tolist()
        output_path.write_text("\n".join(ids) + "\n")
        print(f"chr{chrom}: wrote {len(ids)} variant IDs to {output_path}")


if __name__ == "__main__":
    main()
