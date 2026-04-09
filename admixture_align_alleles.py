#!/usr/bin/env python3
"""
Align alleles between the reference allele frequency TSV and the extracted
PLINK bim, then produce a correctly ordered .P file for ADMIXTURE projection.

The .P file must have rows in the same order as the .bim, with each row
containing allele frequencies for the bim's A1 allele across K populations.

Reads from the working directory:
    admixture_allele_freqs.tsv  — reference allele frequencies
    ukb_extracted.bim           — extracted PLINK bim

Writes to the working directory:
    ref_aligned.P               — allele frequency matrix (space-separated, no header)
    snps_aligned.txt            — rsIDs of retained SNPs (for plink2 --extract)
    admixture_align_log.txt     — alignment summary
"""
import csv

POP_COLS = ["European", "East Asian", "American", "African", "South Asian", "Oceanian"]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Load reference TSV into dict keyed by snp_id
tsv = {}
with open("admixture_allele_freqs.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        tsv[row["snp_id"]] = row

# Load extracted bim
bim_rows = []
with open("ukb_extracted.bim") as f:
    for line in f:
        parts = line.strip().split("\t")
        bim_rows.append({
            "chr": parts[0], "rsid": parts[1], "gendist": parts[2],
            "pos": parts[3], "a1": parts[4], "a2": parts[5],
        })

# Align alleles
p_rows = []
kept_snps = []
n_same, n_flipped, n_excluded = 0, 0, 0

for row in bim_rows:
    rsid = row["rsid"]
    if rsid not in tsv:
        n_excluded += 1
        continue

    ref = tsv[rsid]
    bim_a1, bim_a2 = row["a1"], row["a2"]
    tsv_a1, tsv_a2 = ref["a1"], ref["a2"]
    freqs = [float(ref[col]) for col in POP_COLS]

    # Skip strand-ambiguous SNPs (A/T or C/G pairs)
    allele_set = {bim_a1, bim_a2}
    if allele_set == {"A", "T"} or allele_set == {"C", "G"}:
        n_excluded += 1
        continue

    if bim_a1 == tsv_a1 and bim_a2 == tsv_a2:
        p_rows.append(freqs)
        kept_snps.append(rsid)
        n_same += 1
    elif bim_a1 == tsv_a2 and bim_a2 == tsv_a1:
        p_rows.append([1.0 - f for f in freqs])
        kept_snps.append(rsid)
        n_flipped += 1
    elif bim_a1 == COMPLEMENT.get(tsv_a1) and bim_a2 == COMPLEMENT.get(tsv_a2):
        p_rows.append(freqs)
        kept_snps.append(rsid)
        n_same += 1
    elif bim_a1 == COMPLEMENT.get(tsv_a2) and bim_a2 == COMPLEMENT.get(tsv_a1):
        p_rows.append([1.0 - f for f in freqs])
        kept_snps.append(rsid)
        n_flipped += 1
    else:
        n_excluded += 1

summary = (
    f"Allele alignment: {n_same} same, {n_flipped} flipped, {n_excluded} excluded\n"
    f"{len(kept_snps)} SNPs retained out of {len(bim_rows)} extracted\n"
)
print(summary, end="")

# Write aligned SNP list
with open("snps_aligned.txt", "w") as f:
    for s in kept_snps:
        f.write(s + "\n")

# Write .P file in bim order
with open("ref_aligned.P", "w") as f:
    for freqs in p_rows:
        f.write(" ".join(str(v) for v in freqs) + "\n")
print(f"Wrote ref_aligned.P with {len(p_rows)} rows x {len(POP_COLS)} populations")

# Write alignment log
with open("admixture_align_log.txt", "w") as f:
    f.write(summary)
    f.write(f"Populations: {', '.join(POP_COLS)}\n")
