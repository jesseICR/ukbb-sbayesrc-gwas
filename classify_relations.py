#!/usr/bin/env python3
# classify_relations.py — Classify close relationships from KING kinship output.
#
# Reads the KING .kin0 file, classifies close relationships, joins with
# birth year/month data to add fractional year-of-birth columns, and
# filters out pairs with large age differences and low IBS0.
#
# Relationship types:
#   - sibling:      first-degree kinship AND IBS0 >= 0.0012
#   - parent_child: first-degree kinship AND IBS0 <  0.0012
#   - identical:    kinship >= 2^{-1.5}  AND IBS0 <  0.0012 (twins/duplicates)
#
# Filter: removes implausible *sibling* pairs (yob_diff > 20 AND IBS0 < 0.002)
#
# Input:
#   /mnt/project${DX_OUTPUT_DIR}/kinship/ukb_all_direct_rel.kin0
#   birth_year_month.csv  (extracted from UKBB fields 34 + 52)
# Output:
#   close_relations.csv   (auto-uploaded by SAK)

import csv
import os
import subprocess
from datetime import date

DX_OUTPUT_DIR = os.environ.get("DX_OUTPUT_DIR", "/sbayesrc_genotypes")
_BASE = f"/mnt/project{DX_OUTPUT_DIR}"

KINSHIP_PATH = f"{_BASE}/kinship/ukb_all_direct_rel.kin0"
BIRTH_PATH = "birth_year_month.csv"
OUTPUT_PATH = "close_relations.csv"
DATASET_NAME_GLOB = "app*.dataset"

# Kinship thresholds (Manichaikul et al. 2010; Bycroft et al. 2018)
FIRST_DEGREE_LOWER = 0.1767       # ~2^{-2.5}
FIRST_DEGREE_UPPER = 0.3535       # ~2^{-1.5}
IBS0_CUTOFF        = 0.0012       # UKBB official processing threshold

# Age-gap filter thresholds
YOB_DIFF_THRESHOLD  = 20
IBS0_FILTER_CUTOFF  = 0.002


def fractional_yob(year, month):
    """Approximate fractional year of birth using the 15th of the month."""
    d = date(year, month, 15)
    days_in_year = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    return year + (d.timetuple().tm_yday - 1) / days_in_year


# ── Discover dataset and extract birth data ──────────────────────────────────
project_id = os.environ["DX_PROJECT_CONTEXT_ID"]
result = subprocess.run(
    ["dx", "find", "data", "--name", DATASET_NAME_GLOB,
     "--project", project_id, "--brief"],
    capture_output=True, text=True, check=True,
)
# --brief returns "project-xxx:record-xxx", use as-is
dataset_id = result.stdout.strip().split("\n")[0]
print(f"Discovered dataset: {dataset_id}")

subprocess.run(
    ["dx", "extract_dataset", dataset_id,
     "--fields", "participant.eid,participant.p34,participant.p52",
     "--delimiter", ",",
     "--output", BIRTH_PATH],
    check=True,
)
print(f"Extracted birth data ({sum(1 for _ in open(BIRTH_PATH)) - 1} participants)")

# ── Load birth data ──────────────────────────────────────────────────────────
birth = {}
with open(BIRTH_PATH) as f:
    reader = csv.DictReader(f)
    for row in reader:
        eid = row["participant.eid"]
        year_str = row.get("participant.p34", "").strip()
        month_str = row.get("participant.p52", "").strip()
        if year_str and month_str:
            birth[eid] = fractional_yob(int(year_str), int(month_str))

print(f"Loaded birth data for {len(birth)} participants")

# ── Classify relationships ───────────────────────────────────────────────────
rows = []
counts_before = {"sibling": 0, "parent_child": 0, "identical": 0}

with open(KINSHIP_PATH) as fin:
    for line in fin:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        fields = line.split()
        # Fields: FID1 IID1 FID2 IID2 NSNP HETHET IBS0 KINSHIP
        iid1    = fields[1]
        iid2    = fields[3]
        hethet  = float(fields[5])
        ibs0    = float(fields[6])
        kinship = float(fields[7])

        if kinship >= FIRST_DEGREE_UPPER and ibs0 < IBS0_CUTOFF:
            rel = "identical"
        elif kinship >= FIRST_DEGREE_LOWER and kinship < FIRST_DEGREE_UPPER:
            rel = "sibling" if ibs0 >= IBS0_CUTOFF else "parent_child"
        else:
            continue

        yob1 = birth.get(iid1)
        yob2 = birth.get(iid2)
        counts_before[rel] += 1
        rows.append((iid1, iid2, kinship, ibs0, hethet, rel, yob1, yob2))

total_before = sum(counts_before.values())
print(f"\nBefore filtering: {total_before} close relationships")
for rel, n in counts_before.items():
    print(f"  {rel}: {n}")

# ── Filter: remove implausible sibling pairs (large age gap + low IBS0) ──────
filtered = []
counts_after = {"sibling": 0, "parent_child": 0, "identical": 0}
removed = 0

for row in rows:
    iid1, iid2, kinship, ibs0, hethet, rel, yob1, yob2 = row
    if rel == "sibling" and yob1 is not None and yob2 is not None:
        if abs(yob1 - yob2) > YOB_DIFF_THRESHOLD and ibs0 < IBS0_FILTER_CUTOFF:
            removed += 1
            continue
    counts_after[rel] += 1
    filtered.append(row)

total_after = sum(counts_after.values())
print(f"\nAfter filtering (removed {removed} pairs with yob_diff > {YOB_DIFF_THRESHOLD} & IBS0 < {IBS0_FILTER_CUTOFF}):")
print(f"  Total: {total_after}")
for rel, n in counts_after.items():
    print(f"  {rel}: {n}")

# ── Write output ─────────────────────────────────────────────────────────────
with open(OUTPUT_PATH, "w", newline="") as fout:
    writer = csv.writer(fout)
    writer.writerow(["eid1", "eid2", "kinship", "ibs0", "hethet", "relationship", "yob1", "yob2"])
    for iid1, iid2, kinship, ibs0, hethet, rel, yob1, yob2 in filtered:
        writer.writerow([
            iid1, iid2, kinship, ibs0, hethet, rel,
            f"{yob1:.4f}" if yob1 is not None else "",
            f"{yob2:.4f}" if yob2 is not None else "",
        ])
