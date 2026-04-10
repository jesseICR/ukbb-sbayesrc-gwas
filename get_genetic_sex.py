"""Build genetic sex covariate file from WGS fam + UKB field 22001.

Algorithm:
  1. Extract fields 22001 (genetic sex) and 22019 (sex chr aneuploidy) from UKBB
  2. Read fam file and imputed-only IID list
  3. Exclude individuals with sex chromosome aneuploidy (field 22019 == 1)
  4. For imputed-only IIDs: derive sex from field 22001 (0=Female, 1=Male)
  5. For WGS IIDs: derive sex from fam (1=male, 2=female)
  6. Write sex_covar.txt with coding: 0=female, 1=male

Input:
  /mnt/project${DX_OUTPUT_DIR}/direct_bfile/chr1_22_merged.fam
  /mnt/project${DX_OUTPUT_DIR}/merge_steps/imputed_only_iids.txt
  UKBB fields: participant.p22001, participant.p22019

Output (auto-uploaded by SAK):
  sex_covar.txt        — FID  IID  sex_01
  readme.txt           — method and coding description
  genetic_sex_log.txt  — counts at each step
"""

import csv
import os
import subprocess
import sys

DATASET_NAME_GLOB = "app*.dataset"

DX_OUTPUT_DIR = os.environ.get("DX_OUTPUT_DIR", "/sbayesrc_genotypes")
_BASE = f"/mnt/project{DX_OUTPUT_DIR}"

FAM_PATH = f"{_BASE}/direct_bfile/chr1_22_merged.fam"
IMPUTED_ONLY_PATH = f"{_BASE}/merge_steps/imputed_only_iids.txt"
EXTRACT_PATH = "ukb_sex_fields.csv"
OUTPUT_PATH = "sex_covar.txt"
README_PATH = "readme.txt"
LOG_PATH = "genetic_sex_log.txt"

log_lines = []


def log(msg):
    print(msg)
    log_lines.append(msg)


# ── Discover dataset and extract fields ──────────────────────────────────────
project_id = os.environ["DX_PROJECT_CONTEXT_ID"]
result = subprocess.run(
    ["dx", "find", "data", "--name", DATASET_NAME_GLOB,
     "--project", project_id, "--brief"],
    capture_output=True, text=True, check=True,
)
dataset_id = result.stdout.strip().split("\n")[0]
log(f"Discovered dataset: {dataset_id}")

subprocess.run(
    ["dx", "extract_dataset", dataset_id,
     "--fields", "participant.eid,participant.p22001,participant.p22019",
     "--delimiter", ",",
     "--output", EXTRACT_PATH],
    check=True,
)
log(f"Extracted fields 22001 + 22019 ({sum(1 for _ in open(EXTRACT_PATH)) - 1} rows)")

# ── Load extracted fields ────────────────────────────────────────────────────
field_22001 = {}  # eid -> int (0=Female, 1=Male)
aneuploidy_eids = set()

with open(EXTRACT_PATH) as f:
    reader = csv.DictReader(f)
    for row in reader:
        eid = row["participant.eid"]
        sex_val = row["participant.p22001"]
        aneuploidy_val = row["participant.p22019"]
        if sex_val != "":
            field_22001[eid] = int(sex_val)
        if aneuploidy_val == "1":
            aneuploidy_eids.add(eid)

log(f"\nField 22001: {len(field_22001)} participants with genetic sex")
log(f"Field 22019: {len(aneuploidy_eids)} participants with sex chromosome aneuploidy")

# Validate field 22001 values
field_22001_values = set(field_22001.values())
assert field_22001_values <= {0, 1}, f"Unexpected values in field 22001: {field_22001_values}"
log(f"Field 22001 unique values: {sorted(field_22001_values)} (validated: only 0, 1)")

# ── Load fam file ────────────────────────────────────────────────────────────
fam_sex = {}  # iid -> int (1=male, 2=female)
with open(FAM_PATH) as f:
    for line in f:
        parts = line.strip().split()
        iid = parts[1]
        sex = int(parts[4])
        fam_sex[iid] = sex

log(f"\nFam file: {len(fam_sex)} total IIDs")

# Validate fam sex values
fam_sex_values = set(fam_sex.values())
assert fam_sex_values <= {1, 2}, f"Unexpected values in fam sex column: {fam_sex_values}"
log(f"Fam sex unique values: {sorted(fam_sex_values)} (validated: only 1, 2)")

# ── Load imputed-only IIDs ───────────────────────────────────────────────────
imputed_only = set()
with open(IMPUTED_ONLY_PATH) as f:
    for line in f:
        parts = line.strip().split()
        if parts:
            imputed_only.add(parts[0])

log(f"Imputed-only IIDs: {len(imputed_only)}")

# ── Exclude aneuploidy ───────────────────────────────────────────────────────
fam_iids = set(fam_sex.keys())
aneuploidy_in_fam = fam_iids & aneuploidy_eids
remaining_iids = fam_iids - aneuploidy_eids

log(f"\nAneuploidy IIDs in fam: {len(aneuploidy_in_fam)}")
log(f"Remaining after aneuploidy exclusion: {len(remaining_iids)}")

# ── Assign sex ───────────────────────────────────────────────────────────────
sex_results = {}  # iid -> int (0=female, 1=male)
from_field = 0
from_fam = 0
missing_field = []

for iid in remaining_iids:
    if iid in imputed_only:
        # Get sex from field 22001
        if iid not in field_22001:
            missing_field.append(iid)
            continue
        sex_results[iid] = field_22001[iid]  # already 0/1
        from_field += 1
    else:
        # Get sex from fam: 1=male -> 1, 2=female -> 0
        fam_val = fam_sex[iid]
        sex_results[iid] = 1 if fam_val == 1 else 0
        from_fam += 1

log(f"\nSex from field 22001 (imputed-only): {from_field}")
log(f"Sex from fam (WGS): {from_fam}")
if missing_field:
    log(f"WARNING: {len(missing_field)} imputed-only IIDs missing from field 22001: {missing_field[:10]}")

males = sum(1 for v in sex_results.values() if v == 1)
females = sum(1 for v in sex_results.values() if v == 0)
log(f"\nFinal counts:")
log(f"  Total: {len(sex_results)}")
log(f"  Male (1): {males}")
log(f"  Female (0): {females}")
log(f"  Male/Female ratio: {males/females:.3f}")

# ── Write sex_covar.txt ─────────────────────────────────────────────────────
sorted_iids = sorted(sex_results.keys(), key=lambda x: int(x))
with open(OUTPUT_PATH, "w") as f:
    f.write("FID\tIID\tsex_01\n")
    for iid in sorted_iids:
        f.write(f"{iid}\t{iid}\t{sex_results[iid]}\n")

log(f"\nWrote {len(sex_results)} rows to {OUTPUT_PATH}")

# ── Write readme.txt ────────────────────────────────────────────────────────
with open(README_PATH, "w") as f:
    f.write("""Genetic Sex Covariate File (sex_covar.txt)
==========================================

Columns: FID  IID  sex_01

Coding:
  0 = Female
  1 = Male

Method:
  Starting from all individuals in the direct_bfile chr1_22_merged.fam:

  1. Exclusion: Individuals flagged in UK Biobank field 22019 (sex chromosome
     aneuploidy) are excluded entirely from this file.

  2. Imputed-only individuals: For individuals present in
     imputed_only_iids.txt (those without WGS data), genetic sex is obtained
     from UK Biobank field 22001 (Genetic sex), which codes 0=Female, 1=Male.
     These individuals do not have reliable sex entries in the fam file.

  3. WGS individuals: For all other individuals, sex is obtained from the
     direct_bfile chr1_22_merged.fam file, which codes 1=male, 2=female
     (converted to 1=male, 0=female for this file).
""")

log(f"Wrote {README_PATH}")

# ── Verification ─────────────────────────────────────────────────────────────
print("\n=== Verification checks ===")
all_passed = True

# Check 1: output unique values
output_values = set(sex_results.values())
if output_values == {0, 1}:
    print("PASS: sex_01 unique values are {0, 1}")
else:
    print(f"FAIL: sex_01 unique values are {output_values}")
    all_passed = False

# Check 2: no aneuploidy IIDs in output
aneuploidy_in_output = aneuploidy_eids & set(sex_results.keys())
if len(aneuploidy_in_output) == 0:
    print("PASS: no aneuploidy IIDs in output")
else:
    print(f"FAIL: {len(aneuploidy_in_output)} aneuploidy IIDs in output")
    all_passed = False

# Check 3: all imputed-only IIDs (minus aneuploidy) in output
imputed_expected = imputed_only - aneuploidy_eids
imputed_in_output = imputed_expected & set(sex_results.keys())
if len(imputed_in_output) == len(imputed_expected):
    print(f"PASS: all {len(imputed_expected)} imputed-only IIDs (excl. aneuploidy) in output")
else:
    print(f"FAIL: {len(imputed_in_output)}/{len(imputed_expected)} imputed-only IIDs in output")
    all_passed = False

# Check 4: all WGS IIDs (minus aneuploidy) in output
wgs_iids = fam_iids - imputed_only
wgs_expected = wgs_iids - aneuploidy_eids
wgs_in_output = wgs_expected & set(sex_results.keys())
if len(wgs_in_output) == len(wgs_expected):
    print(f"PASS: all {len(wgs_expected)} WGS IIDs (excl. aneuploidy) in output")
else:
    print(f"FAIL: {len(wgs_in_output)}/{len(wgs_expected)} WGS IIDs in output")
    all_passed = False

# Check 5: total count
expected_total = len(fam_iids) - len(aneuploidy_in_fam)
if len(sex_results) == expected_total:
    print(f"PASS: total {len(sex_results)} == fam ({len(fam_iids)}) - aneuploidy ({len(aneuploidy_in_fam)})")
else:
    print(f"FAIL: total {len(sex_results)} != expected {expected_total}")
    all_passed = False

# Check 6: no missing field lookups
if len(missing_field) == 0:
    print("PASS: all imputed-only IIDs found in field 22001")
else:
    print(f"FAIL: {len(missing_field)} imputed-only IIDs missing from field 22001")
    all_passed = False

if all_passed:
    print("\nAll verification checks passed.")
else:
    print("\nWARNING: Some verification checks failed!")
    sys.exit(1)

# ── Write log ────────────────────────────────────────────────────────────────
with open(LOG_PATH, "w") as f:
    f.write("\n".join(log_lines) + "\n")

# ── Cleanup ──────────────────────────────────────────────────────────────────
os.remove(EXTRACT_PATH)
os.remove("get_genetic_sex.py")
