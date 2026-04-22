#!/usr/bin/env python3
# sample_ld_cohort.py - Sample LD-reference cohort for custom SBayesRC LD.
#
# Algorithm:
#   1. Load fit_pca_iids.txt (unrelated Europeans for PCA fitting)
#   2. Load sex_covar.txt (IIDs with a genetic-sex call; aneuploidy excluded)
#   3. Extract White British (UKBB field 22006 == 1) via dx extract_dataset
#   4. Intersect all three sets
#   5. Random-sample LD_COHORT_SIZE IIDs with a fixed seed for reproducibility
#
# Input:
#   /mnt/project${DX_OUTPUT_DIR}/pca_eur/fit_pca_iids.txt
#   /mnt/project${DX_OUTPUT_DIR}/genetic_sex/sex_covar.txt
#   UK Biobank field 22006 (extracted via dx extract_dataset)
#
# Output:
#   ld_ref_40k_iids.txt    (FID IID, space-separated, FID==IID, sorted numerically)
#   ld_ref_cohort_log.txt  (step-by-step counts)

import csv
import os
import random
import subprocess
import sys

DX_OUTPUT_DIR = os.environ.get("DX_OUTPUT_DIR", "/sbayesrc_genotypes")
LD_COHORT_SIZE = int(os.environ["LD_COHORT_SIZE"])
RANDOM_SEED = int(os.environ["RANDOM_SEED"])

_BASE = f"/mnt/project{DX_OUTPUT_DIR}"
FIT_PCA_PATH = f"{_BASE}/pca_eur/fit_pca_iids.txt"
SEX_COVAR_PATH = f"{_BASE}/genetic_sex/sex_covar.txt"

DATASET_NAME_GLOB = "app*.dataset"

OUTPUT_FILE = "ld_ref_40k_iids.txt"
LOG_FILE = "ld_ref_cohort_log.txt"

log_lines = []


def log(msg):
    print(msg)
    log_lines.append(msg)


def write_iid_file(path, iids):
    """Write FID IID file (FID == IID), sorted numerically."""
    sorted_iids = sorted(iids, key=lambda x: int(x))
    with open(path, "w") as f:
        for iid in sorted_iids:
            f.write(f"{iid} {iid}\n")


# ── 1. Load fit_pca_iids.txt ───────────────────────────────────────────────
fit_pca_iids = set()
with open(FIT_PCA_PATH) as f:
    for line in f:
        parts = line.strip().split()
        if parts:
            fit_pca_iids.add(parts[1])
log(f"fit_pca_iids.txt: {len(fit_pca_iids)} IIDs")

# ── 2. Load sex_covar.txt (header: FID IID sex_01, tab-separated) ──────────
sex_covar_iids = set()
with open(SEX_COVAR_PATH) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sex_covar_iids.add(row["IID"])
log(f"sex_covar.txt: {len(sex_covar_iids)} IIDs")

# ── 3. Extract White British (UKBB field 22006) ────────────────────────────
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
     "--fields", "participant.eid,participant.p22006",
     "--delimiter", ",",
     "--output", "white_british.csv"],
    check=True,
)

white_british = set()
with open("white_british.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row.get("participant.p22006", "").strip() == "1":
            white_british.add(row["participant.eid"])
log(f"White British (field 22006 == 1): {len(white_british)}")

# ── 4. Intersect ───────────────────────────────────────────────────────────
log("")
log("=== Intersections ===")
log(f"fit_pca & sex_covar:       {len(fit_pca_iids & sex_covar_iids)}")
log(f"fit_pca & white_british:   {len(fit_pca_iids & white_british)}")
log(f"sex_covar & white_british: {len(sex_covar_iids & white_british)}")

eligible = fit_pca_iids & sex_covar_iids & white_british
log(f"Eligible (all three):      {len(eligible)}")

if len(eligible) < LD_COHORT_SIZE:
    log(f"\nERROR: eligible pool ({len(eligible)}) smaller than requested "
        f"cohort size ({LD_COHORT_SIZE})")
    with open(LOG_FILE, "w") as f:
        f.write("\n".join(log_lines) + "\n")
    sys.exit(1)

# ── 5. Random sample with fixed seed ───────────────────────────────────────
log("")
log(f"=== Random sample (seed={RANDOM_SEED}, size={LD_COHORT_SIZE}) ===")
rng = random.Random(RANDOM_SEED)
sampled = set(rng.sample(sorted(eligible, key=int), LD_COHORT_SIZE))
log(f"Sampled: {len(sampled)} / {len(eligible)}")

# ── 6. Write output ────────────────────────────────────────────────────────
write_iid_file(OUTPUT_FILE, sampled)
log(f"\nWrote {OUTPUT_FILE} ({len(sampled)} IIDs)")

with open(LOG_FILE, "w") as f:
    f.write("\n".join(log_lines) + "\n")

# ── 7. Verification ────────────────────────────────────────────────────────
print("\n=== Verification checks ===")
all_passed = True

check_fit_pca = sampled - fit_pca_iids
passed = len(check_fit_pca) == 0
all_passed &= passed
print(f"Verification - subset_of_fit_pca: {'PASS' if passed else 'FAIL'}")
print(f"  {len(sampled)} IIDs checked, {len(check_fit_pca)} not in fit_pca_iids")

check_sex = sampled - sex_covar_iids
passed = len(check_sex) == 0
all_passed &= passed
print(f"Verification - subset_of_sex_covar: {'PASS' if passed else 'FAIL'}")
print(f"  {len(sampled)} IIDs checked, {len(check_sex)} not in sex_covar")

check_wb = sampled - white_british
passed = len(check_wb) == 0
all_passed &= passed
print(f"Verification - all_white_british: {'PASS' if passed else 'FAIL'}")
print(f"  {len(sampled)} IIDs checked, {len(check_wb)} not White British")

passed = len(sampled) == LD_COHORT_SIZE
all_passed &= passed
print(f"Verification - cohort_size: {'PASS' if passed else 'FAIL'}")
print(f"  {len(sampled)} == {LD_COHORT_SIZE}")

if not all_passed:
    print("\nERROR: one or more verification checks failed")
    sys.exit(1)

# ── 8. Cleanup ─────────────────────────────────────────────────────────────
os.remove("white_british.csv")
os.remove("sample_ld_cohort.py")
