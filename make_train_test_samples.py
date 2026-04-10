#!/usr/bin/env python3
# make_train_test_samples.py — Build train/test sample split for REGENIE.
#
# Splits European-ancestry individuals into training and test samples,
# ensuring no cryptic relatedness (3rd degree or closer) between sets.
#
# Algorithm:
#   1. Initial test  = European & White British (field 22006) & siblings
#   2. Initial train = European - initial test
#   3. Build kinship graph from .kin0 (pairs with KINSHIP >= 0.0441941)
#   4. Secondary train = initial train + European neighbours of initial train
#   5. Secondary test  = initial test - secondary train
#   6. Final train     = secondary train - neighbours of secondary test
#      Final test      = secondary test (unchanged)
#
# Input:
#   /mnt/project${DX_OUTPUT_DIR}/europeans/classified_european_iids.txt
#   /mnt/project${DX_OUTPUT_DIR}/kinship/close_relations.csv
#   /mnt/project${DX_OUTPUT_DIR}/kinship/ukb_all_direct_rel.kin0
#   UK Biobank field 22006 (extracted via dx extract_dataset)
#
# Output:
#   final_train_iids.txt  (auto-uploaded by SAK)
#   final_test_iids.txt   (auto-uploaded by SAK)
#   train_test_log.txt    (uploaded to train_test/)
#   scrap/verification/*  (uploaded to train_test/scrap/verification/)

import csv
import os
import subprocess
from collections import defaultdict

DX_OUTPUT_DIR = os.environ.get("DX_OUTPUT_DIR", "/sbayesrc_genotypes")
_BASE = f"/mnt/project{DX_OUTPUT_DIR}"

EUROPEANS_PATH = f"{_BASE}/europeans/classified_european_iids.txt"
CLOSE_RELATIONS_PATH = f"{_BASE}/kinship/close_relations.csv"
KIN0_PATH = f"{_BASE}/kinship/ukb_all_direct_rel.kin0"

KINSHIP_THRESHOLD = 0.0441941  # 0.5^(9/2), lower bound for 3rd-degree relatives

DATASET_NAME_GLOB = "app*.dataset"

TRAIN_OUTPUT = "final_train_iids.txt"
TEST_OUTPUT = "final_test_iids.txt"
LOG_FILE = "train_test_log.txt"

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


# ── 1. Extract White British field 22006 ────────────────────────────────────
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
log(f"White British (field 22006): {len(white_british)}")

# ── 2. Load European IIDs ───────────────────────────────────────────────────
european_iids = set()
with open(EUROPEANS_PATH) as f:
    for line in f:
        parts = line.strip().split()
        if parts:
            european_iids.add(parts[1])
log(f"European IIDs (our classification): {len(european_iids)}")

# ── 3. Load siblings from close_relations.csv ──────────────────────────────
siblings = set()
with open(CLOSE_RELATIONS_PATH) as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row["relationship"] == "sibling":
            siblings.add(row["eid1"])
            siblings.add(row["eid2"])
log(f"Siblings (from close_relations.csv): {len(siblings)}")

# ── 4. Build kinship graph from .kin0 ───────────────────────────────────────
kinship_graph = defaultdict(set)
n_edges = 0
with open(KIN0_PATH) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        iid1 = fields[1]
        iid2 = fields[3]
        kinship = float(fields[7])
        if kinship >= KINSHIP_THRESHOLD:
            kinship_graph[iid1].add(iid2)
            kinship_graph[iid2].add(iid1)
            n_edges += 1
log(f"Kinship graph: {len(kinship_graph)} nodes, {n_edges} edges (KINSHIP >= {KINSHIP_THRESHOLD})")

# ── 5. Algorithm ────────────────────────────────────────────────────────────
log("")
log("=== Train/test split algorithm ===")

# Step 1: Initial test = European & White British & siblings
initial_test = european_iids & white_british & siblings
log(f"Step 1 — Initial test (European & WB & siblings): {len(initial_test)}")

eur_and_wb = european_iids & white_british
eur_and_sib = european_iids & siblings
wb_and_sib = white_british & siblings
log(f"  European & WB: {len(eur_and_wb)}")
log(f"  European & siblings: {len(eur_and_sib)}")
log(f"  WB & siblings: {len(wb_and_sib)}")

# Step 2: Initial train = European - initial_test
initial_train = european_iids - initial_test
log(f"Step 2 — Initial train (European - initial_test): {len(initial_train)}")

# Step 4: Secondary train = initial_train + European neighbours
secondary_train = set(initial_train)
added_via_expansion = 0
for iid in initial_train:
    for neighbor in kinship_graph.get(iid, set()):
        if neighbor in european_iids and neighbor not in secondary_train:
            secondary_train.add(neighbor)
            added_via_expansion += 1
log(f"Step 4 — Secondary train: {len(secondary_train)} (+{added_via_expansion} Europeans added via kinship expansion)")

from_initial_test_absorbed = initial_test & secondary_train
log(f"  Of which {len(from_initial_test_absorbed)} were in the initial test set (absorbed into train)")

# Step 5: Secondary test = initial_test - secondary_train
secondary_test = initial_test - secondary_train
removed_from_test = len(initial_test) - len(secondary_test)
log(f"Step 5 — Secondary test: {len(secondary_test)} (removed {removed_from_test} from initial test)")

# Step 6: Final train = secondary_train - neighbours of secondary_test
test_neighbors = set()
for iid in secondary_test:
    test_neighbors.update(kinship_graph.get(iid, set()))
removed_from_train = secondary_train & test_neighbors
final_train = secondary_train - test_neighbors
log(f"Step 6 — Final train: {len(final_train)} (removed {len(removed_from_train)} neighbours of secondary test)")

final_test = secondary_test
log("")
log("=== Final results ===")
log(f"Final train: {len(final_train)}")
log(f"Final test:  {len(final_test)}")

in_neither = european_iids - final_train - final_test
log(f"Europeans in neither set: {len(in_neither)}")

# ── 6. Write output files ──────────────────────────────────────────────────
write_iid_file(TRAIN_OUTPUT, final_train)
write_iid_file(TEST_OUTPUT, final_test)
log(f"\nWrote {TRAIN_OUTPUT} ({len(final_train)} IIDs)")
log(f"Wrote {TEST_OUTPUT} ({len(final_test)} IIDs)")

# ── 7. Write log ─────────────────────────────────────────────────────────────
# Log file stays in working directory — SAK auto-uploads it to --destination
with open(LOG_FILE, "w") as f:
    f.write("\n".join(log_lines) + "\n")

# ── 8. Verification (printed to stdout / job log) ───────────────────────────
print("\n=== Verification checks ===")

all_passed = True

# Check 1: Every final_test IID is a sibling
not_siblings = final_test - siblings
passed = len(not_siblings) == 0
all_passed &= passed
print(f"Verification — siblings: {'PASS' if passed else 'FAIL'}")
print(f"  {len(final_test)} test IIDs checked, {len(not_siblings)} not in siblings")
if not_siblings:
    print(f"  First 10: {sorted(list(not_siblings))[:10]}")

# Check 2: Every final_test IID is White British
not_wb = final_test - white_british
passed = len(not_wb) == 0
all_passed &= passed
print(f"Verification — white_british: {'PASS' if passed else 'FAIL'}")
print(f"  {len(final_test)} test IIDs checked, {len(not_wb)} not White British")
if not_wb:
    print(f"  First 10: {sorted(list(not_wb))[:10]}")

# Check 3: Every final_test IID is European
not_eur = final_test - european_iids
passed = len(not_eur) == 0
all_passed &= passed
print(f"Verification — european: {'PASS' if passed else 'FAIL'}")
print(f"  {len(final_test)} test IIDs checked, {len(not_eur)} not European")
if not_eur:
    print(f"  First 10: {sorted(list(not_eur))[:10]}")

# Check 4: No final_test IID has a 3rd+ degree relative in final_train
violations = []
for iid in final_test:
    related_in_train = kinship_graph.get(iid, set()) & final_train
    if related_in_train:
        violations.append((iid, len(related_in_train)))
passed = len(violations) == 0
all_passed &= passed
print(f"Verification — kinship_separation: {'PASS' if passed else 'FAIL'}")
print(f"  {len(final_test)} test IIDs checked, {len(violations)} have relatives in train")
if violations:
    print(f"  First 10: {violations[:10]}")

if not all_passed:
    print("\nWARNING: Some verification checks failed!")

# ── 9. Cleanup ──────────────────────────────────────────────────────────────
# Remove files that should NOT be auto-uploaded by SAK.
# Keep final_train_iids.txt, final_test_iids.txt, train_test_log.txt for auto-upload.
os.remove("white_british.csv")
os.remove("make_train_test_samples.py")
