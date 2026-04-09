#!/usr/bin/env python3
# select_pca_europeans.py — Select unrelated European IIDs for fitting PCA.
#
# Algorithm:
#   1. Start with classified European IIDs
#   2. Identify White British siblings (field 22006 == 1 AND relationship == "sibling")
#   3. Expand WB siblings to include all their relatives (kinship >= 0.0441941)
#   4. Remove expanded set from Europeans
#   5. Remove imputed-only IIDs (lower-quality samples)
#   6. Apply plink2 --king-cutoff-table to get maximal unrelated set
#
# Input:
#   /mnt/project/sbayesrc_genotypes/europeans/classified_european_iids.txt
#   /mnt/project/sbayesrc_genotypes/kinship/close_relations.csv
#   /mnt/project/sbayesrc_genotypes/kinship/ukb_all_direct_rel.kin0
#   /mnt/project/sbayesrc_genotypes/merge_steps/imputed_only_iids.txt
#   /mnt/project/sbayesrc_genotypes/direct_bfile/chr1_22_merged.{bed,bim,fam}
#   UK Biobank field 22006 (extracted via dx extract_dataset)
#
# Output:
#   fit_pca_iids.txt   (auto-uploaded by SAK)
#   pca_eur_log.txt    (auto-uploaded by SAK)
#   scrap/*            (auto-uploaded by SAK)

import csv
import os
import subprocess
from collections import defaultdict

EUROPEANS_PATH = "/mnt/project/sbayesrc_genotypes/europeans/classified_european_iids.txt"
CLOSE_RELATIONS_PATH = "/mnt/project/sbayesrc_genotypes/kinship/close_relations.csv"
KIN0_PATH = "/mnt/project/sbayesrc_genotypes/kinship/ukb_all_direct_rel.kin0"
IMPUTED_ONLY_PATH = "/mnt/project/sbayesrc_genotypes/merge_steps/imputed_only_iids.txt"
DIRECT_BFILE_PREFIX = "/mnt/project/sbayesrc_genotypes/direct_bfile/chr1_22_merged"

KINSHIP_THRESHOLD = 0.0441941  # 0.5^(9/2), lower bound for 3rd-degree relatives

DATASET_NAME_GLOB = "app*.dataset"

OUTPUT_FILE = "fit_pca_iids.txt"
LOG_FILE = "pca_eur_log.txt"

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
log(f"European IIDs (classified): {len(european_iids)}")

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

# ── 5. Identify WB siblings and expand to all relatives ────────────────────
log("")
log("=== Sample filtering ===")

wb_siblings = siblings & white_british
log(f"Step 1 — WB siblings (siblings & White British): {len(wb_siblings)}")

# Expand WB siblings to include all their relatives
expanded_exclusions = set(wb_siblings)
for iid in wb_siblings:
    expanded_exclusions.update(kinship_graph.get(iid, set()))
log(f"Step 2 — Expanded exclusions (WB siblings + all relatives): {len(expanded_exclusions)}")
log(f"  Added {len(expanded_exclusions) - len(wb_siblings)} relatives via kinship expansion")

# ── 6. Remove expanded exclusions from Europeans ───────────────────────────
excluded_europeans = european_iids & expanded_exclusions
pca_europeans = european_iids - expanded_exclusions
log(f"Step 3 — PCA Europeans after exclusion: {len(pca_europeans)}")
log(f"  Removed {len(excluded_europeans)} Europeans who were in expanded exclusion set")

# ── 7. Remove imputed-only IIDs ────────────────────────────────────────────
imputed_only = set()
with open(IMPUTED_ONLY_PATH) as f:
    for line in f:
        parts = line.strip().split()
        if parts:
            imputed_only.add(parts[0])
log(f"Imputed-only IIDs: {len(imputed_only)}")

removed_imputed = pca_europeans & imputed_only
pca_europeans = pca_europeans - imputed_only
log(f"Step 4 — PCA Europeans after removing imputed-only: {len(pca_europeans)}")
log(f"  Removed {len(removed_imputed)} imputed-only individuals")

# ── 8. Write intermediate IID file for plink2 ──────────────────────────────
os.makedirs("scrap", exist_ok=True)
pca_europeans_file = "scrap/pca_europeans_iids.txt"
write_iid_file(pca_europeans_file, pca_europeans)
log(f"\nWrote {pca_europeans_file} ({len(pca_europeans)} IIDs)")

# ── 9. Apply plink2 --king-cutoff-table for maximal unrelated set ──────────
log("")
log("=== KING cutoff (maximal unrelated set) ===")

plink2_cmd = [
    "plink2",
    "--bfile", DIRECT_BFILE_PREFIX,
    "--keep", pca_europeans_file,
    "--king-cutoff-table", KIN0_PATH, str(KINSHIP_THRESHOLD),
    "--out", "scrap/fit_pca",
]
log(f"Running: {' '.join(plink2_cmd)}")
subprocess.run(plink2_cmd, check=True)

# ── 10. Read plink2 output and write final IID file ────────────────────────
fit_pca_iids = set()
with open("scrap/fit_pca.king.cutoff.in.id") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        fit_pca_iids.add(parts[1])  # IID column

removed_by_king = pca_europeans - fit_pca_iids
log(f"Step 5 — Maximal unrelated set: {len(fit_pca_iids)}")
log(f"  Removed {len(removed_by_king)} related individuals via KING cutoff")

write_iid_file(OUTPUT_FILE, fit_pca_iids)

log("")
log("=== Final result ===")
log(f"fit_pca_iids: {len(fit_pca_iids)}")
log(f"\nWrote {OUTPUT_FILE} ({len(fit_pca_iids)} IIDs)")

# ── 11. Write log ──────────────────────────────────────────────────────────
with open(LOG_FILE, "w") as f:
    f.write("\n".join(log_lines) + "\n")

# ── 12. Verification ───────────────────────────────────────────────────────
print("\n=== Verification checks ===")

all_passed = True

# Check 1: Every fit_pca IID is European
not_eur = fit_pca_iids - european_iids
passed = len(not_eur) == 0
all_passed &= passed
print(f"Verification — european: {'PASS' if passed else 'FAIL'}")
print(f"  {len(fit_pca_iids)} IIDs checked, {len(not_eur)} not European")
if not_eur:
    print(f"  First 10: {sorted(list(not_eur))[:10]}")

# Check 2: No fit_pca IID is in expanded exclusions
in_excluded = fit_pca_iids & expanded_exclusions
passed = len(in_excluded) == 0
all_passed &= passed
print(f"Verification — not_in_exclusions: {'PASS' if passed else 'FAIL'}")
print(f"  {len(fit_pca_iids)} IIDs checked, {len(in_excluded)} in expanded exclusion set")
if in_excluded:
    print(f"  First 10: {sorted(list(in_excluded))[:10]}")

# Check 3: No fit_pca IID is imputed-only
in_imputed = fit_pca_iids & imputed_only
passed = len(in_imputed) == 0
all_passed &= passed
print(f"Verification — not_imputed_only: {'PASS' if passed else 'FAIL'}")
print(f"  {len(fit_pca_iids)} IIDs checked, {len(in_imputed)} in imputed-only set")
if in_imputed:
    print(f"  First 10: {sorted(list(in_imputed))[:10]}")

# Check 4: No pair of fit_pca IIDs has kinship >= threshold
violations = []
for iid in fit_pca_iids:
    related_in_set = kinship_graph.get(iid, set()) & fit_pca_iids
    if related_in_set:
        violations.append((iid, len(related_in_set)))
passed = len(violations) == 0
all_passed &= passed
print(f"Verification — unrelated: {'PASS' if passed else 'FAIL'}")
print(f"  {len(fit_pca_iids)} IIDs checked, {len(violations)} have relatives in set")
if violations:
    print(f"  First 10: {violations[:10]}")

if not all_passed:
    print("\nWARNING: Some verification checks failed!")

# ── 13. Cleanup ────────────────────────────────────────────────────────────
os.remove("white_british.csv")
os.remove("select_pca_europeans.py")
