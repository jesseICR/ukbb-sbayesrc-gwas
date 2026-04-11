"""Set up height GWAS example: phenotype, covariates, and training sample files.

Algorithm:
  1. Intersect final_train_iids with sex_covar to get eligible IIDs
     (train IIDs not in sex_covar have sex chromosome aneuploidy)
  2. Extract UK Biobank fields 50 (height) and 21003 (age) across instances 0-3
  3. For each eligible IID: pair height+age by instance, drop height < 140 cm,
     compute median(height) and mean(age) across remaining measurements
  4. Center covariates: age_c = age - mean(age), sex_c = sex - 0.5,
     interaction = age_c * sex_c
  5. Merge with PC scores 1-10 from ukb_projected.sscore
  6. Write output files with validation checks

Input:
  /mnt/project${DX_OUTPUT_DIR}/train_test/final_train_iids.txt
  /mnt/project${DX_OUTPUT_DIR}/genetic_sex/sex_covar.txt
  /mnt/project${DX_OUTPUT_DIR}/pca_eur/ukb_projected.sscore
  UK Biobank fields: participant.p50_i0..i3, participant.p21003_i0..i3

Output (auto-uploaded by SAK):
  training_iids.txt      -- FID IID (no header, space-separated)
  phen.txt               -- FID IID height (tab-separated, with header)
  base_covar.txt         -- FID IID age_c sex_c age_c_sex_c_inter (tab-separated)
  covar.txt              -- base_covar + PC1_AVG..PC10_AVG (tab-separated)
  height_gwas_log.txt    -- step-by-step log with counts
"""

import csv
import os
import statistics
import subprocess
import sys

DATASET_NAME_GLOB = "app*.dataset"

DX_OUTPUT_DIR = os.environ.get("DX_OUTPUT_DIR", "/sbayesrc_genotypes")
_BASE = f"/mnt/project{DX_OUTPUT_DIR}"

TRAIN_IID_PATH = f"{_BASE}/train_test/final_train_iids.txt"
SEX_COVAR_PATH = f"{_BASE}/genetic_sex/sex_covar.txt"
SSCORE_PATH = f"{_BASE}/pca_eur/ukb_projected.sscore"

EXTRACT_PATH = "ukb_height_age.csv"
HEIGHT_THRESHOLD = 140  # cm -- measurements below this are excluded

TRAINING_IIDS_PATH = "training_iids.txt"
PHEN_PATH = "phen.txt"
BASE_COVAR_PATH = "base_covar.txt"
COVAR_PATH = "covar.txt"
LOG_PATH = "height_gwas_log.txt"

log_lines = []


def log(msg):
    print(msg)
    log_lines.append(msg)


# =============================================================================
# Step A: Determine eligible IIDs (train intersect sex_covar)
# =============================================================================

# Load training IIDs
train_iids = set()
with open(TRAIN_IID_PATH) as f:
    for line in f:
        parts = line.strip().split()
        if parts:
            train_iids.add(parts[1])

# Load sex covariate
sex_map = {}  # iid -> int (0=female, 1=male)
with open(SEX_COVAR_PATH) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        sex_map[row["IID"]] = int(row["sex_01"])

eligible_iids = train_iids & set(sex_map.keys())
excluded_aneuploidy = train_iids - eligible_iids

log("=== Step A: Eligible IIDs ===")
log(f"  Training IIDs: {len(train_iids)}")
log(f"  Sex covar IIDs: {len(sex_map)}")
log(f"  Eligible (train intersect sex_covar): {len(eligible_iids)}")
log(f"  Excluded (likely aneuploidy): {len(excluded_aneuploidy)}")

# =============================================================================
# Step B: Extract UK Biobank height and age fields
# =============================================================================

project_id = os.environ["DX_PROJECT_CONTEXT_ID"]
result = subprocess.run(
    ["dx", "find", "data", "--name", DATASET_NAME_GLOB,
     "--project", project_id, "--brief"],
    capture_output=True, text=True, check=True,
)
dataset_id = result.stdout.strip().split("\n")[0]
log(f"\n=== Step B: Extract UKBB fields ===")
log(f"  Discovered dataset: {dataset_id}")

height_fields = [f"participant.p50_i{i}" for i in range(4)]
age_fields = [f"participant.p21003_i{i}" for i in range(4)]
all_fields = ",".join(["participant.eid"] + height_fields + age_fields)

subprocess.run(
    ["dx", "extract_dataset", dataset_id,
     "--fields", all_fields,
     "--delimiter", ",",
     "--output", EXTRACT_PATH],
    check=True,
)
extract_rows = sum(1 for _ in open(EXTRACT_PATH)) - 1
log(f"  Extracted fields 50 + 21003 across instances 0-3 ({extract_rows} rows)")

# =============================================================================
# Step C: Process phenotype -- pair by instance, filter, collapse
# =============================================================================

log(f"\n=== Step C: Phenotype processing ===")

phenotype = {}        # iid -> (median_height, mean_age)
eligible_in_extract = set()
no_paired_data = set()
all_below_threshold = set()

with open(EXTRACT_PATH) as f:
    reader = csv.DictReader(f)
    for row in reader:
        eid = row["participant.eid"]
        if eid not in eligible_iids:
            continue
        eligible_in_extract.add(eid)

        # Collect (height, age) pairs where both are present at the same instance
        pairs = []
        for i in range(4):
            h_val = row.get(f"participant.p50_i{i}", "").strip()
            a_val = row.get(f"participant.p21003_i{i}", "").strip()
            if h_val and a_val:
                pairs.append((float(h_val), float(a_val)))

        if not pairs:
            no_paired_data.add(eid)
            continue

        # Filter out measurements with height < threshold
        valid_pairs = [(h, a) for h, a in pairs if h >= HEIGHT_THRESHOLD]

        if not valid_pairs:
            all_below_threshold.add(eid)
            continue

        heights = [h for h, a in valid_pairs]
        ages = [a for h, a in valid_pairs]
        phenotype[eid] = (statistics.median(heights), statistics.mean(ages))

not_in_extract = eligible_iids - eligible_in_extract

log(f"  Eligible IIDs found in extract: {len(eligible_in_extract)}")
log(f"  Eligible IIDs not in extract (no height/age data): {len(not_in_extract)}")
log(f"  Eligible IIDs with no paired height+age at any instance: {len(no_paired_data)}")
log(f"  Eligible IIDs with all heights < {HEIGHT_THRESHOLD} cm: {len(all_below_threshold)}")
log(f"  Final GWAS individuals with valid phenotype: {len(phenotype)}")

# =============================================================================
# Step D: Build covariates (center age and sex, compute interaction)
# =============================================================================

gwas_iids = sorted(phenotype.keys(), key=lambda x: int(x))

all_ages = [phenotype[iid][1] for iid in gwas_iids]
mean_age = statistics.mean(all_ages)

covar_data = {}  # iid -> (age_c, sex_c, age_c_sex_c_inter)
for iid in gwas_iids:
    _, ind_age = phenotype[iid]
    age_c = ind_age - mean_age
    sex_c = sex_map[iid] - 0.5
    inter = age_c * sex_c
    covar_data[iid] = (age_c, sex_c, inter)

log(f"\n=== Step D: Covariate centering ===")
log(f"  Mean age across {len(gwas_iids)} GWAS individuals: {mean_age:.4f}")
log(f"  Mean age_c (should be ~0): {statistics.mean([c[0] for c in covar_data.values()]):.6f}")

males = sum(1 for iid in gwas_iids if sex_map[iid] == 1)
females = sum(1 for iid in gwas_iids if sex_map[iid] == 0)
log(f"  Males: {males}, Females: {females}, Ratio: {males / females:.4f}")

# =============================================================================
# Step E: Merge with PC scores 1-10
# =============================================================================

pc_data = {}  # iid -> list of 10 PC value strings
gwas_set = set(gwas_iids)

with open(SSCORE_PATH) as f:
    header = f.readline().strip().split("\t")
    if header[0] == "#FID":
        header[0] = "FID"
    iid_idx = header.index("IID")
    pc_indices = [header.index(f"PC{i}_AVG") for i in range(1, 11)]

    for line in f:
        fields = line.strip().split("\t")
        iid = fields[iid_idx]
        if iid in gwas_set:
            pc_data[iid] = [fields[idx] for idx in pc_indices]

missing_pcs = gwas_set - set(pc_data.keys())

log(f"\n=== Step E: PC merge ===")
log(f"  GWAS IIDs with PC data: {len(pc_data)}")
if missing_pcs:
    log(f"  ERROR: {len(missing_pcs)} GWAS IIDs missing from sscore file")
    sys.exit(1)
else:
    log(f"  All GWAS IIDs have PC data")

# =============================================================================
# Write output files
# =============================================================================

log(f"\n=== Writing output files ===")

# training_iids.txt -- FID IID, no header, space-separated
with open(TRAINING_IIDS_PATH, "w") as f:
    for iid in gwas_iids:
        f.write(f"{iid} {iid}\n")
log(f"  Wrote {TRAINING_IIDS_PATH} ({len(gwas_iids)} rows)")

# phen.txt -- FID IID height, tab-separated with header
with open(PHEN_PATH, "w") as f:
    f.write("FID\tIID\theight\n")
    for iid in gwas_iids:
        height, _ = phenotype[iid]
        f.write(f"{iid}\t{iid}\t{height}\n")
log(f"  Wrote {PHEN_PATH} ({len(gwas_iids)} rows)")

# base_covar.txt -- FID IID age_c sex_c age_c_sex_c_inter, tab-separated
with open(BASE_COVAR_PATH, "w") as f:
    f.write("FID\tIID\tage_c\tsex_c\tage_c_sex_c_inter\n")
    for iid in gwas_iids:
        age_c, sex_c, inter = covar_data[iid]
        f.write(f"{iid}\t{iid}\t{age_c}\t{sex_c}\t{inter}\n")
log(f"  Wrote {BASE_COVAR_PATH} ({len(gwas_iids)} rows)")

# covar.txt -- base_covar columns + PC1_AVG..PC10_AVG, tab-separated
pc_headers = "\t".join(f"PC{i}_AVG" for i in range(1, 11))
with open(COVAR_PATH, "w") as f:
    f.write(f"FID\tIID\tage_c\tsex_c\tage_c_sex_c_inter\t{pc_headers}\n")
    for iid in gwas_iids:
        age_c, sex_c, inter = covar_data[iid]
        pcs = "\t".join(pc_data[iid])
        f.write(f"{iid}\t{iid}\t{age_c}\t{sex_c}\t{inter}\t{pcs}\n")
log(f"  Wrote {COVAR_PATH} ({len(gwas_iids)} rows)")

# =============================================================================
# Verification checks
# =============================================================================

print("\n=== Verification checks ===")
all_passed = True


def check(name, condition, detail=""):
    global all_passed
    if condition:
        print(f"PASS: {name}" + (f" ({detail})" if detail else ""))
    else:
        print(f"FAIL: {name}" + (f" ({detail})" if detail else ""))
        all_passed = False


# Check 1: Read back IIDs from all output files and compare
def read_iids(path, has_header, sep="\t"):
    iids = []
    with open(path) as f:
        if has_header:
            f.readline()
        for line in f:
            parts = line.strip().split(sep)
            iids.append(parts[1])
    return iids

iids_training = read_iids(TRAINING_IIDS_PATH, has_header=False, sep=" ")
iids_phen = read_iids(PHEN_PATH, has_header=True)
iids_base = read_iids(BASE_COVAR_PATH, has_header=True)
iids_covar = read_iids(COVAR_PATH, has_header=True)

check("Same IIDs in training_iids and phen", iids_training == iids_phen)
check("Same IIDs in training_iids and base_covar", iids_training == iids_base)
check("Same IIDs in training_iids and covar", iids_training == iids_covar)

# Check 2: All IIDs are in final_train_iids.txt
check("All GWAS IIDs in final_train_iids",
      set(iids_training) <= train_iids,
      f"{len(iids_training)} GWAS IIDs, all in training set")

# Check 3: Row counts match
n = len(iids_training)
check("Row counts match across files",
      len(iids_phen) == n and len(iids_base) == n and len(iids_covar) == n,
      f"all have {n} rows")

# Check 4: No missing values in phen.txt
missing_phen = 0
with open(PHEN_PATH) as f:
    f.readline()
    for line in f:
        if any(v.strip() == "" for v in line.strip().split("\t")):
            missing_phen += 1
check("No missing values in phen.txt", missing_phen == 0)

# Check 5: No missing values in covar.txt
missing_covar = 0
with open(COVAR_PATH) as f:
    f.readline()
    for line in f:
        if any(v.strip() == "" for v in line.strip().split("\t")):
            missing_covar += 1
check("No missing values in covar.txt", missing_covar == 0)

# Check 6: All heights > threshold
heights_ok = all(phenotype[iid][0] >= HEIGHT_THRESHOLD for iid in gwas_iids)
check(f"All heights >= {HEIGHT_THRESHOLD} cm", heights_ok)

# Check 7: Mean of age_c approximately 0
mean_age_c = statistics.mean([covar_data[iid][0] for iid in gwas_iids])
check("Mean of age_c ~ 0", abs(mean_age_c) < 0.001, f"mean = {mean_age_c:.8f}")

# Check 8: sex_c values are {-0.5, 0.5}
sex_c_vals = set(covar_data[iid][1] for iid in gwas_iids)
check("sex_c values are {-0.5, 0.5}", sex_c_vals == {-0.5, 0.5})

# Check 9: FID == IID in all files
fid_eq_iid = True
for path, has_header, sep in [
    (TRAINING_IIDS_PATH, False, " "),
    (PHEN_PATH, True, "\t"),
    (BASE_COVAR_PATH, True, "\t"),
    (COVAR_PATH, True, "\t"),
]:
    with open(path) as f:
        if has_header:
            f.readline()
        for line in f:
            parts = line.strip().split(sep)
            if parts[0] != parts[1]:
                fid_eq_iid = False
                break
check("FID == IID in all output files", fid_eq_iid)

if all_passed:
    print("\nAll verification checks passed.")
else:
    print("\nWARNING: Some verification checks failed!")
    sys.exit(1)

# =============================================================================
# Write log file
# =============================================================================

with open(LOG_PATH, "w") as f:
    f.write("\n".join(log_lines) + "\n")

# =============================================================================
# Cleanup (prevent SAK from auto-uploading intermediates)
# =============================================================================

os.remove(EXTRACT_PATH)
os.remove("setup_height_gwas.py")
