#!/bin/bash
# get_genotypes.sh — Orchestrate genotype extraction pipeline.
#
# Steps:
#   1. Generate per-chromosome DRAGEN variant ID files (local)
#   2. Upload DRAGEN ID files to DNAnexus
#   3. Extract matching WGS variants from UKBB PLINK files (DNAnexus)
#   4. Generate per-chromosome TopMed variant ID files (local)
#   5. Upload TopMed ID files to DNAnexus
#   6. Extract matching imputed variants from UKBB BGEN files (DNAnexus)
#   7. Back up original pvar files on DNAnexus
#   8. Standardize pvar files (rsid mapping + column trimming)
#   9. Find IIDs in imputed data but not WGS (for merge)
#  10. Merge WGS + imputed-only individuals into per-chromosome bfiles
#  11. QC validation of merged pfiles (sample + variant count checks)
#  12. Convert merged pfiles to BGEN format + index
#  13. Extract direct SNPs from per-chromosome pfiles
#  14. Merge per-chromosome direct-SNP pfiles into one bfile
#  15. Subset direct SNPs to kinship-relevant SNPs (local + upload)
#  16. Run KING kinship estimation on DNAnexus
#  17. QC kinship: compare KING results against UKB reference
#  18. Classify close relationships from kinship data
#  19. Prepare ADMIXTURE projection inputs (download freqs, align alleles, build .P)
#  20. Split aligned bfile into batches for ADMIXTURE
#  21. Run ADMIXTURE K=6 projection per batch + concatenate into results TSV
#  22. Classify European ancestry individuals from ADMIXTURE results
#  23. Build train/test sample split from European ancestry + kinship data
#  24. Select unrelated European IIDs for fitting PCA
#  25. QC SNPs for PCA: MAF filter, LD region exclusion, LD pruning
#  26. Fit PCA on unrelated Europeans, project onto all samples
#  27. Build genetic sex covariate file (aneuploidy exclusion + sex assignment)
#
# Idempotent: each step checks for existing outputs and skips if already done.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
export GENO=0.03        # for WGS DRAGEN extraction: missingness threshold
export MAC=1000         # for WGS DRAGEN extraction: minor allele count threshold
export INSTANCE_TYPE="mem2_ssd1_v2_x16"  # DNAnexus instance for plink2 jobs
export DX_PRIORITY="high"               # DNAnexus job priority (low, normal, high)

# Local paths
export LOCAL_DRAGEN_ID_DIR="${SCRIPT_DIR}/data/dragen_ids"
export LOCAL_TOPMED_ID_DIR="${SCRIPT_DIR}/data/topmed_ids"
export ALIGNMENT_FILE="${SCRIPT_DIR}/data/support/sbayesrc_hg38.csv"

# DNAnexus output directory — all pipeline output lives here.
# Change this to write output to a different location. Must start with "/".
# Inside SAK job commands, files are accessed at /mnt/project${DX_OUTPUT_DIR}/...
export DX_OUTPUT_DIR="/sbayesrc_genotypes"

# DNAnexus paths (all derived from DX_OUTPUT_DIR)
export DX_DRAGEN_ID_DIR="${DX_OUTPUT_DIR}/dragen_ids"
export DX_WGS_PFILE_DIR="${DX_OUTPUT_DIR}/wgs_pfiles"
export WGS_PLINK_DIR="/mnt/project/Bulk/DRAGEN WGS/DRAGEN population level WGS variants, PLINK format [500k release]"
export DX_TOPMED_ID_DIR="${DX_OUTPUT_DIR}/topmed_ids"
export DX_IMPUTED_PFILE_DIR="${DX_OUTPUT_DIR}/imputed_pfiles"
export DX_BACKUP_DIR="${DX_OUTPUT_DIR}/backups"
export DX_MERGE_DIR="${DX_OUTPUT_DIR}/merge_steps"
export DX_MERGED_PFILE_DIR="${DX_OUTPUT_DIR}/merged_pfiles"
export DX_MERGED_BGEN_DIR="${DX_OUTPUT_DIR}/merged_bgens"
export MERGE_INSTANCE_TYPE="mem2_ssd1_v2_x16"  # merge jobs; increase for larger chromosomes
export IMPUTED_BGEN_DIR="/mnt/project/Bulk/Imputation/Imputation from genotype (TOPmed)"
export DX_DIRECT_PFILE_DIR="${DX_OUTPUT_DIR}/direct_pfiles"
export DX_DIRECT_BFILE_DIR="${DX_OUTPUT_DIR}/direct_bfile"
export LOCAL_DIRECT_SNPS_FILE="${SCRIPT_DIR}/data/support/direct_snps/ukbb_500k_qc_pass_direct_snps.txt"
export DX_DIRECT_SNPS_FILE="${DX_OUTPUT_DIR}/ukbb_500k_qc_pass_direct_snps.txt"

# Kinship estimation
export LOCAL_KINSHIP_SNPS_FILE="${SCRIPT_DIR}/data/support/ukbb_500k_qc_pass_direct_kinship_subsetted_snps.txt"
export LOCAL_SNP_QC_FILE="${SCRIPT_DIR}/data/support/ukb_snp_qc.txt"
export DX_KINSHIP_DIR="${DX_OUTPUT_DIR}/kinship"
export DX_KINSHIP_SNPS_FILE="${DX_OUTPUT_DIR}/kinship/ukbb_500k_qc_pass_direct_kinship_subsetted_snps.txt"
export KINSHIP_INSTANCE_TYPE="mem2_ssd1_v2_x64"

# ADMIXTURE K=6 projection
export ADMIXTURE_TSV_URL="https://raw.githubusercontent.com/jesseICR/public-statgen/main/outputs/admixture-global-6/admixture_allele_freqs.tsv"
export ADMIXTURE_DOWNLOAD_URL="https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz"
export ADMIXTURE_K=6
export ADMIXTURE_BATCH_SIZE=20000
export DX_STATGEN_DIR="${DX_OUTPUT_DIR}/statgen"
export DX_ADMIXTURE_SCRAP_DIR="${DX_OUTPUT_DIR}/statgen/scrap"
export DX_ADMIXTURE_BATCH_DIR="${DX_OUTPUT_DIR}/statgen/scrap/batches"
export ADMIXTURE_INSTANCE_TYPE="mem2_ssd1_v2_x2"
export ADMIXTURE_PREP_INSTANCE_TYPE="mem2_ssd1_v2_x16"

# European ancestry classification
export DX_EUROPEANS_DIR="${DX_OUTPUT_DIR}/europeans"

# Train/test sample split
export DX_TRAIN_TEST_DIR="${DX_OUTPUT_DIR}/train_test"

# PCA European sample selection
export DX_PCA_EUR_DIR="${DX_OUTPUT_DIR}/pca_eur"

# Genetic sex covariate
export DX_GENETIC_SEX_DIR="${DX_OUTPUT_DIR}/genetic_sex"

# ---------------------------------------------------------------------------
# Setup: Install Python dependencies
# ---------------------------------------------------------------------------
echo "=== Setup: Python dependencies ==="
pip install -r "${SCRIPT_DIR}/requirements.txt" --quiet

# ---------------------------------------------------------------------------
# Setup: Download SBayesRC alignment file (if not already cached)
# ---------------------------------------------------------------------------
echo "=== Setup: SBayesRC alignment file ==="
if [[ -s "${ALIGNMENT_FILE}" ]]; then
    echo "  Already cached at ${ALIGNMENT_FILE} — skipping download"
else
    echo "  Downloading sbayesrc_hg38.csv ..."
    curl -fsSL -o "${ALIGNMENT_FILE}" \
        "https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/sbayesrc_hg38.csv"
    echo "  Downloaded ($(wc -l < "${ALIGNMENT_FILE}") lines)"
fi

# ---------------------------------------------------------------------------
# Step 1: Generate DRAGEN variant IDs (local)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 1: Generate DRAGEN variant IDs ==="

python3 "${SCRIPT_DIR}/store_dragen_ids.py" \
    --input-file "${ALIGNMENT_FILE}" \
    --output-dir "${LOCAL_DRAGEN_ID_DIR}"
echo "Done."

# ---------------------------------------------------------------------------
# Step 2: Upload DRAGEN IDs to DNAnexus
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Upload DRAGEN IDs to DNAnexus ==="
bash "${SCRIPT_DIR}/upload_dragen_ids.sh"

# ---------------------------------------------------------------------------
# Step 3: Extract WGS variants from PLINK files (DNAnexus)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 3: Extract WGS variants per chromosome ==="
bash "${SCRIPT_DIR}/wgs_extract_variants.sh"

# ---------------------------------------------------------------------------
# Step 4: Generate TopMed variant IDs (local)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 4: Generate TopMed variant IDs ==="

python3 "${SCRIPT_DIR}/store_topmed_ids.py" \
    --input-file "${ALIGNMENT_FILE}" \
    --output-dir "${LOCAL_TOPMED_ID_DIR}"
echo "Done."

# ---------------------------------------------------------------------------
# Step 5: Upload TopMed IDs to DNAnexus
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 5: Upload TopMed IDs to DNAnexus ==="
bash "${SCRIPT_DIR}/upload_topmed_ids.sh"

# ---------------------------------------------------------------------------
# Step 6: Extract imputed variants from BGEN files (DNAnexus)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 6: Extract imputed variants per chromosome ==="
bash "${SCRIPT_DIR}/imputed_extract_variants.sh"

# ---------------------------------------------------------------------------
# Step 7: Back up original pvar files on DNAnexus
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 7: Back up original pvar files ==="
bash "${SCRIPT_DIR}/backup_pvars.sh"

# ---------------------------------------------------------------------------
# Step 8: Standardize pvar files (rsid mapping + column trimming)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 8: Standardize pvar files ==="
bash "${SCRIPT_DIR}/standardize_pvars.sh"

# ---------------------------------------------------------------------------
# Step 9: Find IIDs in imputed data but not WGS (for merge)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 9: Find imputed-only IIDs ==="
bash "${SCRIPT_DIR}/find_imputed_only_iids.sh"

# ---------------------------------------------------------------------------
# Step 10: Merge WGS + imputed-only individuals into per-chromosome bfiles
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 10: Merge WGS + imputed genotypes ==="
bash "${SCRIPT_DIR}/merge_wgs_imputed.sh"

# ---------------------------------------------------------------------------
# Step 11: QC — Validate merged pfiles
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 11: QC validation of merged pfiles ==="
bash "${SCRIPT_DIR}/validate_merged_pfiles.sh"

# ---------------------------------------------------------------------------
# Step 12: Convert merged pfiles to BGEN format + index
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 12: Convert merged pfiles to BGEN ==="
bash "${SCRIPT_DIR}/convert_to_bgens.sh"

# ---------------------------------------------------------------------------
# Step 13: Extract direct SNPs from per-chromosome pfiles
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 13: Extract direct SNPs per chromosome ==="
bash "${SCRIPT_DIR}/extract_direct_snps.sh"

# ---------------------------------------------------------------------------
# Step 14: Merge per-chromosome direct-SNP pfiles into one bfile
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 14: Merge direct-SNP pfiles into bfile ==="
bash "${SCRIPT_DIR}/make_direct_bfile.sh"

# ---------------------------------------------------------------------------
# Step 15: Subset direct SNPs to kinship-relevant SNPs
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 15: Subset direct SNPs to kinship-relevant SNPs ==="
bash "${SCRIPT_DIR}/subset_kinship_snps.sh"

# ---------------------------------------------------------------------------
# Step 16: Run KING kinship estimation on DNAnexus
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 16: Run KING kinship estimation ==="
bash "${SCRIPT_DIR}/run_king_kinship.sh"

# ---------------------------------------------------------------------------
# Step 17: QC kinship — compare KING results against UKB reference
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 17: Kinship QC ==="
bash "${SCRIPT_DIR}/kinship_qc.sh"

# ---------------------------------------------------------------------------
# Step 18: Classify close relationships from kinship data
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 18: Classify close relationships ==="
bash "${SCRIPT_DIR}/classify_relations.sh"

# ---------------------------------------------------------------------------
# Step 19: Prepare ADMIXTURE projection inputs
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 19: Prepare ADMIXTURE inputs ==="
bash "${SCRIPT_DIR}/admixture_prep.sh"

# ---------------------------------------------------------------------------
# Step 20: Split into ADMIXTURE batches
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 20: Split ADMIXTURE batches ==="
bash "${SCRIPT_DIR}/admixture_split_batches.sh"

# ---------------------------------------------------------------------------
# Step 21: Run ADMIXTURE projection and build results TSV
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 21: Run ADMIXTURE projection ==="
bash "${SCRIPT_DIR}/admixture_run_projection.sh"

# ---------------------------------------------------------------------------
# Step 22: Classify European ancestry individuals
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 22: Classify European ancestry individuals ==="
bash "${SCRIPT_DIR}/classify_europeans.sh"

# ---------------------------------------------------------------------------
# Step 23: Build train/test sample split
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 23: Build train/test sample split ==="
bash "${SCRIPT_DIR}/make_train_test_samples.sh"

# ---------------------------------------------------------------------------
# Step 24: Select unrelated European IIDs for fitting PCA
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 24: Select PCA European IIDs ==="
bash "${SCRIPT_DIR}/select_pca_europeans.sh"

# ---------------------------------------------------------------------------
# Step 25: QC SNPs for PCA (MAF, LD region exclusion, LD pruning)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 25: PCA SNP QC ==="
bash "${SCRIPT_DIR}/pca_snp_qc.sh"

# ---------------------------------------------------------------------------
# Step 26: Fit PCA on unrelated Europeans, project onto all samples
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 26: Fit PCA & project onto all samples ==="
bash "${SCRIPT_DIR}/fit_project_pca.sh"

# ---------------------------------------------------------------------------
# Step 27: Build genetic sex covariate file
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 27: Build genetic sex covariate ==="
bash "${SCRIPT_DIR}/get_genetic_sex.sh"

echo ""
echo "=== Pipeline complete ==="
