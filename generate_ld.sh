#!/bin/bash
# generate_ld.sh — Optional pipeline, run AFTER get_genotypes.sh.
#
# Produces a QC-passed SNP list that will feed downstream LD-matrix generation.
# Each step is idempotent. Only Step 1 writes to DNAnexus; Steps 2–3 are local.
#
# Steps:
#   1. Compute WGS vs TopMed-imputed allele-frequency comparison (DNAnexus, EUR samples)
#   2. Download the per-variant frequency-comparison CSV from DNAnexus
#   3. QC-filter variants against three thresholds:
#        - MAF < MAF_THRESHOLD in WGS                                         (rare-variant noise)
#        - |alt_freq_wgs - alt_freq_topmed_imputed| > FREQ_DIFF_THRESHOLD     (imputation disagrees with WGS in UKB)
#        - |alt_freq_wgs - alt_freq_hrc_sbayesrc| > SBAYESRC_FREQ_DIFF_THRESHOLD (disagrees with SBayesRC's own reference panel)
#      Reference: Zheng et al. 2024, Nat Genet (SBayesRC paper).
#   4. Sample 40k LD-reference cohort (unrelated European White British, fixed seed)
#   5. Derive hg38 LD-block boundaries from SBayesRC's per-SNP Block assignments
#   6. Build per-chromosome WGS bfiles filtered to the 40k cohort x QC-passed SNPs
#   7. Initialize the LD-matrix workspace on DNAnexus (SBayesRC::LDstep1)
#   8. Build per-block LD matrices + eigen decomposition (LDstep2 + LDstep3, per-chrom)
#   9. Merge per-block LD info into snp.info (LDstep4)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# QC thresholds — modifiable
# ---------------------------------------------------------------------------
export MAF_THRESHOLD=0.009                  # min minor-allele frequency in WGS
export FREQ_DIFF_THRESHOLD=0.025            # max |WGS − TopMed-imputed| allele-freq diff
export SBAYESRC_FREQ_DIFF_THRESHOLD=0.03    # max |WGS − HRC-imputed (SBayesRC)| allele-freq diff

# ---------------------------------------------------------------------------
# DNAnexus config
# ---------------------------------------------------------------------------
export DX_PRIORITY="${DX_PRIORITY:-normal}"

export DX_OUTPUT_DIR="/sbayesrc_genotypes"
export DX_FREQ_COMPARE_DIR="${DX_OUTPUT_DIR}/freq_compare"
export DX_FREQ_COMPARE_CSV="${DX_FREQ_COMPARE_DIR}/wgs_vs_imputed_freq.csv"

# Source WGS pfiles (produced by get_genotypes.sh Step 3 + standardized in Step 8).
export DX_WGS_PFILE_DIR="${DX_OUTPUT_DIR}/wgs_pfiles"

# LD-reference outputs
export DX_LD_REF_DIR="${DX_OUTPUT_DIR}/ld_reference"
export DX_LD_COHORT_FILE="${DX_LD_REF_DIR}/ld_ref_40k_iids.txt"
export DX_LD_SNPS_FILE="${DX_LD_REF_DIR}/ld_ref_snps.txt"
export DX_LD_BFILE_DIR="${DX_LD_REF_DIR}/bfiles"
export DX_LD_BLOCKS_FILE="${DX_LD_REF_DIR}/ref4cM_hg38.pos"

# LD-matrix outputs (SBayesRC LDstep1–4)
export DX_LDM_DIR="${DX_LD_REF_DIR}/ld_matrices"
export DX_LDM_MA_FILE="${DX_LDM_DIR}/dummy_ma.ma"
export DX_LDM_INFO_FILE="${DX_LDM_DIR}/ldm.info"
export DX_LDM_SNP_INFO_FILE="${DX_LDM_DIR}/snp.info"

# ---------------------------------------------------------------------------
# Local paths
# ---------------------------------------------------------------------------
export LOCAL_FREQ_COMPARE_DIR="${SCRIPT_DIR}/data/freq_compare"
export LOCAL_FREQ_COMPARE_CSV="${LOCAL_FREQ_COMPARE_DIR}/wgs_vs_imputed_freq.csv"
export LOCAL_SBAYESRC_LIFTOVER_CSV="${LOCAL_FREQ_COMPARE_DIR}/sbayesrc_liftover_results.csv"
export LOCAL_QC_PASSED_CSV="${LOCAL_FREQ_COMPARE_DIR}/wgs_vs_imputed_freq_qc_passed.csv"

export LOCAL_LD_REF_DIR="${SCRIPT_DIR}/data/ld_reference"
export LOCAL_BLOCKS_HG38="${LOCAL_LD_REF_DIR}/ref4cM_hg38.pos"

# Canonical hg38 alignment file (downloaded by get_genotypes.sh).
# Used by Step 5's round-trip QC as the position source.
export ALIGNMENT_FILE="${SCRIPT_DIR}/data/support/sbayesrc_hg38.csv"

export SBAYESRC_LIFTOVER_URL="https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/sbayesrc_liftover_results.csv"

# ---------------------------------------------------------------------------
# LD-reference cohort + bfile config
# ---------------------------------------------------------------------------
export LD_COHORT_SIZE=40000
export RANDOM_SEED=0
export LD_REF_INSTANCE_TYPE="mem2_ssd1_v2_x16"

# LD-matrix generation (Steps 7–9)
export LD_MATRICES_INSTANCE_TYPE="mem2_ssd1_v2_x16"
export LD_MATRICES_OMP_THREADS=16

# ---------------------------------------------------------------------------
# Logging — all subsequent output goes to both terminal and log file
# ---------------------------------------------------------------------------
mkdir -p "${SCRIPT_DIR}/logs"
LOG_FILE="${SCRIPT_DIR}/logs/generate_ld_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "Generate-LD Pipeline — $(date)"
echo "============================================"
echo "  LOG_FILE=${LOG_FILE}"
echo "  MAF_THRESHOLD=${MAF_THRESHOLD}"
echo "  FREQ_DIFF_THRESHOLD=${FREQ_DIFF_THRESHOLD}"
echo "  SBAYESRC_FREQ_DIFF_THRESHOLD=${SBAYESRC_FREQ_DIFF_THRESHOLD}"
echo ""

# ---------------------------------------------------------------------------
# Setup: Install Python dependencies
# ---------------------------------------------------------------------------
echo "=== Setup: Python dependencies ==="
pip install -r "${SCRIPT_DIR}/requirements.txt" --quiet

# ---------------------------------------------------------------------------
# Step 1: Compute WGS vs TopMed-imputed allele-frequency comparison (DNAnexus)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 1: Compute WGS vs TopMed allele-frequency comparison ==="
bash "${SCRIPT_DIR}/compute_freq_compare.sh"

# ---------------------------------------------------------------------------
# Step 2: Download the frequency-comparison CSV from DNAnexus
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 2: Download frequency-comparison CSV ==="
bash "${SCRIPT_DIR}/download_freq_compare.sh"

# ---------------------------------------------------------------------------
# Step 3: QC-filter variants (MAF + TopMed-freq-diff + HRC-freq-diff)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 3: QC-filter variants ==="
bash "${SCRIPT_DIR}/qc_snps.sh"

# ---------------------------------------------------------------------------
# Step 4: Sample 40k LD-reference cohort (DNAnexus)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 4: Sample 40k LD-reference cohort ==="
bash "${SCRIPT_DIR}/sample_ld_cohort.sh"

# ---------------------------------------------------------------------------
# Step 5: Derive hg38 LD-block boundaries from per-SNP Block assignments
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 5: Derive hg38 LD-block boundaries ==="
bash "${SCRIPT_DIR}/build_hg38_blocks.sh"

# ---------------------------------------------------------------------------
# Step 6: Build per-chromosome LD-reference bfiles (DNAnexus)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 6: Build per-chromosome LD-reference bfiles ==="
bash "${SCRIPT_DIR}/build_ld_ref_bfiles.sh"

# ---------------------------------------------------------------------------
# Step 7: Initialize LD-matrix workspace (SBayesRC::LDstep1)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 7: Initialize LD-matrix workspace (LDstep1) ==="
bash "${SCRIPT_DIR}/init_ld_matrices.sh"

# ---------------------------------------------------------------------------
# Step 8: Build per-chromosome LD matrices + eigen (LDstep2 + LDstep3)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 8: Build LD matrices + eigen decomposition (LDstep2 + LDstep3) ==="
bash "${SCRIPT_DIR}/build_ld_matrices.sh"

# ---------------------------------------------------------------------------
# Step 9: Merge per-block LD info (SBayesRC::LDstep4)
# ---------------------------------------------------------------------------
echo ""
echo "=== Step 9: Merge per-block LD info (LDstep4) ==="
bash "${SCRIPT_DIR}/merge_ld_info.sh"

echo ""
echo "=== Generate-LD pipeline complete ==="
echo "    Final QC-passed CSV: ${LOCAL_QC_PASSED_CSV}"
echo "    LD-reference cohort: ${DX_LD_COHORT_FILE}"
echo "    LD-reference blocks: ${LOCAL_BLOCKS_HG38} (also at ${DX_LD_BLOCKS_FILE})"
echo "    LD-reference bfiles: ${DX_LD_BFILE_DIR}/chr{1..22}.{bed,bim,fam,log}"
echo "    LD matrices (eigen): ${DX_LDM_DIR}/block{i}.eigen.bin"
echo "    LD matrices (info):  ${DX_LDM_INFO_FILE} and ${DX_LDM_SNP_INFO_FILE}"
