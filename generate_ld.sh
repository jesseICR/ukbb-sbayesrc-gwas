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

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# QC thresholds — modifiable
# ---------------------------------------------------------------------------
export MAF_THRESHOLD=0.007                  # min minor-allele frequency in WGS
export FREQ_DIFF_THRESHOLD=0.025            # max |WGS − TopMed-imputed| allele-freq diff
export SBAYESRC_FREQ_DIFF_THRESHOLD=0.03    # max |WGS − HRC-imputed (SBayesRC)| allele-freq diff

# ---------------------------------------------------------------------------
# DNAnexus config
# ---------------------------------------------------------------------------
export DX_PRIORITY="${DX_PRIORITY:-normal}"

export DX_OUTPUT_DIR="/sbayesrc_genotypes"
export DX_FREQ_COMPARE_DIR="${DX_OUTPUT_DIR}/freq_compare"
export DX_FREQ_COMPARE_CSV="${DX_FREQ_COMPARE_DIR}/wgs_vs_imputed_freq.csv"

# ---------------------------------------------------------------------------
# Local paths
# ---------------------------------------------------------------------------
export LOCAL_FREQ_COMPARE_DIR="${SCRIPT_DIR}/data/freq_compare"
export LOCAL_FREQ_COMPARE_CSV="${LOCAL_FREQ_COMPARE_DIR}/wgs_vs_imputed_freq.csv"
export LOCAL_SBAYESRC_LIFTOVER_CSV="${LOCAL_FREQ_COMPARE_DIR}/sbayesrc_liftover_results.csv"
export LOCAL_QC_PASSED_CSV="${LOCAL_FREQ_COMPARE_DIR}/wgs_vs_imputed_freq_qc_passed.csv"

export SBAYESRC_LIFTOVER_URL="https://github.com/jesseICR/sbayesrc-liftover/releases/download/v1.0/sbayesrc_liftover_results.csv"

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

echo ""
echo "=== Generate-LD pipeline complete ==="
echo "    Final QC-passed CSV: ${LOCAL_QC_PASSED_CSV}"
