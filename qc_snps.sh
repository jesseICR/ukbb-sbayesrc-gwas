#!/bin/bash
# qc_snps.sh — Filter the WGS vs TopMed frequency-comparison CSV down to a QC-
# passed SNP set suitable for LD-matrix generation.
#
# Three filters are applied:
#   (1) MAF >= MAF_THRESHOLD in WGS
#   (2) |alt_freq_wgs − alt_freq_topmed_imputed| <= FREQ_DIFF_THRESHOLD
#   (3) |alt_freq_wgs − alt_freq_hrc_sbayesrc|   <= SBAYESRC_FREQ_DIFF_THRESHOLD
#
# Filter (3) requires the SBayesRC liftover CSV (hg19→hg38 allele mapping with
# HRC-imputed white-British allele frequencies). Downloaded once and cached.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -s "${LOCAL_QC_PASSED_CSV}" ]]; then
    echo "Already exists at ${LOCAL_QC_PASSED_CSV} — skipping"
    exit 0
fi

mkdir -p "${LOCAL_FREQ_COMPARE_DIR}"

# --- Cache-download SBayesRC liftover CSV ---
if [[ -s "${LOCAL_SBAYESRC_LIFTOVER_CSV}" ]]; then
    echo "SBayesRC liftover CSV already cached at ${LOCAL_SBAYESRC_LIFTOVER_CSV}"
else
    echo "Downloading SBayesRC liftover CSV ..."
    curl -fsSL -o "${LOCAL_SBAYESRC_LIFTOVER_CSV}" "${SBAYESRC_LIFTOVER_URL}"
    echo "Downloaded ($(wc -l < "${LOCAL_SBAYESRC_LIFTOVER_CSV}") lines)"
fi

# --- Run Python QC ---
python3 "${SCRIPT_DIR}/qc_snps.py" \
    --freq-compare "${LOCAL_FREQ_COMPARE_CSV}" \
    --sbayesrc-liftover "${LOCAL_SBAYESRC_LIFTOVER_CSV}" \
    --maf-threshold "${MAF_THRESHOLD}" \
    --freq-diff-threshold "${FREQ_DIFF_THRESHOLD}" \
    --sbayesrc-freq-diff-threshold "${SBAYESRC_FREQ_DIFF_THRESHOLD}" \
    --output "${LOCAL_QC_PASSED_CSV}"
