#!/bin/bash
# kinship_qc.sh — QC kinship: compare KING results against UKB reference.
#
# Submits a Swiss Army Knife job that compares our KING kinship output
# against UKB's ukb_rel.dat, producing summary statistics and diagnostic
# plots for both kinship coefficients and IBS0 values.
#
# Output: ${DX_KINSHIP_DIR}/qc/
#   - kinship_comparison_summary.txt, kinship_comparison_plots.png, kinship_bland_altman.png
#   - ibs0_comparison_summary.txt, ibs0_comparison_plots.png, ibs0_bland_altman.png
#
# Expects env vars: DX_KINSHIP_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Idempotency: skip if output already exists
if dx ls "${DX_KINSHIP_DIR}/qc/kinship_comparison_summary.txt" &>/dev/null; then
    echo "  kinship_comparison_summary.txt already exists — skipping"
    exit 0
fi

# Create output directory
dx mkdir -p "${DX_KINSHIP_DIR}/qc"

# Upload Python script as SAK input
script_id=$(dx upload "${SCRIPT_DIR}/kinship_qc.py" \
    --destination "${DX_KINSHIP_DIR}/qc/" --brief --no-progress)

cmd="export DX_OUTPUT_DIR='${DX_OUTPUT_DIR}' && set -eo pipefail && \
echo '--- Kinship QC: comparing WGS KING results against UKB ---' && \
pip install matplotlib > /dev/null 2>&1 || true && \
python3 kinship_qc.py && \
rm -f kinship_qc.py && \
echo '--- Done ---'"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_KINSHIP_DIR}/qc/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "kinship_qc" \
    --ignore-reuse \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Clean up the uploaded Python script from DNAnexus
dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
