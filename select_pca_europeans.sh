#!/bin/bash
# select_pca_europeans.sh — Select unrelated European IIDs for fitting PCA.
#
# Submits a Swiss Army Knife job that:
#   1. Identifies WB siblings and expands to all their relatives
#   2. Removes expanded exclusions and imputed-only IIDs from Europeans
#   3. Applies plink2 --king-cutoff-table for maximal unrelated set
#
# Output:       ${DX_PCA_EUR_DIR}/fit_pca_iids.txt
#               ${DX_PCA_EUR_DIR}/pca_eur_log.txt
# Verification: printed to stdout (visible in job log via `dx watch`)
#
# Expects env vars: DX_PCA_EUR_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Idempotency: skip if output file already exists
if dx ls "${DX_PCA_EUR_DIR}/fit_pca_iids.txt" &>/dev/null; then
    echo "  fit_pca_iids.txt already exists — skipping"
    exit 0
fi

# Create output directory
dx mkdir -p "${DX_PCA_EUR_DIR}"

# Upload Python script as SAK input
script_id=$(dx upload "${SCRIPT_DIR}/select_pca_europeans.py" \
    --destination "${DX_PCA_EUR_DIR}/" --brief --no-progress)

cmd="set -eo pipefail && \
echo '--- Selecting PCA European IIDs ---' && \
python3 select_pca_europeans.py && \
echo '--- Done ---'"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_PCA_EUR_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "select_pca_europeans" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Clean up the uploaded Python script from DNAnexus
dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
