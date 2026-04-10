#!/bin/bash
# make_train_test_samples.sh — Build train/test sample split for REGENIE.
#
# Submits a Swiss Army Knife job that:
#   1. Extracts White British classification from UKBB field 22006
#   2. Builds kinship graph from KING .kin0 output (3rd degree and above)
#   3. Splits European individuals into train and test samples
#   4. Runs verification checks on the final split
#
# Output:       ${DX_TRAIN_TEST_DIR}/final_train_iids.txt
#               ${DX_TRAIN_TEST_DIR}/final_test_iids.txt
#               ${DX_TRAIN_TEST_DIR}/train_test_log.txt
# Verification: printed to stdout (visible in job log via `dx watch`)
#
# Expects env vars: DX_TRAIN_TEST_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Idempotency: skip if both output files already exist
if dx ls "${DX_TRAIN_TEST_DIR}/final_train_iids.txt" &>/dev/null && \
   dx ls "${DX_TRAIN_TEST_DIR}/final_test_iids.txt" &>/dev/null; then
    echo "  final_train_iids.txt and final_test_iids.txt already exist — skipping"
    exit 0
fi

# Create output directory
dx mkdir -p "${DX_TRAIN_TEST_DIR}"

# Upload Python script as SAK input
script_id=$(dx upload "${SCRIPT_DIR}/make_train_test_samples.py" \
    --destination "${DX_TRAIN_TEST_DIR}/" --brief --no-progress)

cmd="export DX_OUTPUT_DIR='${DX_OUTPUT_DIR}' && set -eo pipefail && \
echo '--- Building train/test sample split ---' && \
python3 make_train_test_samples.py && \
echo '--- Done ---'"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_TRAIN_TEST_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "make_train_test_samples" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Clean up the uploaded Python script from DNAnexus
dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
