#!/bin/bash
# sample_ld_cohort.sh - Sample LD-reference cohort for custom SBayesRC LD.
#
# Submits a Swiss Army Knife job that intersects three sources and random-
# samples LD_COHORT_SIZE individuals with a fixed seed:
#   1. fit_pca_iids.txt     (unrelated Europeans for PCA fitting)
#   2. sex_covar.txt        (IIDs with genetic-sex call; aneuploidy excluded)
#   3. UKBB field 22006     (White British)
#
# Output:       ${DX_LD_COHORT_FILE}
#               ${DX_LD_REF_DIR}/ld_ref_cohort_log.txt
# Verification: printed to stdout (visible via `dx watch`)
#
# Expects env vars: DX_OUTPUT_DIR, DX_LD_REF_DIR, DX_LD_COHORT_FILE,
#                   LD_COHORT_SIZE, RANDOM_SEED,
#                   LD_REF_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Idempotency: skip if output file already exists
if dx ls "${DX_LD_COHORT_FILE}" &>/dev/null; then
    echo "  ld_ref_40k_iids.txt already exists - skipping"
    exit 0
fi

dx mkdir -p "${DX_LD_REF_DIR}"

script_id=$(dx upload "${SCRIPT_DIR}/sample_ld_cohort.py" \
    --destination "${DX_LD_REF_DIR}/" --brief --no-progress)

cmd="export DX_OUTPUT_DIR='${DX_OUTPUT_DIR}' \
LD_COHORT_SIZE='${LD_COHORT_SIZE}' \
RANDOM_SEED='${RANDOM_SEED}' && \
set -eo pipefail && \
echo '--- Sampling LD-reference cohort ---' && \
python3 sample_ld_cohort.py && \
echo '--- Done ---'"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_LD_REF_DIR}/" \
    --instance-type "${LD_REF_INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "sample_ld_cohort" \
    --ignore-reuse \
    -y --brief)

echo "  Submitted ${job_id} - waiting for completion ..."
dx wait "${job_id}"

dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
