#!/bin/bash
# get_genetic_sex.sh — Build genetic sex covariate file.
#
# Submits a Swiss Army Knife job that:
#   1. Extracts fields 22001 (genetic sex) and 22019 (sex chr aneuploidy)
#   2. Excludes individuals with sex chromosome aneuploidy
#   3. Assigns sex from field 22001 for imputed-only IIDs
#   4. Assigns sex from fam file for WGS IIDs
#   5. Writes sex_covar.txt (FID IID sex_01) with 0=female, 1=male
#
# Output:       ${DX_GENETIC_SEX_DIR}/sex_covar.txt
#               ${DX_GENETIC_SEX_DIR}/readme.txt
#               ${DX_GENETIC_SEX_DIR}/genetic_sex_log.txt
# Verification: printed to stdout (visible in job log via `dx watch`)
#
# Expects env vars: DX_GENETIC_SEX_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Idempotency: skip if output file already exists
if dx ls "${DX_GENETIC_SEX_DIR}/sex_covar.txt" &>/dev/null; then
    echo "  sex_covar.txt already exists — skipping"
    exit 0
fi

# Create output directory
dx mkdir -p "${DX_GENETIC_SEX_DIR}"

# Upload Python script as SAK input
script_id=$(dx upload "${SCRIPT_DIR}/get_genetic_sex.py" \
    --destination "${DX_GENETIC_SEX_DIR}/" --brief --no-progress)

cmd="export DX_OUTPUT_DIR='${DX_OUTPUT_DIR}' && set -eo pipefail && \
echo '--- Building genetic sex covariate file ---' && \
python3 get_genetic_sex.py && \
echo '--- Done ---'"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_GENETIC_SEX_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "get_genetic_sex" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Clean up the uploaded Python script from DNAnexus
dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
