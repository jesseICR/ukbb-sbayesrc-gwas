#!/bin/bash
# setup_height_gwas.sh — Set up height GWAS example.
#
# Submits a Swiss Army Knife job that:
#   1. Extracts fields 50 (height) and 21003 (age) across instances 0-3
#   2. Intersects training IIDs with sex covariate (excludes aneuploidy)
#   3. Pairs height + age by instance, filters height < 140 cm
#   4. Collapses to median height and mean age per individual
#   5. Centers covariates (age_c, sex_c, interaction)
#   6. Merges with PC scores 1-10
#
# Output:       ${DX_HEIGHT_GWAS_DIR}/training_iids.txt
#               ${DX_HEIGHT_GWAS_DIR}/phen.txt
#               ${DX_HEIGHT_GWAS_DIR}/base_covar.txt
#               ${DX_HEIGHT_GWAS_DIR}/covar.txt
#               ${DX_HEIGHT_GWAS_DIR}/height_gwas_log.txt
# Verification: printed to stdout (visible in job log via `dx watch`)
#
# Expects env vars: DX_HEIGHT_GWAS_DIR, DX_OUTPUT_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Idempotency: skip if final output file already exists
if dx ls "${DX_HEIGHT_GWAS_DIR}/covar.txt" &>/dev/null; then
    echo "  covar.txt already exists — skipping"
    exit 0
fi

# Create output directory
dx mkdir -p "${DX_HEIGHT_GWAS_DIR}"

# Upload Python script as SAK input
script_id=$(dx upload "${SCRIPT_DIR}/setup_height_gwas.py" \
    --destination "${DX_HEIGHT_GWAS_DIR}/" --brief --no-progress)

cmd="export DX_OUTPUT_DIR='${DX_OUTPUT_DIR}' && set -eo pipefail && \
echo '--- Setting up height GWAS example ---' && \
python3 setup_height_gwas.py && \
echo '--- Done ---'"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_HEIGHT_GWAS_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "setup_height_gwas" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Clean up the uploaded Python script from DNAnexus
dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
