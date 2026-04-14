#!/bin/bash
# classify_relations.sh — Classify close relationships from KING kinship output.
#
# Submits a Swiss Army Knife job that:
#   1. Extracts birth year/month from UKBB tabular data (fields 34 + 52)
#   2. Classifies close relationships from KING .kin0 output
#   3. Joins with birth data, adds fractional year-of-birth columns
#   4. Filters out pairs with large age gap + low IBS0
#
# Output:       ${DX_KINSHIP_DIR}/close_relations.csv
# Intermediate: ${DX_KINSHIP_DIR}/yob_calc/birth_year_month.csv
#
# Expects env vars: DX_KINSHIP_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Idempotency: skip if output already exists
if dx ls "${DX_KINSHIP_DIR}/close_relations.csv" &>/dev/null; then
    echo "  close_relations.csv already exists — skipping"
    exit 0
fi

# Upload Python script as SAK input (use file ID to avoid path ambiguity)
script_id=$(dx upload "${SCRIPT_DIR}/classify_relations.py" \
    --destination "${DX_KINSHIP_DIR}/" --brief --no-progress)

cmd="export DX_OUTPUT_DIR='${DX_OUTPUT_DIR}' && set -eo pipefail && \
echo '--- Classifying and filtering close relationships ---' && \
python3 classify_relations.py && \
echo '--- Saving intermediate files ---' && \
dx mkdir -p '${DX_KINSHIP_DIR}/yob_calc/' && \
dx upload birth_year_month.csv --destination '${DX_KINSHIP_DIR}/yob_calc/' --brief --no-progress && \
rm -f birth_year_month.csv classify_relations.py"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_KINSHIP_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "classify_close_relations" \
    --ignore-reuse \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Clean up the uploaded Python script from DNAnexus
dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
