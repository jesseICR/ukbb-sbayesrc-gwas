#!/bin/bash
# run_king_kinship.sh — Run KING kinship estimation on DNAnexus.
#
# Submits a Swiss Army Knife job that runs plink2 --make-king-table on the
# direct bfile, extracting only the kinship-relevant SNP subset (~85k SNPs).
# KING estimates kinship coefficients for all sample pairs above the
# specified threshold.
#
# Output in kinship/:
#   ukb_all_direct_rel.kin0, ukb_all_direct_rel.log
#
# Expects env vars: DX_DIRECT_BFILE_DIR, DX_KINSHIP_DIR,
#                   DX_KINSHIP_SNPS_FILE, KINSHIP_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

# Idempotency: skip if .kin0 already exists
if dx ls "${DX_KINSHIP_DIR}/ukb_all_direct_rel.kin0" &>/dev/null; then
    echo "  ukb_all_direct_rel.kin0 already exists — skipping"
    exit 0
fi

dx mkdir -p "${DX_KINSHIP_DIR}"

bfile="/mnt/project${DX_DIRECT_BFILE_DIR}/chr1_22_merged"
snps="/mnt/project${DX_KINSHIP_SNPS_FILE}"

cmd="set -eo pipefail && \
echo '--- Running KING kinship estimation ---' && \
echo 'Input bfile: ${bfile}' && \
echo 'SNP extract list: ${snps}' && \
plink2 --bfile '${bfile}' \
    --extract '${snps}' \
    --make-king-table \
    --king-table-filter 0.03 \
    --out ukb_all_direct_rel && \
echo '--- Summary ---' && \
wc -l ukb_all_direct_rel.kin0 && \
echo 'Done: ukb_all_direct_rel.kin0 created'"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_KINSHIP_DIR}/" \
    --instance-type "${KINSHIP_INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "king_kinship_estimation" \
    --ignore-reuse \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"
echo "  Done."
