#!/bin/bash
# find_imputed_only_iids.sh — Identify positive IIDs for merge: WGS and imputed-only.
#
# Submits a Swiss Army Knife job that compares chr1 psam files to find:
#   1. Positive (non-negative) WGS IIDs → wgs_positive_iids.txt
#   2. Positive IIDs in imputed but not WGS → imputed_only_iids.txt
# Both are two-column FID/IID files for use with plink2 --keep.
#
# Expects env vars: DX_WGS_PFILE_DIR, DX_IMPUTED_PFILE_DIR, DX_MERGE_DIR,
#                   INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

# Idempotency: skip if output already exists
if dx ls "${DX_MERGE_DIR}/imputed_only_iids.txt" &>/dev/null && \
   dx ls "${DX_MERGE_DIR}/wgs_positive_iids.txt" &>/dev/null; then
    echo "imputed_only_iids.txt and wgs_positive_iids.txt already exist — skipping"
    exit 0
fi

dx mkdir -p "${DX_MERGE_DIR}"

wgs_psam="/mnt/project/${DX_WGS_PFILE_DIR}/chr1.psam"
imputed_psam="/mnt/project/${DX_IMPUTED_PFILE_DIR}/chr1.psam"

cmd="tail -n +2 '${wgs_psam}' | awk '\$2 > 0 {print \$2}' | sort > wgs_pos.tmp && \
tail -n +2 '${imputed_psam}' | awk '\$2 > 0 {print \$2}' | sort > imp_pos.tmp && \
comm -23 imp_pos.tmp wgs_pos.tmp | awk '{print \$1, \$1}' > imputed_only_iids.txt && \
awk '{print \$1, \$1}' wgs_pos.tmp > wgs_positive_iids.txt && \
rm -f wgs_pos.tmp imp_pos.tmp"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_MERGE_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "find_imputed_only_iids" \
    -y --brief)

echo "Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"
echo "Done."
