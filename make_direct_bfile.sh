#!/bin/bash
# make_direct_bfile.sh — Merge per-chromosome direct-SNP pfiles into one bfile.
#
# Submits a Swiss Army Knife job that uses plink2 --pmerge-list to merge
# chr1..chr22 direct-SNP pfiles into a single bfile.
#
# Output in direct_bfile/:
#   chr1_22_merged.bed, chr1_22_merged.bim, chr1_22_merged.fam, chr1_22_merged.log
#
# Expects env vars: DX_DIRECT_PFILE_DIR, DX_DIRECT_BFILE_DIR,
#                   INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

# Idempotency: skip if bfile already exists
if dx ls "${DX_DIRECT_BFILE_DIR}/chr1_22_merged.bed" &>/dev/null; then
    echo "  chr1_22_merged.bed already exists — skipping"
    exit 0
fi

dx mkdir -p "${DX_DIRECT_BFILE_DIR}"

input_dir="/mnt/project${DX_DIRECT_PFILE_DIR}"

cmd="set -eo pipefail && \
echo '--- Creating merge list ---' && \
for i in \$(seq 2 22); do echo '${input_dir}/chr'\${i}; done > merge_list.txt && \
cat merge_list.txt && \
echo '--- Merging chr1..22 direct-SNP pfiles into bfile ---' && \
plink2 --pfile '${input_dir}/chr1' --pmerge-list merge_list.txt pfile --make-bed --out chr1_22_merged && \
echo '--- Cleanup ---' && \
rm -f merge_list.txt chr1_22_merged.pgen chr1_22_merged.pvar chr1_22_merged.psam && \
echo 'Done: chr1_22_merged.bed/.bim/.fam created'"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_DIRECT_BFILE_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "make_direct_bfile" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"
echo "  Done."
