#!/bin/bash
# convert_to_bgens.sh — Convert merged pfiles to BGEN v1.2 format and index.
#
# For each chromosome, submits a Swiss Army Knife job that:
#   1. Converts the merged pfile to BGEN v1.2 (8-bit, ref-first) with plink2
#   2. Indexes the BGEN with bgenix
#
# Outputs per chromosome in merged_bgens/:
#   chr{N}.bgen, chr{N}.bgen.bgi, chr{N}.sample, chr{N}.log
#
# Expects env vars: DX_MERGED_PFILE_DIR, DX_MERGED_BGEN_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

CHROMS=$(seq 1 22)

echo "  Converting merged pfiles to BGEN + indexing ..."

dx mkdir -p "${DX_MERGED_BGEN_DIR}"

submitted=0
skipped=0
job_ids=()

for chrom in ${CHROMS}; do
    # Idempotency: skip if bgen already exists
    if dx ls "${DX_MERGED_BGEN_DIR}/chr${chrom}.bgen" &>/dev/null; then
        echo "    chr${chrom}: bgen already exists — skipping"
        skipped=$((skipped + 1))
        continue
    fi

    pfile="/mnt/project/${DX_MERGED_PFILE_DIR}/chr${chrom}"

    cmd="set -eo pipefail && \
echo '--- Converting chr${chrom} pfile to BGEN ---' && \
plink2 --pfile '${pfile}' --export bgen-1.2 bits=8 ref-first --out chr${chrom} && \
echo '--- Indexing chr${chrom}.bgen ---' && \
bgenix -g chr${chrom}.bgen -index -clobber && \
echo 'Done: chr${chrom}.bgen + chr${chrom}.bgen.bgi created'"

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_MERGED_BGEN_DIR}/" \
        --instance-type "${INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "convert_bgen_chr${chrom}" \
        -y --brief)

    echo "    chr${chrom}: submitted ${job_id}"
    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
done

echo "  Submitted: ${submitted}, Skipped: ${skipped}"

if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo "  Waiting for BGEN conversion jobs ..."
    for job_id in "${job_ids[@]}"; do
        dx wait "${job_id}"
    done
    echo "  All BGEN conversions complete."
fi
