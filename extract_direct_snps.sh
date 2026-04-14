#!/bin/bash
# extract_direct_snps.sh — Extract direct SNPs from each per-chromosome pfile.
#
# Uploads the direct SNP list to DNAnexus (if not already present), then
# submits 22 parallel Swiss Army Knife jobs (one per chromosome) to extract
# matching SNPs from the merged pfiles into smaller per-chromosome pfiles.
#
# Output per chromosome in direct_pfiles/:
#   chr{N}.pgen, chr{N}.pvar, chr{N}.psam, chr{N}.log
#
# Expects env vars: DX_MERGED_PFILE_DIR, DX_DIRECT_PFILE_DIR,
#                   LOCAL_DIRECT_SNPS_FILE, DX_DIRECT_SNPS_FILE,
#                   INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

CHROMS=$(seq 1 22)

# Upload direct SNP list to DNAnexus if not already present
if dx ls "${DX_DIRECT_SNPS_FILE}" &>/dev/null; then
    echo "  Direct SNP list already uploaded — skipping upload"
else
    echo "  Uploading direct SNP list to DNAnexus ..."
    dx upload "${LOCAL_DIRECT_SNPS_FILE}" --destination "${DX_DIRECT_SNPS_FILE}" --brief
    echo "  Upload complete."
fi

echo "  Extracting direct SNPs per chromosome ..."

dx mkdir -p "${DX_DIRECT_PFILE_DIR}"

submitted=0
skipped=0
job_ids=()

for chrom in ${CHROMS}; do
    # Idempotency: skip if pfile already exists
    if dx ls "${DX_DIRECT_PFILE_DIR}/chr${chrom}.pgen" &>/dev/null; then
        echo "    chr${chrom}: direct pfile already exists — skipping"
        skipped=$((skipped + 1))
        continue
    fi

    pfile="/mnt/project${DX_MERGED_PFILE_DIR}/chr${chrom}"
    snps="/mnt/project${DX_DIRECT_SNPS_FILE}"

    cmd="set -eo pipefail && \
echo '--- Extracting direct SNPs for chr${chrom} ---' && \
plink2 --pfile '${pfile}' --extract '${snps}' --make-pgen --out chr${chrom} && \
echo 'Done: chr${chrom}.pgen/.pvar/.psam created'"

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_DIRECT_PFILE_DIR}/" \
        --instance-type "${INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "extract_direct_chr${chrom}" \
        --ignore-reuse \
        -y --brief)

    echo "    chr${chrom}: submitted ${job_id}"
    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
done

echo "  Submitted: ${submitted}, Skipped: ${skipped}"

if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo "  Waiting for extraction jobs ..."
    for job_id in "${job_ids[@]}"; do
        dx wait "${job_id}"
    done
    echo "  All extractions complete."
fi
