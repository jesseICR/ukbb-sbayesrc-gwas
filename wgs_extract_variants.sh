#!/bin/bash
# wgs_extract_variants.sh — Submit per-chromosome WGS plink2 extraction jobs on DNAnexus.
#
# Submits 22 Swiss Army Knife jobs (one per chromosome) that run in parallel.
# Expects env vars: WGS_PLINK_DIR, DX_DRAGEN_ID_DIR, DX_WGS_PFILE_DIR, GENO, MAC, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

dx mkdir -p "${DX_WGS_PFILE_DIR}"

submitted=0
skipped=0
job_ids=()

for chrom in $(seq 1 22); do
    # Skip if output already exists on DNAnexus
    if dx ls "${DX_WGS_PFILE_DIR}/chr${chrom}.pgen" &>/dev/null; then
        echo "chr${chrom}: skipping — pfiles already exist"
        skipped=$((skipped + 1))
        continue
    fi

    plink_prefix="${WGS_PLINK_DIR}/ukb24308_c${chrom}_b0_v1"
    variant_file="/mnt/project${DX_DRAGEN_ID_DIR}/chr${chrom}.txt"

    cmd="plink2 \
  --pfile '${plink_prefix}' \
  --no-pheno \
  --extract '${variant_file}' \
  --geno ${GENO} \
  --mac ${MAC} \
  --var-filter PASS \
  --make-pgen \
  --out chr${chrom}"

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_WGS_PFILE_DIR}/" \
        --instance-type "${INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "wgs_extract_chr${chrom}" \
        -y --brief)

    echo "chr${chrom}: submitted ${job_id}"
    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
done

echo ""
echo "Submitted: ${submitted}, Skipped (already exist): ${skipped}"

if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo ""
    echo "Waiting for WGS extraction jobs ..."
    for job_id in "${job_ids[@]}"; do
        dx wait "${job_id}"
    done
    echo "All WGS extractions complete."
fi
