#!/bin/bash
# imputed_extract_variants.sh — Submit per-chromosome imputed plink2 extraction jobs on DNAnexus.
#
# Submits 22 Swiss Army Knife jobs (one per chromosome) that run in parallel.
# Expects env vars: IMPUTED_BGEN_DIR, DX_TOPMED_ID_DIR, DX_IMPUTED_PFILE_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

dx mkdir -p "${DX_IMPUTED_PFILE_DIR}"

submitted=0
skipped=0
job_ids=()

for chrom in $(seq 1 22); do
    # Skip if output already exists on DNAnexus
    if dx ls "${DX_IMPUTED_PFILE_DIR}/chr${chrom}.pgen" &>/dev/null; then
        echo "chr${chrom}: skipping — pfiles already exist"
        skipped=$((skipped + 1))
        continue
    fi

    bgen_in="${IMPUTED_BGEN_DIR}/ukb21007_c${chrom}_b0_v1.bgen"
    sample_in="${IMPUTED_BGEN_DIR}/ukb21007_c${chrom}_b0_v1.sample"
    snp_list="/mnt/project${DX_TOPMED_ID_DIR}/chr${chrom}.txt"

    cmd="plink2 --bgen '${bgen_in}' ref-first \
  --sample '${sample_in}' \
  --set-all-var-ids '@:#:\$r:\$a' \
  --new-id-max-allele-len 7000 \
  --extract '${snp_list}' \
  --make-pgen \
  --out chr${chrom}"

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_IMPUTED_PFILE_DIR}/" \
        --instance-type "${INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "imputed_extract_chr${chrom}" \
        -y --brief)

    echo "chr${chrom}: submitted ${job_id}"
    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
done

echo ""
echo "Submitted: ${submitted}, Skipped (already exist): ${skipped}"

if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo ""
    echo "Job IDs:"
    printf '  %s\n' "${job_ids[@]}"
    echo ""
    echo "Monitor with: dx watch <job_id>"
fi
