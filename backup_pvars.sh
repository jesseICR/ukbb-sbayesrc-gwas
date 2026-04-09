#!/bin/bash
# backup_pvars.sh — Back up original pvar files on DNAnexus before standardization.
#
# Submits DNAnexus jobs to copy WGS and imputed pvar files to backup directories,
# then waits for all jobs to complete so downstream steps can proceed.
#
# Expects env vars: DX_WGS_PFILE_DIR, DX_IMPUTED_PFILE_DIR, DX_BACKUP_DIR,
#                   INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

dx mkdir -p "${DX_BACKUP_DIR}/wgs_pvars"
dx mkdir -p "${DX_BACKUP_DIR}/imputed_pvars"

job_ids=()

for genotype_type in wgs imputed; do
    if [[ "${genotype_type}" == "wgs" ]]; then
        source_dir="${DX_WGS_PFILE_DIR}"
    else
        source_dir="${DX_IMPUTED_PFILE_DIR}"
    fi

    backed_up=0
    skipped=0

    for chrom in $(seq 1 22); do
        backup_path="${DX_BACKUP_DIR}/${genotype_type}_pvars/chr${chrom}.pvar"

        if dx ls "${backup_path}" &>/dev/null; then
            skipped=$((skipped + 1))
            continue
        fi

        cmd="cp /mnt/project/${source_dir}/chr${chrom}.pvar chr${chrom}.pvar"

        job_id=$(dx run swiss-army-knife \
            -icmd="${cmd}" \
            --destination "${DX_BACKUP_DIR}/${genotype_type}_pvars/" \
            --instance-type "${INSTANCE_TYPE}" \
            --priority "${DX_PRIORITY}" \
            --name "backup_${genotype_type}_pvar_chr${chrom}" \
            -y --brief)

        echo "chr${chrom} (${genotype_type}): submitted ${job_id}"
        job_ids+=("${job_id}")
        backed_up=$((backed_up + 1))
    done

    echo "${genotype_type}: Submitted: ${backed_up}, Skipped (already exist): ${skipped}"
done

if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo ""
    echo "Waiting for ${#job_ids[@]} backup jobs to complete..."
    for job_id in "${job_ids[@]}"; do
        dx wait "${job_id}"
    done
    echo "All backups complete."
fi
