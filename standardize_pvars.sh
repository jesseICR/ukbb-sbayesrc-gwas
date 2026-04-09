#!/bin/bash
# standardize_pvars.sh — Standardize pvar files to #CHROM/POS/ID/REF/ALT with rsid lookups.
#
# Downloads each pvar from DNAnexus, remaps variant IDs to rsids using
# alignment CSVs, uploads the standardized pvar back, and cleans up locally.
#
# Expects env vars: DX_WGS_PFILE_DIR, DX_IMPUTED_PFILE_DIR, DX_BACKUP_DIR, ALIGNMENT_FILE

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

WORK_DIR=$(mktemp -d)
trap 'rm -rf "${WORK_DIR}"' EXIT

for genotype_type in wgs imputed; do
    if [[ "${genotype_type}" == "wgs" ]]; then
        dx_pfile_dir="${DX_WGS_PFILE_DIR}"
    else
        dx_pfile_dir="${DX_IMPUTED_PFILE_DIR}"
    fi

    sentinel="${dx_pfile_dir}/.pvars_standardized"

    # Skip entirely if sentinel exists — all pvars already standardized
    if dx ls "${sentinel}" &>/dev/null; then
        echo "${genotype_type}: all pvars already standardized, skipping"
        continue
    fi

    for chrom in $(seq 1 22); do
        backup_path="${DX_BACKUP_DIR}/${genotype_type}_pvars/chr${chrom}.pvar"
        remote_pvar="${dx_pfile_dir}/chr${chrom}.pvar"
        local_pvar="${WORK_DIR}/chr${chrom}.pvar"
        local_output="${WORK_DIR}/chr${chrom}_standardized.pvar"

        # Verify backup exists — step 7 must have run first
        if ! dx ls "${backup_path}" &>/dev/null; then
            echo "ERROR: Backup not found at ${backup_path}. Run step 7 first."
            exit 1
        fi

        # Download and transform
        echo "chr${chrom} (${genotype_type}): standardizing..."
        dx download "${remote_pvar}" -o "${local_pvar}" -f
        python3 "${SCRIPT_DIR}/standardize_pvar.py" \
            --pvar "${local_pvar}" \
            --alignment "${ALIGNMENT_FILE}" \
            --chrom "${chrom}" \
            --output "${local_output}"

        # Rename locally to match the original filename, then replace on DNAnexus
        mv "${local_output}" "${local_pvar}"
        dx rm "${remote_pvar}"
        dx upload "${local_pvar}" --destination "${dx_pfile_dir}/" --brief

        # Clean up local files
        rm -f "${local_pvar}"

    done

    echo "${genotype_type}: all 22 chromosomes standardized"

    # Write sentinel so future runs skip instantly
    echo "standardized on $(date -Iseconds)" | dx upload - --destination "${sentinel}" --brief
done
