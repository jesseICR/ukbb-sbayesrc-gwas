#!/bin/bash
# upload_dragen_ids.sh — Upload per-chromosome DRAGEN ID files to DNAnexus.
#
# Expects env vars: LOCAL_DRAGEN_ID_DIR, DX_DRAGEN_ID_DIR

set -euo pipefail

dx mkdir -p "${DX_DRAGEN_ID_DIR}"

uploaded=0
skipped=0

for chrom in $(seq 1 22); do
    filename="chr${chrom}.txt"
    local_path="${LOCAL_DRAGEN_ID_DIR}/${filename}"

    if ! [[ -f "${local_path}" ]]; then
        echo "ERROR: ${local_path} does not exist. Run step 1 first."
        exit 1
    fi

    if dx ls "${DX_DRAGEN_ID_DIR}/${filename}" &>/dev/null; then
        skipped=$((skipped + 1))
    else
        dx upload "${local_path}" --destination "${DX_DRAGEN_ID_DIR}/" --brief
        uploaded=$((uploaded + 1))
    fi
done

echo "Uploaded: ${uploaded}, Skipped (already exist): ${skipped}"
