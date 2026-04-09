#!/bin/bash
# upload_topmed_ids.sh — Upload per-chromosome TopMed ID files to DNAnexus.
#
# Expects env vars: LOCAL_TOPMED_ID_DIR, DX_TOPMED_ID_DIR

set -euo pipefail

dx mkdir -p "${DX_TOPMED_ID_DIR}"

uploaded=0
skipped=0

for chrom in $(seq 1 22); do
    filename="chr${chrom}.txt"
    local_path="${LOCAL_TOPMED_ID_DIR}/${filename}"

    if ! [[ -f "${local_path}" ]]; then
        echo "ERROR: ${local_path} does not exist. Run step 4 first."
        exit 1
    fi

    if dx ls "${DX_TOPMED_ID_DIR}/${filename}" &>/dev/null; then
        skipped=$((skipped + 1))
    else
        dx upload "${local_path}" --destination "${DX_TOPMED_ID_DIR}/" --brief
        uploaded=$((uploaded + 1))
    fi
done

echo "Uploaded: ${uploaded}, Skipped (already exist): ${skipped}"
