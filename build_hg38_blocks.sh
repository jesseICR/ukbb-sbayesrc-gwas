#!/bin/bash
# build_hg38_blocks.sh - Derive hg38 LD-block boundaries from
# sbayesrc_liftover_results.csv (no external liftOver needed).
#
# See build_hg38_blocks.py for the derivation logic and block-boundary
# convention. The liftover CSV is produced by Step 3 (qc_snps.sh) and
# already contains a per-SNP Block id and hg38 position; we just group
# by Block and take min/max.
#
# Purely local derivation + a single dx upload. Idempotent.
#
# Expects env vars: LOCAL_SBAYESRC_LIFTOVER_CSV, ALIGNMENT_FILE,
#                   LOCAL_LD_REF_DIR, LOCAL_BLOCKS_HG38,
#                   DX_LD_REF_DIR, DX_LD_BLOCKS_FILE

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -s "${LOCAL_BLOCKS_HG38}" ]]; then
    echo "  ${LOCAL_BLOCKS_HG38} already exists - skipping derivation"
else
    if [[ ! -s "${LOCAL_SBAYESRC_LIFTOVER_CSV}" ]]; then
        echo "ERROR: ${LOCAL_SBAYESRC_LIFTOVER_CSV} not found - run Step 3 first"
        exit 1
    fi
    if [[ ! -s "${ALIGNMENT_FILE}" ]]; then
        echo "ERROR: ${ALIGNMENT_FILE} not found - run get_genotypes.sh setup first"
        exit 1
    fi
    mkdir -p "${LOCAL_LD_REF_DIR}"
    python3 "${SCRIPT_DIR}/build_hg38_blocks.py" \
        --liftover "${LOCAL_SBAYESRC_LIFTOVER_CSV}" \
        --alignment "${ALIGNMENT_FILE}" \
        --output "${LOCAL_BLOCKS_HG38}"
fi

if dx ls "${DX_LD_BLOCKS_FILE}" &>/dev/null; then
    echo "  hg38 block file already on DNAnexus - skipping upload"
else
    dx mkdir -p "${DX_LD_REF_DIR}"
    echo "  Uploading ${LOCAL_BLOCKS_HG38} to DNAnexus ..."
    dx upload "${LOCAL_BLOCKS_HG38}" --destination "${DX_LD_BLOCKS_FILE}" \
        --brief --no-progress
    echo "  Upload complete."
fi
