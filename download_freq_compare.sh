#!/bin/bash
# download_freq_compare.sh — Download per-variant WGS vs TopMed-imputed frequency-
# comparison CSV from DNAnexus to the local data/freq_compare/ directory.
#
# This file is variant-level summary data (per-variant allele counts and sample
# totals), not individual-level — safe to download per project policy.

set -euo pipefail

mkdir -p "${LOCAL_FREQ_COMPARE_DIR}"

if [[ -s "${LOCAL_FREQ_COMPARE_CSV}" ]]; then
    echo "Already cached at ${LOCAL_FREQ_COMPARE_CSV} — skipping"
    exit 0
fi

echo "Downloading ${DX_FREQ_COMPARE_CSV} → ${LOCAL_FREQ_COMPARE_CSV}"
dx download "${DX_FREQ_COMPARE_CSV}" -o "${LOCAL_FREQ_COMPARE_CSV}"
echo "Downloaded ($(wc -l < "${LOCAL_FREQ_COMPARE_CSV}") lines)"
