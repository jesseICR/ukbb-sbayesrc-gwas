#!/bin/bash
# validate_merged_pfiles.sh — QC validation of merged pfiles.
#
# Checks that all 22 merged psam files have the expected sample count
# (WGS + imputed-only IIDs) and that variant counts in each merged pvar
# match the WGS variant counts from the per-chromosome merge logs
# (tolerance: ≤100 fewer).
#
# Uploads a QC log to sbayesrc_genotypes/merge_steps/pfile_merge_log.txt.
#
# Expects env vars: DX_MERGE_DIR, DX_MERGED_PFILE_DIR

set -euo pipefail

QC_LOG_DX="${DX_MERGE_DIR}/pfile_merge_log.txt"

# Idempotency: skip if QC log already exists
if dx ls "${QC_LOG_DX}" &>/dev/null; then
    echo "QC log already exists at ${QC_LOG_DX} — skipping"
    exit 0
fi

QC_LOG=$(mktemp)

{
    echo "=== Merged pfile QC log ==="
    echo "Date: $(date -Iseconds)"
    echo ""

    # --- Sample count check ---
    echo "--- Sample counts ---"
    wgs_n=$(dx cat "${DX_MERGE_DIR}/wgs_positive_iids.txt" | wc -l)
    imputed_only_n=$(dx cat "${DX_MERGE_DIR}/imputed_only_iids.txt" | wc -l)
    expected_n=$((wgs_n + imputed_only_n))
    echo "WGS IIDs: ${wgs_n}"
    echo "Imputed-only IIDs: ${imputed_only_n}"
    echo "Expected total: ${expected_n}"
    echo ""

    sample_ok=true
    for chrom in $(seq 1 22); do
        psam_n=$(dx cat "${DX_MERGED_PFILE_DIR}/chr${chrom}.psam" | tail -n +2 | wc -l)
        if [[ "${psam_n}" -ne "${expected_n}" ]]; then
            echo "chr${chrom}: ${psam_n} samples — MISMATCH (expected ${expected_n})"
            sample_ok=false
        else
            echo "chr${chrom}: ${psam_n} samples — OK"
        fi
    done
    echo ""

    # --- Variant count check ---
    echo "--- Variant counts (merged pvar vs WGS from merge log) ---"
    variant_ok=true
    for chrom in $(seq 1 22); do
        pvar_n=$(dx cat "${DX_MERGED_PFILE_DIR}/chr${chrom}.pvar" | tail -n +2 | wc -l)
        wgs_v=$(dx cat "${DX_MERGE_DIR}/bfiles/merge_chr${chrom}/chr${chrom}_merge_log.txt" | grep "^wgs_variants:" | awk '{print $2}')
        diff=$((wgs_v - pvar_n))
        if [[ "${diff}" -gt 100 ]]; then
            echo "chr${chrom}: pvar=${pvar_n} wgs=${wgs_v} fewer=${diff} — WARNING (>100 fewer)"
            variant_ok=false
        else
            echo "chr${chrom}: pvar=${pvar_n} wgs=${wgs_v} fewer=${diff} — OK"
        fi
    done
    echo ""

    # --- Summary ---
    echo "--- Summary ---"
    if ${sample_ok} && ${variant_ok}; then
        echo "QC PASSED: all sample counts match, all variant counts within tolerance"
    else
        ${sample_ok}  || echo "QC FAILED: sample count mismatch detected"
        ${variant_ok} || echo "QC FAILED: variant count deviation >100 detected"
    fi
} | tee "${QC_LOG}"

echo ""
echo "Uploading QC log to ${QC_LOG_DX} ..."
dx upload "${QC_LOG}" --destination "${QC_LOG_DX}" --brief
rm -f "${QC_LOG}"
