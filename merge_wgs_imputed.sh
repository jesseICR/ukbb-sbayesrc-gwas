#!/bin/bash
# merge_wgs_imputed.sh — Merge WGS + imputed-only individuals per chromosome.
#
# For each chromosome, submits a Swiss Army Knife job that:
#   1. Finds common variants (matching ID + REF + ALT) between WGS and imputed pvars
#   2. Extracts WGS bfile (common variants, positive IIDs only)
#   3. Extracts imputed bfile (common variants, positive imputed-only IIDs only)
#   4. Merges the two bfiles with plink1 (trial merge + missnp handling)
#   5. Writes a merge log with key counts
#
# Then converts each merged bfile to a pfile in merged_pfiles/.
#
# Outputs per chromosome:
#   merge_steps/bfiles/merge_chr{N}/  — common_variants.txt, chr{N}.{bed,bim,fam}, log
#   merged_pfiles/                    — chr{N}.{pgen,pvar,psam}
#
# Expects env vars: DX_WGS_PFILE_DIR, DX_IMPUTED_PFILE_DIR, DX_MERGE_DIR,
#                   DX_MERGED_PFILE_DIR, MERGE_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

# ---------------------------------------------------------------------------
# Dependency check
# ---------------------------------------------------------------------------
if ! dx ls "${DX_MERGE_DIR}/imputed_only_iids.txt" &>/dev/null; then
    echo "ERROR: ${DX_MERGE_DIR}/imputed_only_iids.txt not found. Run step 9 first."
    exit 1
fi
if ! dx ls "${DX_MERGE_DIR}/wgs_positive_iids.txt" &>/dev/null; then
    echo "ERROR: ${DX_MERGE_DIR}/wgs_positive_iids.txt not found. Run step 9 first."
    exit 1
fi

# ---------------------------------------------------------------------------
# Phase 1: Merge WGS + imputed-only per chromosome
# ---------------------------------------------------------------------------
echo "  Phase 1: Merging WGS + imputed-only ..."

p1_submitted=0
p1_skipped=0
p1_job_ids=()

for chrom in $(seq 1 22); do
    dest="${DX_MERGE_DIR}/bfiles/merge_chr${chrom}"

    # Idempotency: skip if merged bfile already exists
    if dx ls "${dest}/chr${chrom}.bed" &>/dev/null; then
        echo "    chr${chrom}: merged bfile already exists — skipping"
        p1_skipped=$((p1_skipped + 1))
        continue
    fi

    dx mkdir -p "${dest}"

    wgs_pvar="/mnt/project${DX_WGS_PFILE_DIR}/chr${chrom}.pvar"
    imp_pvar="/mnt/project${DX_IMPUTED_PFILE_DIR}/chr${chrom}.pvar"
    wgs_pfile="/mnt/project${DX_WGS_PFILE_DIR}/chr${chrom}"
    imp_pfile="/mnt/project${DX_IMPUTED_PFILE_DIR}/chr${chrom}"
    keep_wgs="/mnt/project${DX_MERGE_DIR}/wgs_positive_iids.txt"
    keep_imp="/mnt/project${DX_MERGE_DIR}/imputed_only_iids.txt"

    cmd="set -eo pipefail && \
\
echo '--- Step 1: Finding common variants ---' && \
wgs_v=\$(tail -n +2 '${wgs_pvar}' | wc -l) && \
imp_v=\$(tail -n +2 '${imp_pvar}' | wc -l) && \
echo \"WGS variants: \${wgs_v}\" && \
echo \"Imputed variants: \${imp_v}\" && \
tail -n +2 '${wgs_pvar}' | awk -F'\t' '{print \$3\"\t\"\$4\"\t\"\$5}' | sort -u > wgs_vars.txt && \
tail -n +2 '${imp_pvar}' | awk -F'\t' '{print \$3\"\t\"\$4\"\t\"\$5}' | sort -u > imp_vars.txt && \
comm -12 wgs_vars.txt imp_vars.txt | awk -F'\t' '{print \$1}' | sort -u > common_variants.txt && \
common_n=\$(wc -l < common_variants.txt) && \
echo \"Common variants: \${common_n}\" && \
rm -f wgs_vars.txt imp_vars.txt && \
\
echo '--- Step 2: Extracting WGS bfile ---' && \
plink2 --pfile '${wgs_pfile}' --extract common_variants.txt --keep '${keep_wgs}' --make-bed --out wgs && \
wgs_n=\$(wc -l < wgs.fam) && \
echo \"WGS samples: \${wgs_n}\" && \
\
echo '--- Step 3: Extracting imputed-only bfile ---' && \
plink2 --pfile '${imp_pfile}' --extract common_variants.txt --keep '${keep_imp}' --make-bed --out imp && \
imp_n=\$(wc -l < imp.fam) && \
echo \"Imputed-only samples: \${imp_n}\" && \
\
echo '--- Step 4: Merging bfiles ---' && \
missnp_n=0 && \
plink --bfile wgs --bmerge imp.bed imp.bim imp.fam --make-bed --out merged || true && \
\
if [ -f merged-merge.missnp ]; then \
    missnp_n=\$(wc -l < merged-merge.missnp) && \
    echo \"Mismatched SNPs: \${missnp_n} — re-merging without them\" && \
    plink --bfile wgs --exclude merged-merge.missnp --make-bed --out wgs_clean && \
    plink --bfile imp --exclude merged-merge.missnp --make-bed --out imp_clean && \
    plink --bfile wgs_clean --bmerge imp_clean.bed imp_clean.bim imp_clean.fam --make-bed --out merged && \
    rm -f wgs_clean.bed wgs_clean.bim wgs_clean.fam wgs_clean.nosex && \
    rm -f imp_clean.bed imp_clean.bim imp_clean.fam imp_clean.nosex; \
fi && \
\
merged_n=\$(wc -l < merged.fam) && \
merged_v=\$(wc -l < merged.bim) && \
echo \"Merged samples: \${merged_n}\" && \
echo \"Merged variants: \${merged_v}\" && \
\
echo '--- Writing log ---' && \
echo \"chr${chrom} merge summary\" > chr${chrom}_merge_log.txt && \
echo \"wgs_variants: \${wgs_v}\" >> chr${chrom}_merge_log.txt && \
echo \"imputed_variants: \${imp_v}\" >> chr${chrom}_merge_log.txt && \
echo \"common_variants: \${common_n}\" >> chr${chrom}_merge_log.txt && \
echo \"wgs_samples: \${wgs_n}\" >> chr${chrom}_merge_log.txt && \
echo \"imputed_only_samples: \${imp_n}\" >> chr${chrom}_merge_log.txt && \
echo \"mismatched_snps_excluded: \${missnp_n}\" >> chr${chrom}_merge_log.txt && \
echo \"merged_samples: \${merged_n}\" >> chr${chrom}_merge_log.txt && \
echo \"merged_variants: \${merged_v}\" >> chr${chrom}_merge_log.txt && \
\
echo '--- Cleanup ---' && \
mv merged.bed chr${chrom}.bed && \
mv merged.bim chr${chrom}.bim && \
mv merged.fam chr${chrom}.fam && \
rm -f wgs.bed wgs.bim wgs.fam wgs.nosex && \
rm -f imp.bed imp.bim imp.fam imp.nosex && \
rm -f merged.bed merged.bim merged.fam merged.nosex merged-merge.missnp"

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${dest}/" \
        --instance-type "${MERGE_INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "merge_wgs_imputed_chr${chrom}" \
        -y --brief)

    echo "    chr${chrom}: submitted ${job_id}"
    p1_job_ids+=("${job_id}")
    p1_submitted=$((p1_submitted + 1))
done

echo "  Phase 1 — Submitted: ${p1_submitted}, Skipped: ${p1_skipped}"

if [[ ${#p1_job_ids[@]} -gt 0 ]]; then
    echo "  Waiting for phase 1 jobs ..."
    for job_id in "${p1_job_ids[@]}"; do
        dx wait "${job_id}"
    done
    echo "  Phase 1 complete."
fi

# ---------------------------------------------------------------------------
# Phase 2: Convert merged bfiles to pfiles
# ---------------------------------------------------------------------------
echo ""
echo "  Phase 2: Converting merged bfiles to pfiles ..."

dx mkdir -p "${DX_MERGED_PFILE_DIR}"

p2_submitted=0
p2_skipped=0
p2_job_ids=()

for chrom in $(seq 1 22); do
    # Idempotency: skip if pfile already exists
    if dx ls "${DX_MERGED_PFILE_DIR}/chr${chrom}.pgen" &>/dev/null; then
        echo "    chr${chrom}: merged pfile already exists — skipping"
        p2_skipped=$((p2_skipped + 1))
        continue
    fi

    bfile_dir="/mnt/project${DX_MERGE_DIR}/bfiles/merge_chr${chrom}"

    cmd="set -eo pipefail && \
plink2 --bfile '${bfile_dir}/chr${chrom}' --make-pgen --out chr${chrom} && \
echo \"Converted chr${chrom}: \$(wc -l < chr${chrom}.psam) samples, \$(wc -l < chr${chrom}.pvar) variants\""

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_MERGED_PFILE_DIR}/" \
        --instance-type "${MERGE_INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "convert_merged_pfile_chr${chrom}" \
        -y --brief)

    echo "    chr${chrom}: submitted ${job_id}"
    p2_job_ids+=("${job_id}")
    p2_submitted=$((p2_submitted + 1))
done

echo "  Phase 2 — Submitted: ${p2_submitted}, Skipped: ${p2_skipped}"

if [[ ${#p2_job_ids[@]} -gt 0 ]]; then
    echo "  Waiting for phase 2 jobs ..."
    for job_id in "${p2_job_ids[@]}"; do
        dx wait "${job_id}"
    done
    echo "  Phase 2 complete."
fi
