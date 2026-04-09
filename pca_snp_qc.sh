#!/bin/bash
# pca_snp_qc.sh — QC SNPs for PCA: MAF filter, LD region exclusion, LD pruning.
#
# Submits a Swiss Army Knife job that:
#   1. Subsets direct bfile to fit_pca IIDs
#   2. Filters SNPs with MAF < 0.01
#   3. Excludes long-range LD regions (Price et al. 2008, hg38)
#   4. LD prunes (window=1000, step=80, r2=0.1)
#
# Output:       ${DX_PCA_EUR_DIR}/pca_ready.{bed,bim,fam,log}
# Intermediate: ${DX_PCA_EUR_DIR}/scrap/ (log files from each filter stage)
#
# Expects env vars: DX_PCA_EUR_DIR, DX_DIRECT_BFILE_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

# Idempotency: skip if final bfile already exists
if dx ls "${DX_PCA_EUR_DIR}/pca_ready.bed" &>/dev/null; then
    echo "  pca_ready.bed already exists — skipping"
    exit 0
fi

dx mkdir -p "${DX_PCA_EUR_DIR}/scrap"

direct_bfile="/mnt/project/${DX_DIRECT_BFILE_DIR}/chr1_22_merged"
fit_pca_iids="/mnt/project/${DX_PCA_EUR_DIR}/fit_pca_iids.txt"

cmd="set -eo pipefail && \
echo '=== PCA SNP QC Pipeline ===' && \
\
echo '--- Step 1/4: Subset direct bfile to fit_pca IIDs ---' && \
plink2 --bfile '${direct_bfile}' --keep '${fit_pca_iids}' --make-bed --out fit_pca_subset && \
SNPS=\$(wc -l < fit_pca_subset.bim) && \
SAMPLES=\$(wc -l < fit_pca_subset.fam) && \
echo \"  Subset: \${SNPS} SNPs, \${SAMPLES} samples\" && \
\
echo '--- Step 2/4: MAF filter (--maf 0.01) ---' && \
plink2 --bfile fit_pca_subset --maf 0.01 --make-bed --out maf_filtered && \
SNPS_AFTER=\$(wc -l < maf_filtered.bim) && \
echo \"  Removed \$((\${SNPS} - \${SNPS_AFTER})) SNPs with MAF < 0.01\" && \
echo \"  Remaining: \${SNPS_AFTER} SNPs\" && \
SNPS=\${SNPS_AFTER} && \
\
echo '--- Step 3/4: Exclude long-range LD regions (Price et al. 2008, hg38) ---' && \
curl -fsSL -o high_ld_regions_hg38.bed 'https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed' && \
awk 'BEGIN{OFS=\"\t\"} { chr=\$1; sub(/^chr/,\"\",chr); print chr, \$2+1, \$3, \$4 }' high_ld_regions_hg38.bed > high_ld_ranges.txt && \
plink2 --bfile maf_filtered --exclude range high_ld_ranges.txt --make-bed --out ldregion_excluded && \
SNPS_AFTER=\$(wc -l < ldregion_excluded.bim) && \
echo \"  Removed \$((\${SNPS} - \${SNPS_AFTER})) SNPs in long-range LD regions\" && \
echo \"  Remaining: \${SNPS_AFTER} SNPs\" && \
SNPS=\${SNPS_AFTER} && \
\
echo '--- Step 4/4: LD pruning (window=1000, step=80, r2=0.1) ---' && \
plink2 --bfile ldregion_excluded --indep-pairwise 1000 80 0.1 --out ld_prune && \
plink2 --bfile ldregion_excluded --extract ld_prune.prune.in --make-bed --out pca_ready && \
SNPS_AFTER=\$(wc -l < pca_ready.bim) && \
echo \"  Removed \$((\${SNPS} - \${SNPS_AFTER})) SNPs by LD pruning\" && \
echo \"  Remaining: \${SNPS_AFTER} SNPs\" && \
\
echo '' && \
echo \"=== Final PCA-ready bfile: \${SNPS_AFTER} SNPs, \$(wc -l < pca_ready.fam) samples ===\" && \
\
echo '--- Cleanup intermediate files ---' && \
rm -f fit_pca_subset.bed fit_pca_subset.bim fit_pca_subset.fam fit_pca_subset.nosex && \
rm -f maf_filtered.bed maf_filtered.bim maf_filtered.fam maf_filtered.nosex && \
rm -f ldregion_excluded.bed ldregion_excluded.bim ldregion_excluded.fam ldregion_excluded.nosex && \
rm -f ld_prune.prune.in ld_prune.prune.out ld_prune.nosex && \
rm -f high_ld_regions_hg38.bed high_ld_ranges.txt && \
echo 'Done.'"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_PCA_EUR_DIR}/scrap/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "pca_snp_qc" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Move final bfile from scrap/ to pca_eur/
dx mv "${DX_PCA_EUR_DIR}/scrap/pca_ready.bed" "${DX_PCA_EUR_DIR}/"
dx mv "${DX_PCA_EUR_DIR}/scrap/pca_ready.bim" "${DX_PCA_EUR_DIR}/"
dx mv "${DX_PCA_EUR_DIR}/scrap/pca_ready.fam" "${DX_PCA_EUR_DIR}/"
dx mv "${DX_PCA_EUR_DIR}/scrap/pca_ready.log" "${DX_PCA_EUR_DIR}/"

echo "  Done."
