#!/bin/bash
# fit_project_pca.sh — Fit PCA on unrelated Europeans, project onto all samples.
#
# Submits a Swiss Army Knife job that:
#   1. Fits PCA on pca_ready bfile (20 PCs, approx, allele-wts, seed 0)
#   2. Computes allele frequency counts for the fitting samples
#   3. Projects PCs onto all samples in direct_bfile/chr1_22_merged
#
# Output:       ${DX_PCA_EUR_DIR}/ukb_projected.sscore (+ PCA reference files)
# Intermediate: ${DX_PCA_EUR_DIR}/scrap/ (log files)
#
# Expects env vars: DX_PCA_EUR_DIR, DX_DIRECT_BFILE_DIR, INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

# Idempotency: skip if projected scores already exist
if dx ls "${DX_PCA_EUR_DIR}/ukb_projected.sscore" &>/dev/null; then
    echo "  ukb_projected.sscore already exists — skipping"
    exit 0
fi

dx mkdir -p "${DX_PCA_EUR_DIR}/scrap"

pca_bfile="/mnt/project${DX_PCA_EUR_DIR}/pca_ready"
direct_bfile="/mnt/project${DX_DIRECT_BFILE_DIR}/chr1_22_merged"

cmd="set -eo pipefail && \
echo '=== PCA Fitting & Projection ===' && \
\
echo '--- Step 1/3: Fit PCA on unrelated Europeans (20 PCs, approx) ---' && \
plink2 --bfile '${pca_bfile}' --pca allele-wts 20 approx --seed 0 --out ukb_pcs && \
echo '  Eigenvalues:' && \
cat ukb_pcs.eigenval && \
FIT_SNPS=\$(tail -n +2 ukb_pcs.eigenvec.allele | wc -l) && \
FIT_SAMPLES=\$(tail -n +2 ukb_pcs.eigenvec | wc -l) && \
echo \"  Fit on \${FIT_SNPS} SNPs, \${FIT_SAMPLES} samples\" && \
\
echo '--- Step 2/3: Compute allele frequency counts ---' && \
plink2 --bfile '${pca_bfile}' --freq counts --out pca_eur_counts && \
\
echo '--- Step 3/3: Project PCs onto all samples ---' && \
echo '  eigenvec.allele header:' && \
head -1 ukb_pcs.eigenvec.allele && \
A1_COL=\$(head -1 ukb_pcs.eigenvec.allele | tr '\t' '\n' | grep -n '^A1\$' | cut -d: -f1) && \
FIRST_PC=\$((A1_COL + 1)) && \
LAST_PC=\$((A1_COL + 20)) && \
echo \"  A1 column: \${A1_COL}, score columns: \${FIRST_PC}-\${LAST_PC}\" && \
awk 'NR>1 {print \$2}' ukb_pcs.eigenvec.allele > pca_snps.txt && \
echo \"  Extracting \$(wc -l < pca_snps.txt) PCA SNPs from direct bfile\" && \
plink2 --bfile '${direct_bfile}' \
    --extract pca_snps.txt \
    --read-freq pca_eur_counts.acount \
    --score ukb_pcs.eigenvec.allele 2 \${A1_COL} header-read no-mean-imputation variance-standardize \
    --score-col-nums \${FIRST_PC}-\${LAST_PC} \
    --out ukb_projected && \
PROJECTED=\$(tail -n +2 ukb_projected.sscore | wc -l) && \
echo \"  Projected \${PROJECTED} samples\" && \
\
echo '' && \
echo \"=== Done: \${FIT_SNPS} SNPs, \${FIT_SAMPLES} fit samples, \${PROJECTED} projected samples ===\" && \
\
echo '--- Cleanup intermediates ---' && \
rm -f pca_snps.txt && \
echo 'Done.'"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_PCA_EUR_DIR}/scrap/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "fit_project_pca" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Move key outputs from scrap/ to pca_eur/
dx mv "${DX_PCA_EUR_DIR}/scrap/ukb_projected.sscore" "${DX_PCA_EUR_DIR}/"
dx mv "${DX_PCA_EUR_DIR}/scrap/ukb_pcs.eigenvec" "${DX_PCA_EUR_DIR}/"
dx mv "${DX_PCA_EUR_DIR}/scrap/ukb_pcs.eigenvec.allele" "${DX_PCA_EUR_DIR}/"
dx mv "${DX_PCA_EUR_DIR}/scrap/ukb_pcs.eigenval" "${DX_PCA_EUR_DIR}/"
dx mv "${DX_PCA_EUR_DIR}/scrap/pca_eur_counts.acount" "${DX_PCA_EUR_DIR}/"

echo "  Done."
