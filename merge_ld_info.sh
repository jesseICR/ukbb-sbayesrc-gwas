#!/bin/bash
# merge_ld_info.sh — Step 9: consolidate per-block LD info into snp.info (LDstep4).
#
# LDstep4 reads ldm.info plus every b{i}.ldm.full.info produced by Step 8 and
# emits a single snp.info file (Chrom ID Index GenPos PhysPos A1 A2 A1Freq N
# Block), which is the per-SNP lookup that downstream SBayesRC calls use.
#
# Upstream inputs:
#   - ${DX_LDM_DIR}/ldm.info                (Step 7)
#   - ${DX_LDM_DIR}/b*.ldm.full.info        (Step 8)
#   - ${DX_LDM_DIR}/sbayesrc_sak_setup.sh   (Step 7)
#
# Outputs:
#   - ${DX_LDM_DIR}/snp.info                (this step)
#
# Expects env vars: DX_LDM_DIR, DX_LDM_INFO_FILE, DX_LDM_SNP_INFO_FILE,
#                   LD_MATRICES_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

DX_SAK_SETUP="${DX_LDM_DIR}/sbayesrc_sak_setup.sh"

# Idempotency
if dx ls "${DX_LDM_SNP_INFO_FILE}" &>/dev/null; then
    echo "  snp.info already exists - skipping LDstep4"
    exit 0
fi

# Preconditions
if ! dx ls "${DX_LDM_INFO_FILE}" &>/dev/null; then
    echo "ERROR: ${DX_LDM_INFO_FILE} not found - run Step 7 first" >&2
    exit 1
fi

cmd=$(cat <<CMDEOF
set -eo pipefail
source /mnt/project${DX_SAK_SETUP}

echo '=== LDstep4: merge per-block info into snp.info ==='

mkdir -p ld_matrices
cp /mnt/project${DX_LDM_DIR}/ldm.info ld_matrices/

# Step 8 uploads per-chrom tarballs (chr{1..22}_meta.tar.gz) containing the
# b{i}.ldm.full.info files that LDstep4 needs. Extract all 22 into ld_matrices/.
echo "  extracting chr*_meta.tar.gz tarballs from Step 8 ..."
n_tars=0
for chrom in \$(seq 1 22); do
    tar_path="/mnt/project${DX_LDM_DIR}/chr\${chrom}_meta.tar.gz"
    if [[ -f "\$tar_path" ]]; then
        tar -xzf "\$tar_path" -C ld_matrices
        n_tars=\$((n_tars + 1))
    fi
done
n_info=\$(find ld_matrices -maxdepth 1 -type f -name 'b*.ldm.full.info' | wc -l)
echo "  extracted \${n_tars} tarballs -> \${n_info} b*.ldm.full.info files"

# Sanity: the number of info files should match the number of blocks in ldm.info.
n_blocks=\$(tail -n +2 ld_matrices/ldm.info | wc -l)
echo "  blocks in ldm.info:     \${n_blocks}"
echo "  b*.ldm.full.info files: \${n_info}"
if [[ "\${n_blocks}" != "\${n_info}" ]]; then
    echo "ERROR: block/info count mismatch - Step 8 is incomplete" >&2
    exit 1
fi

Rscript --no-save -e "SBayesRC::LDstep4(outDir='ld_matrices', log2file=FALSE)"

echo ''
echo '=== LDstep4 verification ==='
n_snps=\$(tail -n +2 ld_matrices/snp.info | wc -l)
n_chroms=\$(tail -n +2 ld_matrices/snp.info | awk '{print \$1}' | sort -u | wc -l)
echo "  SNPs in snp.info:       \${n_snps}"
echo "  unique chromosomes:     \${n_chroms}"
echo ''
echo 'snp.info header + first 3 rows:'
head -4 ld_matrices/snp.info

# Flatten output for SAK auto-upload. Keep only snp.info (the new output)
# and the LDstep4 log; everything else was already on DNAnexus.
mv ld_matrices/snp.info .
shopt -s nullglob
for f in ld_matrices/*.log; do
    mv "\$f" .
done
rm -rf ld_matrices
ls -lh
CMDEOF
)

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_LDM_DIR}/" \
    --instance-type "${LD_MATRICES_INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "ldstep4_merge" \
    --ignore-reuse \
    -y --brief)

echo "  Submitted ${job_id} - waiting for LDstep4 ..."
dx wait "${job_id}"
echo "  LDstep4 done."
