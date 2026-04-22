#!/bin/bash
# build_ld_ref_bfiles.sh - Build per-chromosome LD-reference bfiles.
#
# Upstream inputs:
#   - ${LOCAL_QC_PASSED_CSV}              (Step 3 output, ~7.35M QC-passed rsids)
#   - ${DX_LD_COHORT_FILE}                (Step 4 output, 40k IIDs)
#   - ${DX_WGS_PFILE_DIR}/chr{1..22}      (WGS pfiles from get_genotypes.sh)
#
# Outputs (per chromosome in ${DX_LD_BFILE_DIR}/):
#   chr{N}.bed, chr{N}.bim, chr{N}.fam, chr{N}.log
#
# Submits 22 parallel Swiss Army Knife jobs. Per-chrom idempotent.
#
# Expects env vars: DX_WGS_PFILE_DIR, DX_LD_REF_DIR, DX_LD_COHORT_FILE,
#                   DX_LD_SNPS_FILE, DX_LD_BFILE_DIR, LOCAL_QC_PASSED_CSV,
#                   LD_REF_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

CHROMS=$(seq 1 22)

# ---- Upload QC-passed SNP list (variant_id column) to DNAnexus ----
if dx ls "${DX_LD_SNPS_FILE}" &>/dev/null; then
    echo "  LD-reference SNP list already on DNAnexus - skipping upload"
else
    if [[ ! -s "${LOCAL_QC_PASSED_CSV}" ]]; then
        echo "ERROR: ${LOCAL_QC_PASSED_CSV} not found - run Step 3 first"
        exit 1
    fi
    echo "  Extracting variant_id column from ${LOCAL_QC_PASSED_CSV} ..."
    tmp_snps=$(mktemp)
    # variant_id is the 2nd column in the QC-passed CSV
    tail -n +2 "${LOCAL_QC_PASSED_CSV}" | cut -d, -f2 > "${tmp_snps}"
    n_snps=$(wc -l < "${tmp_snps}")
    echo "  SNP list: ${n_snps} rsids"

    dx mkdir -p "${DX_LD_REF_DIR}"
    echo "  Uploading SNP list to DNAnexus ..."
    dx upload "${tmp_snps}" --destination "${DX_LD_SNPS_FILE}" \
        --brief --no-progress
    rm -f "${tmp_snps}"
    echo "  Upload complete."
fi

dx mkdir -p "${DX_LD_BFILE_DIR}"

echo "  Building per-chromosome LD-reference bfiles ..."

submitted=0
skipped=0
job_ids=()

for chrom in ${CHROMS}; do
    if dx ls "${DX_LD_BFILE_DIR}/chr${chrom}.bed" &>/dev/null; then
        echo "    chr${chrom}: bfile already exists - skipping"
        skipped=$((skipped + 1))
        continue
    fi

    pfile="/mnt/project${DX_WGS_PFILE_DIR}/chr${chrom}"
    keep="/mnt/project${DX_LD_COHORT_FILE}"
    snps="/mnt/project${DX_LD_SNPS_FILE}"

    cmd="set -eo pipefail && \
echo '--- chr${chrom}: filter WGS pfile to 40k cohort x QC-passed SNPs ---' && \
plink2 --pfile '${pfile}' \
  --keep '${keep}' \
  --extract '${snps}' \
  --make-bed \
  --out chr${chrom} && \
echo \"chr${chrom}: \$(wc -l < chr${chrom}.bim) SNPs, \$(wc -l < chr${chrom}.fam) samples\""

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_LD_BFILE_DIR}/" \
        --instance-type "${LD_REF_INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "ld_ref_bfile_chr${chrom}" \
        --ignore-reuse \
        -y --brief)

    echo "    chr${chrom}: submitted ${job_id}"
    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
done

echo "  Submitted: ${submitted}, Skipped: ${skipped}"

if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo "  Waiting for bfile jobs ..."
    for job_id in "${job_ids[@]}"; do
        dx wait "${job_id}"
    done
    echo "  All bfile jobs complete."
fi
