#!/bin/bash
# admixture_prep.sh — Prepare inputs for ADMIXTURE K=6 projection.
#
# 1. Downloads ADMIXTURE binary to tools/ if not already present
# 2. Uploads ADMIXTURE binary to DNAnexus scrap dir
# 3. Submits a SAK job that:
#    a. Downloads reference allele frequency TSV from GitHub
#    b. Extracts matching SNPs from the direct bfile
#    c. Aligns alleles (Python script) → .P file + aligned SNP list
#    d. Re-extracts bfile with aligned SNPs only
#
# Output on DNAnexus (in ${DX_ADMIXTURE_SCRAP_DIR}/):
#   admixture_allele_freqs.tsv, ref_aligned.P, admixture_align_log.txt,
#   ukb_admixture_aligned.{bed,bim,fam,log}, ukb_extracted.log, admixture
#
# Expects env vars: DX_DIRECT_BFILE_DIR, DX_ADMIXTURE_SCRAP_DIR,
#                   ADMIXTURE_TSV_URL, ADMIXTURE_DOWNLOAD_URL,
#                   ADMIXTURE_PREP_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="${SCRIPT_DIR}/tools"
ADMIXTURE_BIN="${TOOLS_DIR}/admixture"

# Idempotency: skip if aligned .P file already exists
if dx ls "${DX_ADMIXTURE_SCRAP_DIR}/ref_aligned.P" &>/dev/null; then
    echo "  ref_aligned.P already exists — skipping prep"
    exit 0
fi

# ---- Download ADMIXTURE binary locally if not present ----
if [[ ! -x "${ADMIXTURE_BIN}" ]]; then
    echo "  Downloading ADMIXTURE binary ..."
    mkdir -p "${TOOLS_DIR}"
    curl -fsSL -o "${TOOLS_DIR}/admixture_linux.tar.gz" "${ADMIXTURE_DOWNLOAD_URL}"
    tar xzf "${TOOLS_DIR}/admixture_linux.tar.gz" -C "${TOOLS_DIR}/"
    mv "${TOOLS_DIR}/dist/admixture_linux-1.3.0/admixture" "${ADMIXTURE_BIN}"
    chmod +x "${ADMIXTURE_BIN}"
    rm -rf "${TOOLS_DIR}/dist" "${TOOLS_DIR}/admixture_linux.tar.gz"
    echo "  Downloaded ADMIXTURE to ${ADMIXTURE_BIN}"
else
    echo "  ADMIXTURE binary already present at ${ADMIXTURE_BIN}"
fi

# ---- Upload ADMIXTURE binary to DNAnexus ----
dx mkdir -p "${DX_ADMIXTURE_SCRAP_DIR}"

if dx ls "${DX_ADMIXTURE_SCRAP_DIR}/admixture" &>/dev/null; then
    echo "  ADMIXTURE binary already on DNAnexus — skipping upload"
else
    echo "  Uploading ADMIXTURE binary to DNAnexus ..."
    dx upload "${ADMIXTURE_BIN}" --destination "${DX_ADMIXTURE_SCRAP_DIR}/" --brief --no-progress
    echo "  Upload complete."
fi

# ---- Upload Python alignment script as SAK input ----
script_id=$(dx upload "${SCRIPT_DIR}/admixture_align_alleles.py" \
    --destination "${DX_ADMIXTURE_SCRAP_DIR}/" --brief --no-progress)

# ---- Build and submit SAK job ----
bfile="/mnt/project${DX_DIRECT_BFILE_DIR}/chr1_22_merged"

cmd="set -eo pipefail && \
echo '--- Downloading reference allele frequencies ---' && \
curl -fsSL -o admixture_allele_freqs.tsv '${ADMIXTURE_TSV_URL}' && \
n_ref=\$(tail -n +2 admixture_allele_freqs.tsv | wc -l) && \
echo \"Reference TSV: \${n_ref} SNPs\" && \
echo '--- Extracting reference SNPs from direct bfile ---' && \
tail -n +2 admixture_allele_freqs.tsv | cut -f2 > ref_snp_ids.txt && \
plink2 --bfile '${bfile}' --extract ref_snp_ids.txt --make-bed --out ukb_extracted && \
echo \"Extracted: \$(wc -l < ukb_extracted.bim) SNPs, \$(wc -l < ukb_extracted.fam) individuals\" && \
echo '--- Aligning alleles ---' && \
python3 admixture_align_alleles.py && \
echo '--- Re-extracting with aligned SNPs ---' && \
plink2 --bfile ukb_extracted --extract snps_aligned.txt --make-bed --out ukb_admixture_aligned && \
echo \"Aligned bfile: \$(wc -l < ukb_admixture_aligned.bim) SNPs, \$(wc -l < ukb_admixture_aligned.fam) individuals\" && \
echo '--- Cleanup ---' && \
rm -f ukb_extracted.bed ukb_extracted.bim ukb_extracted.fam ukb_extracted.nosex \
      ref_snp_ids.txt snps_aligned.txt admixture_align_alleles.py"

job_id=$(dx run swiss-army-knife \
    -iin="${script_id}" \
    -icmd="${cmd}" \
    --destination "${DX_ADMIXTURE_SCRAP_DIR}/" \
    --instance-type "${ADMIXTURE_PREP_INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "admixture_prep_align" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"

# Clean up the uploaded Python script from DNAnexus
dx rm "${script_id}" 2>/dev/null || true

echo "  Done."
