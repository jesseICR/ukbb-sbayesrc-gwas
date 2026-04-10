#!/bin/bash
# subset_kinship_snps.sh — Subset direct SNPs to kinship-relevant SNPs.
#
# Downloads ukb_snp_qc.txt (if not cached locally), filters to SNPs where
# in_Relatedness == 1, intersects with the existing QC-pass direct SNP list,
# and uploads the result to DNAnexus.
#
# Output:
#   Local:  data/support/ukbb_500k_qc_pass_direct_kinship_subsetted_snps.txt
#   Remote: ${DX_KINSHIP_SNPS_FILE}
#
# Expects env vars: LOCAL_SNP_QC_FILE, LOCAL_DIRECT_SNPS_FILE,
#                   LOCAL_KINSHIP_SNPS_FILE, DX_KINSHIP_DIR, DX_KINSHIP_SNPS_FILE

set -euo pipefail

# Idempotency: skip if already uploaded to DNAnexus
if dx ls "${DX_KINSHIP_SNPS_FILE}" &>/dev/null; then
    echo "  Kinship SNP list already uploaded — skipping"
    exit 0
fi

# Create temp directory with cleanup trap
tmpdir=$(mktemp -d)
trap 'rm -rf "${tmpdir}"' EXIT

# Download ukb_snp_qc.txt if not cached locally
if [[ -s "${LOCAL_SNP_QC_FILE}" ]]; then
    echo "  ukb_snp_qc.txt already cached locally — skipping download"
else
    echo "  Downloading ukb_snp_qc.txt ..."
    curl -fsSL -o "${LOCAL_SNP_QC_FILE}" \
        "https://biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_qc.txt"
    echo "  Downloaded ($(wc -l < "${LOCAL_SNP_QC_FILE}") lines)"
fi

# Find the in_Relatedness column index dynamically from the header
# (ukb_snp_qc.txt is space-delimited)
header=$(head -1 "${LOCAL_SNP_QC_FILE}")
rel_col=$(echo "${header}" | tr ' ' '\n' | grep -n '^in_Relatedness$' | cut -d: -f1)

if [[ -z "${rel_col}" ]]; then
    echo "ERROR: in_Relatedness column not found in ukb_snp_qc.txt header" >&2
    exit 1
fi
echo "  in_Relatedness is column ${rel_col}"

# Extract rsIDs where in_Relatedness == 1, then sort
awk -v col="${rel_col}" 'NR > 1 && $col == 1 {print $1}' \
    "${LOCAL_SNP_QC_FILE}" | sort > "${tmpdir}/relatedness_rsids.txt"
echo "  SNPs with in_Relatedness == 1: $(wc -l < "${tmpdir}/relatedness_rsids.txt")"

# Sort existing direct SNP list
sort "${LOCAL_DIRECT_SNPS_FILE}" > "${tmpdir}/direct_snps_sorted.txt"

# Intersect
comm -12 "${tmpdir}/relatedness_rsids.txt" "${tmpdir}/direct_snps_sorted.txt" \
    > "${LOCAL_KINSHIP_SNPS_FILE}"

n_kinship=$(wc -l < "${LOCAL_KINSHIP_SNPS_FILE}")
echo "  Kinship-relevant direct SNPs: ${n_kinship}"

# Upload to DNAnexus
dx mkdir -p "${DX_KINSHIP_DIR}"
dx upload "${LOCAL_KINSHIP_SNPS_FILE}" \
    --destination "${DX_KINSHIP_SNPS_FILE}" --brief
echo "  Uploaded to ${DX_KINSHIP_SNPS_FILE}"
