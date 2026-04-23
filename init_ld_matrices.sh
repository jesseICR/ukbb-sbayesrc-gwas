#!/bin/bash
# init_ld_matrices.sh — Step 7: initialize SBayesRC LD-matrix workspace (LDstep1).
#
# Synthesizes a dummy .ma file from the QC-passed variant list (LDstep1 requires
# an ma file even for LD-reference building; only SNP/A1/A2 are actually used
# for allele alignment — freq/b/se/p/N are present but immaterial for LD).
# Submits a single SAK job that runs SBayesRC::LDstep1, which reads the per-
# chromosome bim files via /mnt/project and writes ld.sh, ldm.info, and
# snplist/{block}.snplist into ${DX_LDM_DIR}.
#
# Upstream inputs:
#   - ${LOCAL_QC_PASSED_CSV}                    (Step 3 output, ~7.35M rows)
#   - ${DX_LD_BFILE_DIR}/chr{1..22}.bim         (Step 6 outputs)
#   - ${DX_LD_BLOCKS_FILE} (ref4cM_hg38.pos)    (Step 5 output)
#
# Outputs (in ${DX_LDM_DIR}):
#   - dummy_ma.ma           (this script synthesizes + uploads)
#   - sbayesrc_sak_setup.sh (this script uploads; used by Steps 7, 8, 9)
#   - ld.sh                 (LDstep1; gctb commands, one per block)
#   - ldm.info              (LDstep1; per-block metadata)
#   - snplist/{block}.snplist  (LDstep1; one file per retained block)
#
# Expects env vars: LOCAL_QC_PASSED_CSV, DX_LD_BFILE_DIR, DX_LD_BLOCKS_FILE,
#                   DX_LDM_DIR, DX_LDM_INFO_FILE, DX_LDM_MA_FILE,
#                   LD_COHORT_SIZE, LD_MATRICES_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Idempotency: skip if ldm.info is already on DNAnexus ----
if dx ls "${DX_LDM_INFO_FILE}" &>/dev/null; then
    echo "  ldm.info already exists - skipping LDstep1"
    exit 0
fi

dx mkdir -p "${DX_LDM_DIR}"

# ---- Ensure the SAK setup helper is uploaded (shared by Steps 7/8/9) ----
# Always re-upload so local changes to the helper propagate on re-runs.
DX_SAK_SETUP="${DX_LDM_DIR}/sbayesrc_sak_setup.sh"
if dx ls "${DX_SAK_SETUP}" &>/dev/null; then
    dx rm "${DX_SAK_SETUP}" 2>/dev/null || true
fi
echo "  Uploading sbayesrc_sak_setup.sh ..."
dx upload "${SCRIPT_DIR}/sbayesrc_sak_setup.sh" \
    --destination "${DX_SAK_SETUP}" --brief --no-progress

# ---- Synthesize + upload dummy_ma.ma ----
if ! dx ls "${DX_LDM_MA_FILE}" &>/dev/null; then
    if [[ ! -s "${LOCAL_QC_PASSED_CSV}" ]]; then
        echo "ERROR: ${LOCAL_QC_PASSED_CSV} not found - run Step 3 first" >&2
        exit 1
    fi
    echo "  Synthesizing dummy_ma.ma from QC-passed CSV ..."
    tmp_ma=$(mktemp --suffix=.ma)
    python3 - "${LOCAL_QC_PASSED_CSV}" "${tmp_ma}" "${LD_COHORT_SIZE}" <<'PYEOF'
import csv
import sys

in_csv, out_ma, n = sys.argv[1], sys.argv[2], int(sys.argv[3])

with open(in_csv) as fin, open(out_ma, "w") as fout:
    reader = csv.DictReader(fin)
    fout.write("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n")
    count = 0
    for row in reader:
        # ma columns: SNP, A1 (effect = alt), A2 (ref), freq (alt freq in WGS),
        # b/se/p dummy (not used by LDstep1), N = LD cohort size.
        fout.write(
            f"{row['variant_id']}\t{row['alt']}\t{row['ref']}\t"
            f"{row['alt_freq_wgs']}\t0\t1\t1\t{n}\n"
        )
        count += 1
print(f"  wrote {count} SNPs to {out_ma}", file=sys.stderr)
PYEOF
    n_ma=$(($(wc -l < "${tmp_ma}") - 1))
    echo "  dummy_ma.ma: ${n_ma} SNPs"

    echo "  Uploading dummy_ma.ma to DNAnexus ..."
    dx upload "${tmp_ma}" --destination "${DX_LDM_MA_FILE}" \
        --brief --no-progress
    rm -f "${tmp_ma}"
fi

# ---- Submit LDstep1 SAK job ----
# The outer heredoc is UNQUOTED so ${DX_*} env vars expand at wrapper time;
# runtime vars are escaped with \$ to pass through to SAK's shell.
cmd=$(cat <<CMDEOF
set -eo pipefail
source /mnt/project${DX_SAK_SETUP}

echo '=== LDstep1: build LD-matrix workspace ==='
echo "  ma file:     /mnt/project${DX_LDM_MA_FILE}"
echo "  bfile prefix: /mnt/project${DX_LD_BFILE_DIR}/chr{CHR}"
echo "  block ref:   /mnt/project${DX_LD_BLOCKS_FILE}"

Rscript --no-save -e "SBayesRC::LDstep1(
    mafile='/mnt/project${DX_LDM_MA_FILE}',
    genoPrefix='/mnt/project${DX_LD_BFILE_DIR}/chr{CHR}',
    outDir='ld_matrices',
    genoCHR='1-22',
    blockRef='/mnt/project${DX_LD_BLOCKS_FILE}',
    log2file=FALSE)"

echo ''
echo '=== LDstep1 verification ==='
n_blocks=\$(tail -n +2 ld_matrices/ldm.info | wc -l)
n_snps=\$(tail -n +2 ld_matrices/ldm.info | awk '{s += \$NF} END {print s}')
n_snplists=\$(ls ld_matrices/snplist | wc -l)
echo "  blocks in ldm.info:     \${n_blocks}"
echo "  total SNPs (sum NumSnps): \${n_snps}"
echo "  snplist files:          \${n_snplists}"
echo ''
echo 'ldm.info header + first 3 rows:'
head -4 ld_matrices/ldm.info
echo ''
echo 'ld.sh first 2 lines:'
head -2 ld_matrices/ld.sh

# Flatten for SAK auto-upload.
# SAK uploads each cwd file via a separate API call, so tar the 591 small
# snplist files into a single snplist.tar.gz (one API call instead of 591).
# Step 8 (build_ld_matrices.sh) extracts it when staging the workspace.
mv ld_matrices/ld.sh ld_matrices/ldm.info .
tar -czf snplist.tar.gz -C ld_matrices snplist
rm -rf ld_matrices
ls -lh
CMDEOF
)

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_LDM_DIR}/" \
    --instance-type "${LD_MATRICES_INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "ldstep1_init" \
    --ignore-reuse \
    -y --brief)

echo "  Submitted ${job_id} - waiting for LDstep1 ..."
dx wait "${job_id}"
echo "  LDstep1 done."
