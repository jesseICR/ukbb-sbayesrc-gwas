#!/bin/bash
# build_ld_matrices.sh — Step 8: per-chromosome LD-matrix construction + eigen.
#
# For each chromosome C (1..22) submits a single SAK job that:
#   - stages the LDstep1 outputs into a writable local dir;
#   - for each block assigned to chromosome C (per ldm.info):
#       * LDstep2  -> gctb --make-full-ldm  (builds b{i}.ldm.full.{bin,info})
#       * LDstep3  -> eigen decomposition    (builds block{i}.eigen.bin)
#       * delete b{i}.ldm.full.bin (~600 MB/block) once eigen is done; retain
#         b{i}.ldm.full.info (needed by LDstep4).
#   - auto-uploads block{i}.eigen.bin, block{i}.eigen.bin.log, and
#     b{i}.ldm.full.info into ${DX_LDM_DIR}.
#
# Per-chromosome idempotency: if every block{i}.eigen.bin for that chromosome
# already exists on DNAnexus, the chromosome is skipped. Per-block idempotency
# is handled inside the SAK job via a `dx ls` check against /mnt/project.
#
# Upstream inputs:
#   - ${DX_LDM_DIR}/{ld.sh, ldm.info, snplist/}   (Step 7 outputs)
#   - ${DX_LDM_DIR}/sbayesrc_sak_setup.sh         (Step 7 upload)
#   - ${DX_LD_BFILE_DIR}/chr{N}.{bed,bim,fam}     (Step 6 outputs)
#
# Expects env vars: DX_LDM_DIR, DX_LDM_INFO_FILE, LD_MATRICES_INSTANCE_TYPE,
#                   LD_MATRICES_OMP_THREADS, DX_PRIORITY

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

DX_SAK_SETUP="${DX_LDM_DIR}/sbayesrc_sak_setup.sh"

# ---- Precondition: ldm.info must exist (Step 7 must have run) ----
if ! dx ls "${DX_LDM_INFO_FILE}" &>/dev/null; then
    echo "ERROR: ${DX_LDM_INFO_FILE} not found - run Step 7 first" >&2
    exit 1
fi

# ---- Fetch ldm.info locally to plan per-chrom block lists ----
LOCAL_LDM_INFO="${SCRIPT_DIR}/data/ld_reference/ldm.info"
mkdir -p "$(dirname "${LOCAL_LDM_INFO}")"
if [[ ! -s "${LOCAL_LDM_INFO}" ]]; then
    echo "  Downloading ldm.info locally for block-planning ..."
    dx download "${DX_LDM_INFO_FILE}" -o "${LOCAL_LDM_INFO}" --no-progress
fi

# Columns of ldm.info: Block  Chrom  StartSnpIdx  StartSnpID  EndSnpIdx  EndSnpID  NumSnps
# We need, per chromosome, the list of (row_idx, Block) pairs where row_idx is
# 1-indexed into the data rows (as used by LDstep2/3's `blockIndex` argument).

# Build a tmp file with 3 cols: chrom row_idx block
tmp_plan=$(mktemp)
tail -n +2 "${LOCAL_LDM_INFO}" | awk '{print $2, NR, $1}' > "${tmp_plan}"

echo "  Per-chromosome block plan:"
for chrom in $(seq 1 22); do
    n=$(awk -v c="$chrom" '$1==c' "${tmp_plan}" | wc -l)
    echo "    chr${chrom}: ${n} blocks"
done

echo ""
echo "  Submitting per-chromosome SAK jobs ..."

submitted=0
skipped=0
job_ids=()

for chrom in $(seq 1 22); do
    # blocks on this chrom: "row_idx:block row_idx:block ..."
    pairs=$(awk -v c="$chrom" '$1==c {printf "%d:%d ", $2, $3}' "${tmp_plan}")
    if [[ -z "${pairs// }" ]]; then
        echo "    chr${chrom}: no blocks - skipping"
        continue
    fi

    # Per-chrom idempotency: all block{N}.eigen.bin exist?
    all_done=1
    blocks_needed=()
    for p in ${pairs}; do
        block="${p#*:}"
        if ! dx ls "${DX_LDM_DIR}/block${block}.eigen.bin" &>/dev/null; then
            all_done=0
            blocks_needed+=("${block}")
        fi
    done
    if [[ "${all_done}" == "1" ]]; then
        echo "    chr${chrom}: all eigen.bin files present - skipping"
        skipped=$((skipped + 1))
        continue
    fi

    echo "    chr${chrom}: ${#blocks_needed[@]} blocks to compute"

    # The idx_block_pairs string is evaluated at wrapper time so the SAK cmd
    # hard-codes the (row_idx,block) list for this chromosome.
    idx_block_pairs="${pairs% }"

    cmd=$(cat <<CMDEOF
set -eo pipefail
source /mnt/project${DX_SAK_SETUP}

export OMP_NUM_THREADS=${LD_MATRICES_OMP_THREADS}

echo ''
echo '=== LDstep2 + LDstep3 for chr${chrom} ==='

# Stage LDstep1 outputs (they need to live under a writable 'ld_matrices/'
# so that LDstep2 can write b{block}.ldm.full.{bin,info} inside it).
# snplist/ was uploaded as a single snplist.tar.gz by Step 7 to avoid SAK's
# per-file upload API tax.
mkdir -p ld_matrices
cp /mnt/project${DX_LDM_DIR}/ld.sh     ld_matrices/
cp /mnt/project${DX_LDM_DIR}/ldm.info  ld_matrices/
tar -xzf /mnt/project${DX_LDM_DIR}/snplist.tar.gz -C ld_matrices

# SBayesRC writes ld.sh with bash-specific '&>' redirection, but R's system()
# dispatches to /bin/sh (dash on Ubuntu) which parses 'cmd &> file' as
# 'cmd &' (background) + '> file' (empty redirect). The consequence is that
# gctb is backgrounded and LDstep2 returns before gctb finishes writing
# b{i}.ldm.full.{bin,info}, causing LDstep3 to fail reading them. Rewrite
# to POSIX-compliant '> file 2>&1' instead.
#
# Also insert '--thread LD_MATRICES_OMP_THREADS' before --out so gctb uses
# all cores on the worker for --make-full-ldm (single-threaded default
# would make each ~10k-SNP block take 30+ min instead of 3-5 min).
sed -i "s| --out | --thread ${LD_MATRICES_OMP_THREADS} --out |; s/ &> \([^[:space:]]*\)\$/ > \1 2>\&1/" ld_matrices/ld.sh

n_processed=0
n_skipped=0
for pair in ${idx_block_pairs}; do
    row_idx=\${pair%:*}
    block=\${pair#*:}

    # Per-block idempotency: skip if the eigen.bin already on DNAnexus.
    if [[ -f /mnt/project${DX_LDM_DIR}/block\${block}.eigen.bin ]]; then
        echo "  chr${chrom} block \${block} (row \${row_idx}): eigen.bin already exists - skipping"
        n_skipped=\$((n_skipped + 1))
        continue
    fi

    echo ''
    echo "--- chr${chrom} block \${block} (row \${row_idx}) ---"

    # LDstep2: build full LD matrix via gctb --make-full-ldm.
    Rscript --no-save -e "SBayesRC::LDstep2(outDir='ld_matrices', blockIndex=\${row_idx}, log2file=FALSE)"

    # LDstep3: eigen decomposition (respects OMP_NUM_THREADS).
    Rscript --no-save -e "SBayesRC::LDstep3(outDir='ld_matrices', blockIndex=\${row_idx}, log2file=FALSE)"

    # Sanity: eigen header (m, k, sumLambda, thresh).
    Rscript --no-save -e "e <- SBayesRC::readEig('ld_matrices', \${block}); cat('  block \${block}: m=', e\\\$m, ' k=', e\\\$k, ' sumLambda=', e\\\$sumLambda, '\n', sep='')"

    # Free the huge full-LD bin; keep .info for LDstep4.
    rm -f ld_matrices/b\${block}.ldm.full.bin

    n_processed=\$((n_processed + 1))
done

echo ''
echo "=== chr${chrom} summary: processed=\${n_processed}, skipped_existing=\${n_skipped} ==="

# Flatten outputs for SAK auto-upload (destination is DX_LDM_DIR).
# We only want to upload NEW files - suppress re-uploading staged inputs
# (ld.sh, ldm.info, snplist/) which came from /mnt/project and are unchanged.
rm -f ld_matrices/ld.sh ld_matrices/ldm.info
rm -rf ld_matrices/snplist

# Each block{i}.eigen.bin is ~hundreds of MB - upload individually. But the
# b{i}.ldm.full.info + *.log sidecars are tiny and number in the dozens per
# chrom; SAK pays one API call per file, so tar them into a single bundle.
shopt -s nullglob
for f in ld_matrices/block*.eigen.bin; do
    mv "\$f" .
done
tar -czf chr${chrom}_meta.tar.gz -C ld_matrices . 2>/dev/null || true
rm -rf ld_matrices
ls -lh | head -40
echo "  total output files: \$(ls -1 | wc -l)"
CMDEOF
)

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_LDM_DIR}/" \
        --instance-type "${LD_MATRICES_INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "ldstep23_chr${chrom}" \
        --ignore-reuse \
        -y --brief)

    echo "    chr${chrom}: submitted ${job_id}"
    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
done

rm -f "${tmp_plan}"

echo ""
echo "  Submitted: ${submitted}, Skipped: ${skipped}"

if [[ ${#job_ids[@]} -gt 0 ]]; then
    echo "  Waiting for LDstep2+3 jobs (this can take hours) ..."
    failed=0
    for job_id in "${job_ids[@]}"; do
        if ! dx wait "${job_id}"; then
            echo "  WARNING: ${job_id} did not exit cleanly" >&2
            failed=$((failed + 1))
        fi
    done
    if (( failed > 0 )); then
        echo "ERROR: ${failed} LDstep2+3 job(s) failed - inspect with 'dx watch <job_id>'" >&2
        exit 1
    fi
    echo "  All LDstep2+3 jobs complete."
fi
