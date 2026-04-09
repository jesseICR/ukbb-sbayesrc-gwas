#!/bin/bash
# admixture_run_projection.sh — Run ADMIXTURE K=6 projection and build results TSV.
#
# 1. Submits one SAK job per batch (~25 parallel jobs) to run ADMIXTURE
#    in projection mode (-P flag).
# 2. After all batch jobs complete, submits a concat SAK job that:
#    - Concatenates all .Q files in sorted order
#    - Validates line counts against .fam
#    - Builds the final TSV with eid + ancestry columns
#
# Output on DNAnexus:
#   ${DX_ADMIXTURE_BATCH_DIR}/batch_NNN.6.Q   (per-batch ADMIXTURE output)
#   ${DX_STATGEN_DIR}/ukb_admixture_k6.tsv    (final output)
#
# Expects env vars: DX_STATGEN_DIR, DX_ADMIXTURE_SCRAP_DIR,
#                   DX_ADMIXTURE_BATCH_DIR, ADMIXTURE_K,
#                   ADMIXTURE_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

K="${ADMIXTURE_K}"

# Idempotency: skip if final output already exists
if dx ls "${DX_STATGEN_DIR}/ukb_admixture_k6.tsv" &>/dev/null; then
    echo "  ukb_admixture_k6.tsv already exists — skipping projection"
    exit 0
fi

# Clean up any .Q files from a previous partial run to avoid duplicates
existing_q=$(dx ls "${DX_ADMIXTURE_BATCH_DIR}/" 2>/dev/null | grep '\.Q$' || true)
if [[ -n "${existing_q}" ]]; then
    echo "  Cleaning up .Q files from previous partial run ..."
    echo "${existing_q}" | while read -r f; do
        dx rm "${DX_ADMIXTURE_BATCH_DIR}/${f}" 2>/dev/null || true
    done
fi

dx mkdir -p "${DX_STATGEN_DIR}"

# ---- Phase 1: Run ADMIXTURE on each batch ----
echo "  Discovering batches ..."
batch_names=$(dx ls "${DX_ADMIXTURE_BATCH_DIR}/" 2>/dev/null \
    | grep '\.bed$' | sed 's/\.bed$//' | sort)

n_batches=$(echo "${batch_names}" | wc -l)
echo "  Found ${n_batches} batches"

admixture_bin="/mnt/project/${DX_ADMIXTURE_SCRAP_DIR}/admixture"
batch_dir="/mnt/project/${DX_ADMIXTURE_BATCH_DIR}"

submitted=0
job_ids=()

for batch in ${batch_names}; do
    cmd="set -eo pipefail && \
echo '--- ADMIXTURE projection: ${batch} ---' && \
cp '${batch_dir}/${batch}.bed' ${batch}.bed && \
cp '${batch_dir}/${batch}.bim' ${batch}.bim && \
cp '${batch_dir}/${batch}.fam' ${batch}.fam && \
cp '${batch_dir}/${batch}.${K}.P.in' ${batch}.${K}.P.in && \
cp '${admixture_bin}' admixture && chmod +x admixture && \
echo \"Individuals: \$(wc -l < ${batch}.fam)\" && \
echo \"SNPs: \$(wc -l < ${batch}.bim)\" && \
./admixture -P ${batch}.bed ${K} && \
echo \"Q file lines: \$(wc -l < ${batch}.${K}.Q)\" && \
rm -f ${batch}.bed ${batch}.bim ${batch}.fam ${batch}.${K}.P.in ${batch}.${K}.P admixture"

    job_id=$(dx run swiss-army-knife \
        -icmd="${cmd}" \
        --destination "${DX_ADMIXTURE_BATCH_DIR}/" \
        --instance-type "${ADMIXTURE_INSTANCE_TYPE}" \
        --priority "${DX_PRIORITY}" \
        --name "admixture_${batch}" \
        -y --brief)

    echo "    ${batch}: submitted ${job_id}"
    job_ids+=("${job_id}")
    submitted=$((submitted + 1))
done

echo "  Submitted ${submitted} ADMIXTURE jobs — waiting for all to complete ..."
for job_id in "${job_ids[@]}"; do
    dx wait "${job_id}"
done
echo "  All ADMIXTURE jobs complete."

# ---- Phase 2: Concatenate results into final TSV ----
echo "  Concatenating results ..."

fam="/mnt/project/${DX_ADMIXTURE_SCRAP_DIR}/ukb_admixture_aligned.fam"

cmd="set -eo pipefail && \
echo '--- Concatenating ADMIXTURE results ---' && \
cat \$(ls '${batch_dir}'/batch_*.${K}.Q | sort) > ukb_all.${K}.Q && \
n_q=\$(wc -l < ukb_all.${K}.Q) && \
n_fam=\$(wc -l < '${fam}') && \
echo \"Q lines: \${n_q}, FAM lines: \${n_fam}\" && \
if [ \"\${n_q}\" -ne \"\${n_fam}\" ]; then echo 'ERROR: line count mismatch' && exit 1; fi && \
awk '{print \$2}' '${fam}' > eids.txt && \
tr ' ' '\t' < ukb_all.${K}.Q > ukb_all.${K}.Q.tab && \
printf 'eid\tEuropean\tEast_Asian\tAmerican\tAfrican\tSouth_Asian\tOceanian\n' > ukb_admixture_k6.tsv && \
paste eids.txt ukb_all.${K}.Q.tab >> ukb_admixture_k6.tsv && \
echo \"Created ukb_admixture_k6.tsv with \${n_q} individuals\" && \
rm -f ukb_all.${K}.Q ukb_all.${K}.Q.tab eids.txt"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_STATGEN_DIR}/" \
    --instance-type "mem1_ssd1_v2_x2" \
    --priority "${DX_PRIORITY}" \
    --name "admixture_concat_results" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"
echo "  Done. Final output: ${DX_STATGEN_DIR}/ukb_admixture_k6.tsv"
