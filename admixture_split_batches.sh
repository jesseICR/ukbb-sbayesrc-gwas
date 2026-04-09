#!/bin/bash
# admixture_split_batches.sh — Split aligned bfile into batches for ADMIXTURE.
#
# Submits a SAK job that reads the aligned bfile from DNAnexus and creates
# per-batch bfiles (batch_001, batch_002, ...) with ADMIXTURE_BATCH_SIZE
# individuals each. Also copies the .P file as batch_NNN.K.P.in for each
# batch (ADMIXTURE requires <prefix>.<K>.P.in co-located with the .bed).
#
# Output on DNAnexus (in ${DX_ADMIXTURE_BATCH_DIR}/):
#   batch_NNN.{bed,bim,fam,log}, batch_NNN.K.P.in  (for each batch)
#
# Expects env vars: DX_ADMIXTURE_SCRAP_DIR, DX_ADMIXTURE_BATCH_DIR,
#                   ADMIXTURE_BATCH_SIZE, ADMIXTURE_K,
#                   ADMIXTURE_PREP_INSTANCE_TYPE, DX_PRIORITY

set -euo pipefail

# Idempotency: skip if first batch already exists
if dx ls "${DX_ADMIXTURE_BATCH_DIR}/batch_001.bed" &>/dev/null; then
    echo "  Batch files already exist — skipping split"
    exit 0
fi

dx mkdir -p "${DX_ADMIXTURE_BATCH_DIR}"

aligned="/mnt/project/${DX_ADMIXTURE_SCRAP_DIR}/ukb_admixture_aligned"
p_file="/mnt/project/${DX_ADMIXTURE_SCRAP_DIR}/ref_aligned.P"
K="${ADMIXTURE_K}"
BATCH_SIZE="${ADMIXTURE_BATCH_SIZE}"

cmd="set -eo pipefail && \
echo '--- Splitting aligned bfile into batches ---' && \
total=\$(wc -l < '${aligned}.fam') && \
n_batches=\$(( (total + ${BATCH_SIZE} - 1) / ${BATCH_SIZE} )) && \
echo \"Total individuals: \${total}, Batch size: ${BATCH_SIZE}, Batches: \${n_batches}\" && \
for i in \$(seq 1 \${n_batches}); do \
    batch_name=\$(printf 'batch_%03d' \$i) && \
    start=\$(( (i-1)*${BATCH_SIZE} + 1 )) && \
    end=\$(( i*${BATCH_SIZE} )) && \
    awk -v s=\${start} -v e=\${end} 'NR>=s && NR<=e {print \$1,\$2}' '${aligned}.fam' > keep_\${batch_name}.txt && \
    plink2 --bfile '${aligned}' --keep keep_\${batch_name}.txt --make-bed --out \${batch_name} && \
    cp '${p_file}' \${batch_name}.${K}.P.in && \
    rm -f keep_\${batch_name}.txt && \
    echo \"  Prepared \${batch_name}: \$(wc -l < \${batch_name}.fam) individuals\"; \
done && \
rm -f *.nosex && \
echo '--- Split complete ---'"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_ADMIXTURE_BATCH_DIR}/" \
    --instance-type "${ADMIXTURE_PREP_INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "admixture_split_batches" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"
echo "  Done."
