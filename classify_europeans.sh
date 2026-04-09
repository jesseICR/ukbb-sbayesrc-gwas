#!/bin/bash
# classify_europeans.sh — Classify European ancestry individuals from ADMIXTURE results.
#
# Submits a Swiss Army Knife job that reads the ADMIXTURE K=6 results TSV
# and applies a threshold classifier to identify European individuals:
#   European >= 0.8, African <= 0.1, American <= 0.1,
#   East_Asian <= 0.1, Oceanian <= 0.1, South_Asian: no cap
#
# Output: ${DX_EUROPEANS_DIR}/classified_european_iids.txt
#         Two-column FID IID file (FID == IID)
#
# Expects env vars: DX_STATGEN_DIR, DX_EUROPEANS_DIR, DX_PRIORITY

set -euo pipefail

# Idempotency: skip if output already exists
if dx ls "${DX_EUROPEANS_DIR}/classified_european_iids.txt" &>/dev/null; then
    echo "  classified_european_iids.txt already exists — skipping"
    exit 0
fi

dx mkdir -p "${DX_EUROPEANS_DIR}"

tsv="/mnt/project/${DX_STATGEN_DIR}/ukb_admixture_k6.tsv"

cmd="set -eo pipefail && \
echo '--- Classifying European ancestry individuals ---' && \
total=\$(tail -n +2 '${tsv}' | wc -l) && \
echo \"Total individuals: \${total}\" && \
awk -F'\t' 'NR>1 && \$2>=0.8 && \$5<=0.1 && \$4<=0.1 && \$3<=0.1 && \$7<=0.1 {print \$1, \$1}' '${tsv}' > classified_european_iids.txt && \
n_eur=\$(wc -l < classified_european_iids.txt) && \
echo \"Classified as European: \${n_eur} (\$(awk \"BEGIN{printf \\\"%.2f\\\", 100*\${n_eur}/\${total}}\")%)\" && \
echo '--- Done ---'"

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_EUROPEANS_DIR}/" \
    --instance-type "mem1_ssd1_v2_x2" \
    --priority "${DX_PRIORITY}" \
    --name "classify_europeans" \
    -y --brief)

echo "  Submitted ${job_id} — waiting for completion ..."
dx wait "${job_id}"
echo "  Done."
