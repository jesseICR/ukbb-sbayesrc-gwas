#!/bin/bash
# compute_freq_compare.sh — One-off QC: compare alternate-allele frequencies
# between WGS (DRAGEN) and imputed (TopMed) pfiles for variants in the merged
# pfiles, computed on the same set of EUR-ancestry individuals present in
# BOTH WGS and imputed pfiles.
#
# Sample set: intersection of pca_eur/fit_pca_iids.txt (all in WGS by
# construction) with positive IIDs from imputed_pfiles/chr1.psam.
# Variant set: IDs in merged_pfiles/chr{1..22}.pvar.
#
# Output:
#   sbayesrc_genotypes/freq_compare/wgs_vs_imputed_freq.csv   (final, per-variant)
#   sbayesrc_genotypes/freq_compare/scrap/
#       intersect_iids.txt
#       {wgs,imp}_chr{1..22}.acount      (raw plink2 allele counts)
#       {wgs,imp}_chr{1..22}.log         (plink2 logs)

set -euo pipefail

DX_WGS_PFILE_DIR="/sbayesrc_genotypes/wgs_pfiles"
DX_IMPUTED_PFILE_DIR="/sbayesrc_genotypes/imputed_pfiles"
DX_MERGED_PFILE_DIR="/sbayesrc_genotypes/merged_pfiles"
DX_PCA_EUR_DIR="/sbayesrc_genotypes/pca_eur"
DX_OUT_DIR="/sbayesrc_genotypes/freq_compare"
DX_SCRAP_DIR="${DX_OUT_DIR}/scrap"
INSTANCE_TYPE="mem2_ssd1_v2_x16"
DX_PRIORITY="${DX_PRIORITY:-normal}"

OUTPUT_CSV="wgs_vs_imputed_freq.csv"

# Idempotency: skip if the final CSV already exists
if dx ls "${DX_OUT_DIR}/${OUTPUT_CSV}" &>/dev/null; then
    echo "${OUTPUT_CSV} already exists in ${DX_OUT_DIR} — skipping"
    exit 0
fi

dx mkdir -p "${DX_SCRAP_DIR}"

# Build the SAK command. Outer heredoc is UNQUOTED so ${DX_*} expand at
# wrapper time; runtime bash vars ($chrom, $EUR_KEEP, $(seq ...), $2, $3)
# are escaped with \$ to pass through literally. Inner python heredoc
# uses a QUOTED delimiter ('PYEOF') so SAK's shell leaves it alone.
cmd=$(cat <<CMDEOF
set -eo pipefail

EUR_KEEP='/mnt/project${DX_PCA_EUR_DIR}/fit_pca_iids.txt'
IMP_PSAM='/mnt/project${DX_IMPUTED_PFILE_DIR}/chr1.psam'
WGS_DIR='/mnt/project${DX_WGS_PFILE_DIR}'
IMP_DIR='/mnt/project${DX_IMPUTED_PFILE_DIR}'
MRG_DIR='/mnt/project${DX_MERGED_PFILE_DIR}'

echo '=== WGS vs Imputed allele-frequency comparison (EUR samples) ==='
echo "  EUR keep:   \$EUR_KEEP"
echo "  WGS pfiles: \$WGS_DIR"
echo "  Imp pfiles: \$IMP_DIR"
echo "  Mrg pfiles: \$MRG_DIR"

# --- Build intersection IID set: fit_pca (EUR, all in WGS) ∩ imputed positives ---
echo ''
echo '--- Building intersection IID set ---'
awk '\$2 > 0 {print \$2}' "\$EUR_KEEP" | sort -u > eur_iids.tmp
tail -n +2 "\$IMP_PSAM" | awk '\$2 > 0 {print \$2}' | sort -u > imp_iids.tmp
comm -12 eur_iids.tmp imp_iids.tmp | awk '{print \$1, \$1}' > intersect_iids.txt
echo "  EUR fit_pca (positive):  \$(wc -l < eur_iids.tmp)"
echo "  Imputed (positive):      \$(wc -l < imp_iids.tmp)"
echo "  Intersection:            \$(wc -l < intersect_iids.txt)"
rm -f eur_iids.tmp imp_iids.tmp

# --- Per-chromosome: extract merged variant IDs, run plink2 --freq counts ---
for chrom in \$(seq 1 22); do
    echo ''
    echo "--- chr\${chrom} ---"
    awk '/^#/ {next} {print \$3}' "\$MRG_DIR/chr\${chrom}.pvar" > variants.txt
    echo "  merged variants: \$(wc -l < variants.txt)"

    plink2 --pfile "\$WGS_DIR/chr\${chrom}" --keep intersect_iids.txt --extract variants.txt --freq counts --out "wgs_chr\${chrom}"
    plink2 --pfile "\$IMP_DIR/chr\${chrom}" --keep intersect_iids.txt --extract variants.txt --freq counts --out "imp_chr\${chrom}"

    rm -f variants.txt
done

# --- Merge WGS + imputed .acount files across all chromosomes into one CSV ---
echo ''
echo '--- Merging .acount files into ${OUTPUT_CSV} ---'
cat > merge_freqs.py <<'PYEOF'
import csv
import sys


def read_acount(path):
    out = {}
    with open(path) as f:
        header = f.readline().rstrip('\n').lstrip('#').split('\t')
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            vid = parts[idx['ID']]
            out[vid] = (
                parts[idx['REF']],
                parts[idx['ALT']],
                parts[idx['ALT_CTS']],
                parts[idx['OBS_CT']],
            )
    return out


total = 0
with open('${OUTPUT_CSV}', 'w', newline='') as fout:
    w = csv.writer(fout)
    w.writerow([
        'chrom', 'variant_id', 'ref', 'alt',
        'wgs_alt_cts', 'wgs_obs_ct',
        'imp_alt_cts', 'imp_obs_ct',
    ])
    for chrom in range(1, 23):
        wgs = read_acount(f'wgs_chr{chrom}.acount')
        imp = read_acount(f'imp_chr{chrom}.acount')
        only_wgs = set(wgs) - set(imp)
        only_imp = set(imp) - set(wgs)
        if only_wgs or only_imp:
            print(
                f'chr{chrom} variant-set mismatch: '
                f'{len(only_wgs)} WGS-only, {len(only_imp)} imputed-only',
                file=sys.stderr,
            )
            sys.exit(1)
        for vid, (wref, walt, wac, woc) in wgs.items():
            iref, ialt, iac, ioc = imp[vid]
            if (wref, walt) != (iref, ialt):
                print(
                    f'allele mismatch at {vid} chr{chrom}: '
                    f'WGS={wref}/{walt} vs IMP={iref}/{ialt}',
                    file=sys.stderr,
                )
                sys.exit(1)
            w.writerow([chrom, vid, wref, walt, wac, woc, iac, ioc])
            total += 1

print(f'wrote {total} rows to ${OUTPUT_CSV}')
PYEOF

python3 merge_freqs.py
rm -f merge_freqs.py

echo ''
echo '--- Output summary ---'
wc -l '${OUTPUT_CSV}'
ls -lh '${OUTPUT_CSV}'
echo '=== Done ==='
CMDEOF
)

job_id=$(dx run swiss-army-knife \
    -icmd="${cmd}" \
    --destination "${DX_SCRAP_DIR}/" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${DX_PRIORITY}" \
    --name "freq_compare_eur_wgs_vs_imputed" \
    --ignore-reuse \
    -y --brief)

echo "Submitted ${job_id} — waiting ..."
dx wait "${job_id}"

# Move the final CSV out of scrap/ up into freq_compare/
echo "Moving ${OUTPUT_CSV} from scrap/ to ${DX_OUT_DIR}/"
dx mv "${DX_SCRAP_DIR}/${OUTPUT_CSV}" "${DX_OUT_DIR}/"

echo "Done. Final CSV at ${DX_OUT_DIR}/${OUTPUT_CSV}"
