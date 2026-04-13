#!/bin/bash
# =============================================================================
# run_continuous_regenie_gwas.sh
# Command-line tool to launch a continuous-trait GWAS via the DNAnexus REGENIE
# app using genotype data produced by the get_genotypes.sh pipeline.
#
# Usage:   bash run_continuous_regenie_gwas.sh <input_name> <output_name> [OPTIONS]
# Example: bash run_continuous_regenie_gwas.sh height_example height_v1
# =============================================================================

set -euo pipefail

# =============================================================================
# Configuration — update these if your pipeline used different paths
# =============================================================================
DX_OUTPUT_DIR="/sbayesrc_genotypes"
REGENIE_APP="app-J1gf4bQ0fjQX3Gff0kzVjxkf"

# =============================================================================
# Usage
# =============================================================================
usage() {
    cat <<EOF
Usage: $(basename "$0") <input_name> <output_name> [OPTIONS]

Launch a continuous-trait GWAS on DNAnexus using the REGENIE app.

Positional arguments:
  input_name           Directory name under regenie_input/ containing the
                       required files: phen.txt, covar.txt, training_iids.txt
  output_name          Directory name for regenie_output/ where results will
                       be written (must not already exist — overwrite protection)

Options:
  --apply-rint         Apply rank inverse normal transformation (default)
  --no-apply-rint      Disable rank inverse normal transformation
  --step1-block-size N Block size for REGENIE Step 1 (default: 1000)
  --step2-block-size N Block size for REGENIE Step 2 (default: 200)
  --priority LEVEL     DNAnexus job priority: low | normal | high (default: normal)
  -h, --help           Show this help message

Examples:
  bash $(basename "$0") height_example height_v1
  bash $(basename "$0") height_example height_v1 --no-apply-rint --priority high
  bash $(basename "$0") height_example height_v1 --step1-block-size 500 --step2-block-size 400
EOF
    exit 0
}

# =============================================================================
# Defaults
# =============================================================================
APPLY_RINT=true
STEP1_BLOCK_SIZE=1000
STEP2_BLOCK_SIZE=200
PRIORITY="normal"

# =============================================================================
# Argument parsing
# =============================================================================
if [[ $# -lt 1 ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    usage
fi

if [[ $# -lt 2 ]]; then
    echo "Error: Two positional arguments required: <input_name> <output_name>"
    echo "Run '$(basename "$0") --help' for usage."
    exit 1
fi

INPUT_NAME="$1"
OUTPUT_NAME="$2"
shift 2

while [[ $# -gt 0 ]]; do
    case "$1" in
        --apply-rint)
            APPLY_RINT=true
            shift
            ;;
        --no-apply-rint)
            APPLY_RINT=false
            shift
            ;;
        --step1-block-size)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --step1-block-size requires a value"
                exit 1
            fi
            STEP1_BLOCK_SIZE="$2"
            shift 2
            ;;
        --step2-block-size)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --step2-block-size requires a value"
                exit 1
            fi
            STEP2_BLOCK_SIZE="$2"
            shift 2
            ;;
        --priority)
            if [[ -z "${2:-}" ]]; then
                echo "Error: --priority requires a value (low, normal, or high)"
                exit 1
            fi
            PRIORITY="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option '$1'"
            echo "Run '$(basename "$0") --help' for usage."
            exit 1
            ;;
    esac
done

# Validate priority
if [[ "$PRIORITY" != "low" && "$PRIORITY" != "normal" && "$PRIORITY" != "high" ]]; then
    echo "Error: --priority must be one of: low, normal, high (got '${PRIORITY}')"
    exit 1
fi

# =============================================================================
# Detect DNAnexus project
# =============================================================================
PROJECT_ID=$(dx env 2>/dev/null | grep 'Current workspace' | head -1 | awk '{print $NF}')
if [[ -z "$PROJECT_ID" ]]; then
    echo "Error: No DNAnexus project selected."
    echo "Select a project with 'dx select' before running this script."
    exit 1
fi

echo "=== REGENIE Continuous-Trait GWAS ==="
echo "DNAnexus project : ${PROJECT_ID}"
echo "Input directory  : ${DX_OUTPUT_DIR}/regenie_input/${INPUT_NAME}/"
echo "Output directory : ${DX_OUTPUT_DIR}/regenie_output/${OUTPUT_NAME}/"
echo "Apply RINT       : ${APPLY_RINT}"
echo "Step 1 block size: ${STEP1_BLOCK_SIZE}"
echo "Step 2 block size: ${STEP2_BLOCK_SIZE}"
echo "Priority         : ${PRIORITY}"
echo ""

# =============================================================================
# Validation
# =============================================================================
echo "--- Validating inputs ---"

# Check input directory exists
if ! dx ls "${DX_OUTPUT_DIR}/regenie_input/${INPUT_NAME}/" &>/dev/null; then
    echo "Error: Input directory not found: ${DX_OUTPUT_DIR}/regenie_input/${INPUT_NAME}/"
    echo "Make sure you have prepared the GWAS input files (see setup_height_gwas.sh for an example)."
    exit 1
fi

# Check required input files
REQUIRED_FILES=("phen.txt" "covar.txt" "training_iids.txt")
for f in "${REQUIRED_FILES[@]}"; do
    if ! dx ls "${DX_OUTPUT_DIR}/regenie_input/${INPUT_NAME}/${f}" &>/dev/null; then
        echo "Error: Required file not found: ${DX_OUTPUT_DIR}/regenie_input/${INPUT_NAME}/${f}"
        exit 1
    fi
done
echo "  Input files OK (phen.txt, covar.txt, training_iids.txt)"

# Overwrite protection — output directory must NOT exist
if dx ls "${DX_OUTPUT_DIR}/regenie_output/${OUTPUT_NAME}/" &>/dev/null; then
    echo "Error: Output directory already exists: ${DX_OUTPUT_DIR}/regenie_output/${OUTPUT_NAME}/"
    echo "This is a safety measure to prevent overwriting existing GWAS results."
    echo "Choose a different output name or remove the existing directory first."
    exit 1
fi
echo "  Output directory OK (does not exist yet)"

# Check Step 1 bfiles
for ext in bed bim fam; do
    if ! dx ls "${DX_OUTPUT_DIR}/direct_bfile/chr1_22_merged.${ext}" &>/dev/null; then
        echo "Error: Step 1 bfile not found: ${DX_OUTPUT_DIR}/direct_bfile/chr1_22_merged.${ext}"
        echo "Run the get_genotypes.sh pipeline first to produce the required genotype files."
        exit 1
    fi
done
echo "  Step 1 bfiles OK (chr1_22_merged.bed/bim/fam)"

# Check Step 2 bgens (chr1 as proxy for all 22)
if ! dx ls "${DX_OUTPUT_DIR}/merged_bgens/chr1.bgen" &>/dev/null; then
    echo "Error: Step 2 BGEN files not found in ${DX_OUTPUT_DIR}/merged_bgens/"
    echo "Run the get_genotypes.sh pipeline first to produce the required genotype files."
    exit 1
fi
echo "  Step 2 BGENs OK"

echo "--- Validation passed ---"
echo ""

# =============================================================================
# Create output directory
# =============================================================================
dx mkdir -p "${DX_OUTPUT_DIR}/regenie_output/${OUTPUT_NAME}"

# =============================================================================
# Build BGEN / BGI / SAMPLE arrays for chromosomes 1-22
# =============================================================================
BGEN_ARGS=()
BGI_ARGS=()
SAMPLE_ARGS=()
for i in {1..22}; do
    BGEN_ARGS+=("-igenotype_bgens=${PROJECT_ID}:${DX_OUTPUT_DIR}/merged_bgens/chr${i}.bgen")
    BGI_ARGS+=("-igenotype_bgis=${PROJECT_ID}:${DX_OUTPUT_DIR}/merged_bgens/chr${i}.bgen.bgi")
    SAMPLE_ARGS+=("-igenotype_samples=${PROJECT_ID}:${DX_OUTPUT_DIR}/merged_bgens/chr${i}.sample")
done

# =============================================================================
# Build optional RINT arguments
# =============================================================================
RINT_ARGS=()
if [[ "$APPLY_RINT" == "true" ]]; then
    RINT_ARGS+=(
        -istep1_extra_cmd_line_args="--apply-rint"
        -istep2_extra_cmd_line_args="--apply-rint"
    )
fi

# =============================================================================
# Launch REGENIE
# =============================================================================
INPUT_DIR="${DX_OUTPUT_DIR}/regenie_input/${INPUT_NAME}"
DEST="${DX_OUTPUT_DIR}/regenie_output/${OUTPUT_NAME}"
JOB_NAME="REGENIE continuous GWAS: ${OUTPUT_NAME}"

echo "Submitting REGENIE job..."

JOB_ID=$(dx run "${REGENIE_APP}" \
    -y \
    -iwgr_genotype_bed="${PROJECT_ID}:${DX_OUTPUT_DIR}/direct_bfile/chr1_22_merged.bed" \
    -iwgr_genotype_bim="${PROJECT_ID}:${DX_OUTPUT_DIR}/direct_bfile/chr1_22_merged.bim" \
    -iwgr_genotype_fam="${PROJECT_ID}:${DX_OUTPUT_DIR}/direct_bfile/chr1_22_merged.fam" \
    "${BGEN_ARGS[@]}" \
    "${BGI_ARGS[@]}" \
    "${SAMPLE_ARGS[@]}" \
    -ipheno_txt="${PROJECT_ID}:${INPUT_DIR}/phen.txt" \
    -icovar_txt="${PROJECT_ID}:${INPUT_DIR}/covar.txt" \
    -istep1_keep_txts="${PROJECT_ID}:${INPUT_DIR}/training_iids.txt" \
    -istep2_keep_txts="${PROJECT_ID}:${INPUT_DIR}/training_iids.txt" \
    "${RINT_ARGS[@]}" \
    -iquant_traits=true \
    -iprs_mode=false \
    -istep1_block_size="${STEP1_BLOCK_SIZE}" \
    -istep2_block_size="${STEP2_BLOCK_SIZE}" \
    -itest_type=additive \
    -istep1_ref_first=true \
    -istep2_ref_first=true \
    --destination "${PROJECT_ID}:${DEST}" \
    --name "${JOB_NAME}" \
    --priority "${PRIORITY}" \
    --brief)

echo ""
echo "Job submitted successfully!"
echo "  Job ID  : ${JOB_ID}"
echo "  Output  : ${PROJECT_ID}:${DEST}"
echo ""
echo "Monitor with:"
echo "  dx watch ${JOB_ID}"
echo "  dx describe ${JOB_ID}"
