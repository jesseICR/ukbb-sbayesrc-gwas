#!/bin/bash
# sbayesrc_sak_setup.sh — install SBayesRC (R package) + GCTB (binary) on a SAK worker.
#
# Sourced at the top of the SAK cmd in init_ld_matrices.sh, build_ld_matrices.sh,
# and merge_ld_info.sh. R 4.4.3 is preinstalled on SAK; we add:
#   - R deps (Rcpp, data.table, stringi, BH, RcppEigen) from CRAN
#   - SBayesRC v0.2.6 from the official GitHub release tarball
#   - GCTB v2.05 binary (only used by LDstep2 via SBayesRC's internal ld.sh)
#
# NOTE: do NOT put `set -u` here (see CLAUDE.md SAK gotchas). Caller runs `set -eo pipefail`.

echo "--- sbayesrc_sak_setup: start $(date -u +%FT%TZ) ---"
t_start=$(date +%s)

NCPU=$(nproc)

# ------------------------------------------------------------
# 1) Install R dependencies + SBayesRC package
# ------------------------------------------------------------
echo "  [1/2] Installing R deps and SBayesRC (Ncpus=${NCPU})..."

R --no-save <<RSCRIPT
options(Ncpus=${NCPU})
repos <- "https://cloud.r-project.org/"
deps <- c("Rcpp", "data.table", "stringi", "BH", "RcppEigen")
missing <- deps[!deps %in% rownames(installed.packages())]
if (length(missing) > 0) {
    install.packages(missing, repos=repos)
}
if (!"SBayesRC" %in% rownames(installed.packages())) {
    install.packages(
        "https://github.com/zhilizheng/SBayesRC/releases/download/v0.2.6/SBayesRC_0.2.6.tar.gz",
        repos=NULL, type="source")
}
stopifnot("SBayesRC" %in% rownames(installed.packages()))
cat("  SBayesRC installed:", as.character(packageVersion("SBayesRC")), "\n")
RSCRIPT

# ------------------------------------------------------------
# 2) Install GCTB 2.04.3 binary (LDstep2 only; harmless in other steps)
#
# The cnsgenomics.com URL 302-redirects to gctbhub.cloud.edu.au. Newer version
# numbers return 404 on the redirect target; 2.04.3 is the current shipped
# Linux build. Pass GCTB_URL to override.
# ------------------------------------------------------------
if ! command -v gctb >/dev/null; then
    echo "  [2/2] Installing GCTB 2.04.3 binary..."
    gctb_url="${GCTB_URL:-https://cnsgenomics.com/software/gctb/download/gctb_2.04.3_Linux.zip}"
    tmpdir=$(mktemp -d)
    if ! wget "${gctb_url}" -O "${tmpdir}/gctb.zip" 2>&1 | tail -5; then
        echo "ERROR: wget failed for ${gctb_url}" >&2
        exit 1
    fi
    unzip -q -o "${tmpdir}/gctb.zip" -d "${tmpdir}"
    gctb_bin=$(find "${tmpdir}" -maxdepth 3 -type f -name gctb -executable | head -1)
    if [[ -z "${gctb_bin}" ]]; then
        echo "ERROR: gctb binary not found inside ${gctb_url}" >&2
        find "${tmpdir}" -maxdepth 3 -type f | head -20
        exit 1
    fi
    cp "${gctb_bin}" /usr/local/bin/gctb
    chmod +x /usr/local/bin/gctb
    rm -rf "${tmpdir}"
else
    echo "  [2/2] GCTB already on PATH: $(command -v gctb)"
fi
echo "  gctb binary: $(command -v gctb)"

t_end=$(date +%s)
echo "--- sbayesrc_sak_setup: done in $((t_end - t_start))s ---"
