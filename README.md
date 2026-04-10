# SBayesRC Genotype Preparation for UK Biobank GWAS with REGENIE

Prepare UK Biobank genotypes for GWAS of ~7 million [SBayesRC](https://github.com/zhilizheng/SBayesRC) SNPs in hg38, using [DRAGEN](https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=185) whole-genome sequencing data on the UK Biobank Research Analysis Platform (DNAnexus). The prepared genotype files are structured for association testing with [REGENIE](https://rgcgithub.github.io/regenie/).

## Overview

This pipeline extracts ~7 million SBayesRC variant genotypes from UK Biobank [DRAGEN WGS](https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=185) data, prepares sample QC covariates (kinship, ancestry, PCA, genetic sex), and builds the REGENIE-ready genotype files — all within a single output directory on your DNAnexus project. It runs locally for lightweight data preparation, then submits jobs to DNAnexus for extraction, merging, and analysis. See [`get_genotypes.sh`](get_genotypes.sh) for the full pipeline orchestration and configuration.

The pipeline produces a training and test set of individuals genetically classified as majority European ancestry via ADMIXTURE K=6 projection. The test set is a subset of [White British](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=22006) siblings. No individuals in the training cohort are related to individuals in the testing cohort above a KING kinship coefficient of 0.0442 (the lower bound of relatedness for third-degree relatives). Individuals with [sex chromosome aneuploidy](https://biobank.ctsu.ox.ac.uk/ukb/field.cgi?id=22019) are excluded from both cohorts. Because DRAGEN WGS includes approximately 5,000 individuals who were not previously genotyped by the UK Biobank, KING kinship is recomputed as part of this pipeline using a subset of the QC SNPs that the UK Biobank used to compute relatedness ([Bycroft et al. 2018](https://www.nature.com/articles/s41586-018-0579-z)). Principal components are fit on European training-set individuals unrelated to all White British siblings, then projected onto the training and testing cohorts for use as covariates in GWAS and PRS validation, respectively.

**Runtime:** approximately 1 day end-to-end. **Cost:** approximately 15–20 GBP on DNAnexus.

The pipeline is fully reproducible and designed to run on **any** UK Biobank RAP user's project. All steps are idempotent — the pipeline can be re-run safely at any point and will skip work that has already been completed.

### What the pipeline reads

The pipeline reads from two sources on your DNAnexus project (both are standard parts of a UK Biobank dispensed dataset and are never modified):

- **`Bulk/`** — [DRAGEN WGS](https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=185) PLINK files and [TOPMed](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=21007)-imputed BGEN files
- **The dispensed dataset record** (`app*.dataset`) — UK Biobank tabular data fields (genetic sex, sex chromosome aneuploidy, White British self-report, kinship reference data)

### What the pipeline writes

All output files are written to a single directory on your DNAnexus project, controlled by `DX_OUTPUT_DIR` in [`get_genotypes.sh`](get_genotypes.sh). By default this is `/sbayesrc_genotypes`. To write output elsewhere, change `DX_OUTPUT_DIR` — all paths are derived from it automatically. No other directories are created or modified.

All DNAnexus paths are absolute (prefixed with `/`), so the pipeline works regardless of your `dx` CLI working directory.

**Important:** Your project must not already have a directory at the configured `DX_OUTPUT_DIR` path (i.e., `/mnt/project/sbayesrc_genotypes/` must not exist if using the default). If it does, rename or remove it before running the pipeline.

Throughout this document, `$DX_OUTPUT_DIR` refers to the configured output directory (default: `/sbayesrc_genotypes`).

### Merging WGS and imputed-only individuals

The primary genotype source is [DRAGEN WGS](https://biobank.ndph.ox.ac.uk/ukb/label.cgi?id=185), which covers the vast majority of UK Biobank participants. However, approximately 2,800 individuals are present in the TOPMed imputation data but absent from DRAGEN WGS. The pipeline rescues these individuals by extracting the same SBayesRC variants from the imputed data and merging them with the WGS genotypes per chromosome, so the final files include all genotyped UK Biobank participants.

### External resources

This pipeline depends on two companion repositories:

- **[sbayesrc-liftover](https://github.com/jesseICR/sbayesrc-liftover)** — Produces the alignment file `sbayesrc_hg38.csv`, which contains hg38 genomic coordinates (chrom, pos, ref, alt, rsid) for ~7.35 million SNPs from the SBayesRC SNP array panel. The liftover pipeline lifts the original hg19 SBayesRC `snp.info` to hg38 with cross-validation against dbSNP, FASTA reference checks, and 1000 Genomes allele frequency QC. This file is downloaded automatically at pipeline startup and defines the set of variants extracted from both DRAGEN WGS and TOPMed imputed data.

- **[public-statgen](https://github.com/jesseICR/public-statgen)** — Provides the reference allele frequency model for ADMIXTURE K=6 ancestry projection. The pipeline downloads a pre-computed allele frequency TSV (6-population global ADMIXTURE model: European, East Asian, American, African, South Asian, Oceanian; 135,020 SNPs) from this repository. These reference frequencies are used to project ancestry proportions onto UK Biobank participants in supervised mode, which are then used to classify European-ancestry individuals for the GWAS sample.

### Pipeline steps

**WGS (DRAGEN)**

1. **Generate DRAGEN variant IDs** (`store_dragen_ids.py`)
   Reads the alignment CSV (`data/support/sbayesrc_hg38.csv`) and writes per-chromosome variant ID files to `data/dragen_ids/chr{1..22}.txt` (git-ignored). Each ID has the format `DRAGEN:chr{chrom}:{pos}:{ref}:{alt}`.

2. **Upload DRAGEN IDs to DNAnexus** (`upload_dragen_ids.sh`)
   Uploads the 22 chromosome ID files to `$DX_OUTPUT_DIR/dragen_ids/` on DNAnexus. Skips files that already exist.

3. **Extract WGS variants from PLINK files** (`wgs_extract_variants.sh`)
   Submits 22 parallel Swiss Army Knife jobs on DNAnexus (one per chromosome). Each job runs `plink2` to extract the selected variants from the UKBB DRAGEN WGS PLINK files, applying QC filters (missingness < 3%, MAC > 1000, PASS variants only). Outputs are stored in `$DX_OUTPUT_DIR/wgs_pfiles/` on DNAnexus.

**Imputed (TOPMed)**

4. **Generate TopMed variant IDs** (`store_topmed_ids.py`)
   Reads the same alignment CSV and writes per-chromosome colon-delimited variant ID files to `data/topmed_ids/chr{1..22}.txt` (git-ignored). Each ID has the format `{chrom}:{pos}:{ref}:{alt}`.

5. **Upload TopMed IDs to DNAnexus** (`upload_topmed_ids.sh`)
   Uploads the 22 chromosome ID files to `$DX_OUTPUT_DIR/topmed_ids/` on DNAnexus. Skips files that already exist.

6. **Extract imputed variants from BGEN files** (`imputed_extract_variants.sh`)
   Submits 22 parallel Swiss Army Knife jobs on DNAnexus (one per chromosome). Each job runs `plink2` to extract the selected variants from the UKBB TOPMed-imputed BGEN files. Outputs are stored in `$DX_OUTPUT_DIR/imputed_pfiles/` on DNAnexus.

**Pvar standardization**

7. **Back up original pvar files** (`backup_pvars.sh`)
   Copies the original WGS and imputed pvar files on DNAnexus to `$DX_OUTPUT_DIR/backups/wgs_pvars/` and `$DX_OUTPUT_DIR/backups/imputed_pvars/`. Skips files that already have backups.

8. **Standardize pvar files** (`standardize_pvars.sh` + `standardize_pvar.py`)
   Downloads each pvar locally, maps variant IDs to rsids using the alignment CSVs, and trims to the 5 core columns (`#CHROM`, `POS`, `ID`, `REF`, `ALT`). Uploads the standardized pvar back to DNAnexus and deletes the local copy. Pvar files contain only variant metadata (not individual-level data), so downloading them locally is safe and consistent with DNAnexus data governance conventions.

**Merge**

9. **Find imputed-only IIDs** (`find_imputed_only_iids.sh`)
   Submits a Swiss Army Knife job that compares WGS and imputed psam files to identify individuals present in the imputed data but absent from WGS. Filters out negative (redacted) IIDs. Outputs a two-column FID/IID file to `$DX_OUTPUT_DIR/merge_steps/imputed_only_iids.txt`.

10. **Merge WGS + imputed-only individuals** (`merge_wgs_imputed.sh`)
    For each chromosome, submits a Swiss Army Knife job that: (a) finds common variants between WGS and imputed pvars by matching on ID + REF + ALT; (b) extracts a WGS bfile with common variants (all WGS samples); (c) extracts an imputed bfile with common variants (imputed-only IIDs only); (d) merges the two bfiles using plink1 `--bmerge` with trial-merge + missnp handling. Intermediate bfiles and merge logs are stored in `$DX_OUTPUT_DIR/merge_steps/bfiles/merge_chr{N}/`. Then converts each merged bfile to a pfile in `$DX_OUTPUT_DIR/merged_pfiles/`.

**QC**

11. **Validate merged pfiles** (`validate_merged_pfiles.sh`)
    Checks that all 22 merged psam files have the expected sample count (WGS + imputed-only IIDs) and that variant counts in each merged pvar match the WGS variant counts from the per-chromosome merge logs (tolerance: ≤100 fewer). Uploads a QC log to `$DX_OUTPUT_DIR/merge_steps/pfile_merge_log.txt`.

**BGEN conversion (for REGENIE step 2)**

12. **Convert merged pfiles to BGEN** (`convert_to_bgens.sh`)
    Submits 22 parallel Swiss Army Knife jobs on DNAnexus (one per chromosome). Each job runs `plink2` to export the merged pfile to BGEN v1.2 format (8-bit dosage, ref-first), then indexes the BGEN with `bgenix`. Outputs are stored in `$DX_OUTPUT_DIR/merged_bgens/` (`.bgen`, `.bgen.bgi`, `.sample`, `.log` per chromosome). These per-chromosome BGENs are used as the genetic data input for REGENIE step 2 (association testing).

**Direct SNP bfile (for REGENIE step 1)**

13. **Extract direct SNPs per chromosome** (`extract_direct_snps.sh`)
    Uploads the direct SNP list (`data/support/direct_snps/ukbb_500k_qc_pass_direct_snps.txt`, ~500k rsids) to DNAnexus, then submits 22 parallel Swiss Army Knife jobs (one per chromosome). Each job runs plink2 `--extract` to subset the merged pfile to directly genotyped SNPs. Outputs are stored in `$DX_OUTPUT_DIR/direct_pfiles/` (.pgen/.pvar/.psam per chromosome).

14. **Merge direct-SNP pfiles into bfile** (`make_direct_bfile.sh`)
    Submits a Swiss Army Knife job that uses plink2 `--pmerge-list` to merge the 22 per-chromosome direct-SNP pfiles into a single bfile at `$DX_OUTPUT_DIR/direct_bfile/chr1_22_merged` (.bed/.bim/.fam). This bfile of directly genotyped SNPs is used as the input for REGENIE step 1 (whole-genome ridge regression).

**Kinship estimation (for sample QC)**

15. **Subset direct SNPs to kinship-relevant SNPs** (`subset_kinship_snps.sh`)
    Downloads the UK Biobank SNP QC file (`ukb_snp_qc.txt`), filters to SNPs where `in_Relatedness == 1` (the SNPs used by Bycroft et al. to compute kinship), and intersects with the QC-pass direct SNP list to produce ~85,681 kinship-relevant SNPs. Uploads the subset to `$DX_OUTPUT_DIR/kinship/` on DNAnexus.

16. **Run KING kinship estimation** (`run_king_kinship.sh`)
    Submits a Swiss Army Knife job that runs plink2 `--make-king-table` on the direct bfile, extracting only the kinship-relevant SNP subset via `--extract`. Uses `--king-table-filter 0.03` to capture pairs near the standard relatedness threshold. Outputs `ukb_all_direct_rel.kin0` to `$DX_OUTPUT_DIR/kinship/`.

17. **QC kinship** (`kinship_qc.sh` + `kinship_qc.py`)
    Compares our KING kinship output against UKB's `ukb_rel.dat` reference, producing summary statistics (mean/median/percentile absolute differences, Pearson correlation, pair counts by kinship bin) and diagnostic plots (scatter, Bland-Altman) for both kinship coefficients and IBS0 values. Outputs to `$DX_OUTPUT_DIR/kinship/qc/`.

18. **Classify close relationships** (`classify_relations.sh` + `classify_relations.py`)
    Filters the KING kinship output to close relationships (kinship >= 0.1767) and classifies each pair as `sibling`, `parent_child`, or `identical` (twins/duplicates) using kinship coefficient and IBS0 thresholds from Bycroft et al. 2018. Outputs `close_relations.csv` to `$DX_OUTPUT_DIR/kinship/`.

**ADMIXTURE K=6 projection (for population stratification)**

19. **Prepare ADMIXTURE inputs** (`admixture_prep.sh` + `admixture_align_alleles.py`)
    Downloads the ADMIXTURE binary to `tools/` (if not already present) and uploads it to DNAnexus. Then submits a Swiss Army Knife job that downloads a reference allele frequency TSV (K=6 global ADMIXTURE model: European, East Asian, American, African, South Asian, Oceanian; 135,020 SNPs), extracts matching SNPs from the direct bfile, aligns alleles between the reference and UKBB genotypes (handling allele swaps and strand flips), and produces an aligned bfile + .P frequency matrix. Outputs to `$DX_OUTPUT_DIR/statgen/scrap/`.

20. **Split into ADMIXTURE batches** (`admixture_split_batches.sh`)
    Submits a Swiss Army Knife job that splits the aligned bfile into batches of 20,000 individuals. Creates per-batch bfiles and copies the .P file as `batch_NNN.K.P.in` (the naming convention ADMIXTURE requires for projection mode). Outputs to `$DX_OUTPUT_DIR/statgen/scrap/batches/`.

21. **Run ADMIXTURE projection and build results TSV** (`admixture_run_projection.sh`)
    Submits ~25 parallel Swiss Army Knife jobs (one per batch) to run ADMIXTURE in projection mode (`-P` flag), which estimates ancestry proportions for each individual while holding allele frequencies fixed to the reference model. After all batches complete, submits a concat job that concatenates the per-batch .Q files and pairs them with individual EIDs to produce the final `ukb_admixture_k6.tsv` at `$DX_OUTPUT_DIR/statgen/`.

**European ancestry classification**

22. **Classify European ancestry individuals** (`classify_europeans.sh`)
    Reads the ADMIXTURE K=6 results TSV and applies a threshold classifier to identify European individuals: European >= 0.8, African <= 0.1, American <= 0.1, East_Asian <= 0.1, Oceanian <= 0.1 (no cap on South Asian). Outputs a two-column FID/IID file to `$DX_OUTPUT_DIR/europeans/classified_european_iids.txt`.

**Train/test sample split (for REGENIE)**

23. **Build train/test sample split** (`make_train_test_samples.sh` + `make_train_test_samples.py`)
    Splits European-ancestry individuals into training and test samples for REGENIE, ensuring no cryptic relatedness (3rd degree or closer, kinship >= 0.0441) between the two sets. The test sample comprises White British (UKBB field 22006) siblings from our European ancestry classification; the training sample comprises all remaining Europeans, expanded to include European relatives up to 3rd degree, then pruned of anyone related to a test individual. Verification checks confirm that all test individuals are siblings, White British, and European, and that no test individual has a 3rd+ degree relative in training. Outputs `final_train_iids.txt` and `final_test_iids.txt` to `$DX_OUTPUT_DIR/train_test/`.

**PCA (for population stratification covariates)**

24. **Select unrelated European IIDs for fitting PCA** (`select_pca_europeans.sh` + `select_pca_europeans.py`)
    Selects unrelated European IIDs for fitting PCA. Starts with classified Europeans, identifies White British siblings and expands to all their relatives (kinship >= 0.0441), removes the expanded set and imputed-only IIDs, then applies plink2 `--king-cutoff-table` to obtain a maximal unrelated subset. Outputs `fit_pca_iids.txt` to `$DX_OUTPUT_DIR/pca_eur/`.

25. **QC SNPs for PCA** (`pca_snp_qc.sh`)
    Subsets the direct bfile to PCA-fitting IIDs, filters SNPs with MAF < 0.01, excludes long-range LD regions (Price et al. 2008, hg38), and LD prunes (window=1000, step=80, r2=0.1). Outputs `pca_ready.{bed,bim,fam}` to `$DX_OUTPUT_DIR/pca_eur/`.

26. **Fit PCA on unrelated Europeans, project onto all samples** (`fit_project_pca.sh`)
    Fits PCA on the PCA-ready bfile (20 PCs, approximate algorithm, allele weights, seed 0), computes allele frequency counts, and projects PCs onto all samples in the direct bfile. Outputs `ukb_projected.sscore` (projected PC scores for all samples) and PCA reference files (eigenvalues, eigenvectors, allele counts) to `$DX_OUTPUT_DIR/pca_eur/`.

**Genetic sex covariate**

27. **Build genetic sex covariate file** (`get_genetic_sex.sh` + `get_genetic_sex.py`)
    Extracts UK Biobank fields 22001 (genetic sex) and 22019 (sex chromosome aneuploidy), excludes individuals with aneuploidy, assigns sex from field 22001 for imputed-only IIDs and from the fam file for WGS IIDs. Outputs `sex_covar.txt` (FID/IID/sex_01, coding 0=female, 1=male) to `$DX_OUTPUT_DIR/genetic_sex/`.

## Setup

### 1. Create a Python virtual environment and install `dxpy`

The DNAnexus Platform SDK (`dxpy`) provides the `dx` command-line tool used throughout this pipeline. Install it in a virtual environment to keep your system Python clean:

```bash
python3 -m venv ~/venvs/dnanexus
source ~/venvs/dnanexus/bin/activate
pip install --upgrade pip
pip install dxpy
```

Verify with `dx --version`.

### 2. Log in to DNAnexus

You must be logged in to a UK Biobank Research Analysis Platform project before running the pipeline:

```bash
dx login
```

This will prompt for your username, password, and 2FA code. Alternatively, use a token: `dx login --token YOUR_TOKEN_HERE`.

### 3. Run the pipeline

```bash
source ~/venvs/dnanexus/bin/activate   # if not already active
bash get_genotypes.sh
```

The pipeline takes approximately 1 day to run end-to-end. If you want to free up your terminal (or keep the pipeline running after closing an SSH session), use `nohup` to run it in the background:

```bash
nohup bash get_genotypes.sh &
```

- `nohup` detaches the process from your terminal session, so it keeps running even if you close the window or disconnect from SSH.
- `&` puts the process in the background and returns your shell prompt immediately.
- When you run this, the shell prints the process ID (PID), e.g. `[1] 48293`. You can check whether the pipeline is still running with `ps -p <PID>`.
- The pipeline writes a timestamped log to `logs/run_YYYYMMDD_HHMMSS.log`. To follow progress live: `tail -f logs/run_*.log`. `nohup` also creates a `nohup.out` file with the same output — this can be ignored.
- Note: `nohup` survives terminal disconnection, but not machine sleep. If you close your laptop lid, the process pauses and network connections (e.g., `dx` commands) may time out. For long unattended runs, use a machine that stays on (e.g., an EC2 instance with `tmux`).

> **⚠️ Caution:** The pipeline automatically runs `pip install -r requirements.txt` at startup to install Python dependencies (currently just `pandas`). This is idempotent — pip does nothing if the packages are already installed. If you do not want these installed globally, make sure you have activated a virtual environment before running the pipeline.

All steps are idempotent — the pipeline can be re-run safely at any point and will skip work that has already been completed.

### Requirements summary

- **Python 3** (>= 3.8) with `dxpy` installed (provides the `dx` CLI)
- **`curl`** and internet access (to download alignment files, ADMIXTURE binary, and reference data)
- Must be **logged in** to a UKBB RAP project via `dx login` before running the pipeline

### QC parameters

Configured at the top of `get_genotypes.sh`. These apply to WGS extraction only (steps 1–3); the imputed extraction does not apply additional QC filters.

| Parameter | Value | Description |
|---|---|---|
| `DX_PRIORITY` | `high` | DNAnexus job priority (`low`, `normal`, or `high`). Set to `normal` to save money at the cost of longer queue times. |

## File descriptions

| File | Description |
|---|---|
| `get_genotypes.sh` | Orchestrator — defines config and delegates to sub-scripts |
| `store_dragen_ids.py` | Generates per-chromosome DRAGEN variant ID files |
| `upload_dragen_ids.sh` | Uploads DRAGEN ID files to DNAnexus |
| `wgs_extract_variants.sh` | Submits 22 parallel WGS plink2 extraction jobs on DNAnexus |
| `store_topmed_ids.py` | Generates per-chromosome TopMed variant ID files |
| `upload_topmed_ids.sh` | Uploads TopMed ID files to DNAnexus |
| `imputed_extract_variants.sh` | Submits 22 parallel imputed plink2 extraction jobs on DNAnexus |
| `backup_pvars.sh` | Backs up original pvar files on DNAnexus before standardization |
| `standardize_pvars.sh` | Orchestrates pvar download, transformation, and re-upload |
| `standardize_pvar.py` | Transforms a single pvar: rsid mapping + column standardization |
| `find_imputed_only_iids.sh` | Identifies individuals in imputed data but not WGS |
| `merge_wgs_imputed.sh` | Merges WGS + imputed-only individuals per chromosome and converts to pfiles |
| `validate_merged_pfiles.sh` | QC validation of merged pfiles (sample + variant count checks) |
| `convert_to_bgens.sh` | Converts merged pfiles to BGEN v1.2 format and indexes with bgenix |
| `extract_direct_snps.sh` | Extracts direct SNPs from each per-chromosome pfile on DNAnexus |
| `make_direct_bfile.sh` | Merges per-chromosome direct-SNP pfiles into one bfile on DNAnexus |
| `subset_kinship_snps.sh` | Subsets direct SNPs to kinship-relevant SNPs and uploads to DNAnexus |
| `run_king_kinship.sh` | Runs KING kinship estimation on the direct bfile on DNAnexus |
| `kinship_qc.sh` | Submits SAK job to compare KING results against UKB reference |
| `kinship_qc.py` | Produces kinship/IBS0 comparison summaries and diagnostic plots |
| `classify_relations.sh` | Submits SAK job to classify close relationships from KING output |
| `classify_relations.py` | Classifies pairs as sibling, parent_child, or identical from kinship data |
| `admixture_prep.sh` | Downloads ADMIXTURE binary, uploads it, submits SAK job for SNP extraction + allele alignment |
| `admixture_align_alleles.py` | Aligns alleles between reference TSV and extracted bim, produces .P file |
| `admixture_split_batches.sh` | Submits SAK job to split aligned bfile into per-batch bfiles |
| `admixture_run_projection.sh` | Submits per-batch ADMIXTURE SAK jobs and a final concat job |
| `classify_europeans.sh` | Classifies European ancestry individuals from ADMIXTURE results |
| `make_train_test_samples.sh` | Submits SAK job to build train/test sample split |
| `make_train_test_samples.py` | Splits Europeans into train/test with kinship separation and verification |
| `select_pca_europeans.sh` | Submits SAK job to select unrelated European IIDs for PCA fitting |
| `select_pca_europeans.py` | Filters Europeans, expands exclusions, applies KING cutoff for unrelated set |
| `pca_snp_qc.sh` | Submits SAK job for MAF filtering, LD region exclusion, and LD pruning |
| `fit_project_pca.sh` | Submits SAK job to fit PCA on unrelated Europeans and project onto all samples |
| `get_genetic_sex.sh` | Submits SAK job to build genetic sex covariate file |
| `get_genetic_sex.py` | Assigns sex from field 22001 / fam file, excludes aneuploidy individuals |
