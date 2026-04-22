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

**Height GWAS example**

28. **Set up height GWAS example** (`setup_height_gwas.sh` + `setup_height_gwas.py`)
    Prepares REGENIE input files for a continuous-trait GWAS of height. Intersects the European training cohort (`final_train_iids.txt` from Step 23) with the genetic sex file (`sex_covar.txt` from Step 27) to exclude individuals with sex chromosome aneuploidy. Extracts UK Biobank fields 50 (standing height) and 21003 (age at assessment) across all assessment instances, pairs measurements by instance, drops heights below 140 cm, and computes median height and mean age per individual. Centers covariates: age is mean-centered, sex is centered at 0.5 (female = -0.5, male = 0.5), and an age-by-sex interaction term is computed. Merges with PC scores 1-10 from Step 26. Outputs `phen.txt`, `covar.txt`, `training_iids.txt`, and a processing log to `$DX_OUTPUT_DIR/regenie_input/height_example/`. See [REGENIE input files](#regenie-input-files) for file format details.

29. **Run height GWAS example** (`run_continuous_regenie_gwas.sh`)
    Launches a continuous-trait REGENIE GWAS using the height example input files from Step 28. This step demonstrates the standalone GWAS runner tool described in [Running a GWAS](#running-a-gwas). Submits the REGENIE app on DNAnexus with default settings (RINT enabled, block sizes 1000/200) and writes results to `$DX_OUTPUT_DIR/regenie_output/height_example/`.

## Optional: Prepare inputs for LD-matrix generation

After `get_genotypes.sh` completes, a second pipeline — [`generate_ld.sh`](generate_ld.sh) — prepares the inputs needed to build a custom SBayesRC LD reference from UK Biobank WGS data. This is an **optional** companion pipeline; it is not required for running a GWAS with REGENIE. It is idempotent and follows the same master-orchestrator + sub-scripts pattern as `get_genotypes.sh`.

Why build a custom LD reference rather than use the published SBayesRC reference panel? Three reasons: (1) subtle hg19↔hg38 drift remains even after careful liftover; (2) the published panel is HRC-imputed, which is lower-quality than TOPmed and far lower-quality than WGS; (3) SBayesRC's own GWAS-beta imputation step performs poorly on SNPs we already drop for large allele-frequency mismatches. Building a UKBB-native, hg38, WGS-based LD reference avoids all three.

### What it does

The pipeline produces three sets of outputs:

- **A QC-filtered SNP list** (Steps 1–3) — compares alternate-allele frequencies between three sources on the same set of EUR-ancestry individuals (intersection of `fit_pca_iids.txt` from Step 24 with positive IIDs from the imputed data), then drops any SNP that fails any of three thresholds (see below).
- **A 40k-individual LD-reference cohort** (Step 4) — a random sample of unrelated European, genetically-sexed, White British UK Biobank participants.
- **An hg38 block-boundary file and per-chromosome WGS bfiles** (Steps 5–6) — the block file derived directly from the per-SNP `Block` assignments in `sbayesrc_liftover_results.csv`, and per-chromosome bfiles filtered to the 40k cohort × QC-passed SNPs. Together these are the direct inputs to SBayesRC's `LDstep1`–`LDstep4` R functions.

The three SNP-QC thresholds applied in Steps 1–3:

1. **Low minor-allele frequency in WGS** (`MAF_THRESHOLD`, default `0.009`) — rare-variant noise.
2. **WGS vs TopMed-imputed allele-frequency disagreement** (`FREQ_DIFF_THRESHOLD`, default `0.025`) — large |WGS − TopMed| indicates the UKB TopMed imputation disagrees with WGS truth in our own EUR samples.
3. **WGS vs HRC-imputed (SBayesRC) allele-frequency disagreement** (`SBAYESRC_FREQ_DIFF_THRESHOLD`, default `0.03`) — large |WGS − HRC| indicates our WGS frequencies disagree with the SBayesRC paper's own HRC-imputed white-British reference panel (hg19; lifted to hg38 and allele-aligned). The SBayesRC paper is [Zheng et al. 2024, *Nat Genet*](https://www.nature.com/articles/s41588-024-01704-y).

All three thresholds, the LD-reference cohort size (`LD_COHORT_SIZE`, default `40000`), and the sampling seed (`RANDOM_SEED`, default `0`) are exposed as modifiable constants at the top of [`generate_ld.sh`](generate_ld.sh).

### Pipeline steps

1. **Compute WGS vs TopMed allele-frequency comparison** (`compute_freq_compare.sh`)
   Submits a Swiss Army Knife job that runs plink2 `--freq counts` on the per-chromosome WGS and imputed pfiles (both restricted to the same EUR sample intersection and to variants in the merged pfiles), then merges the per-chromosome `.acount` outputs into a single per-variant CSV. Outputs `wgs_vs_imputed_freq.csv` to `$DX_OUTPUT_DIR/freq_compare/` on DNAnexus.

2. **Download frequency-comparison CSV** (`download_freq_compare.sh`)
   Downloads the DNAnexus CSV from Step 1 to `data/freq_compare/wgs_vs_imputed_freq.csv`. This is variant-level summary data (per-variant allele counts and sample totals), not individual-level.

3. **QC-filter variants** (`qc_snps.sh` + `qc_snps.py`)
   Downloads the SBayesRC liftover CSV (from the [sbayesrc-liftover](https://github.com/jesseICR/sbayesrc-liftover) release, contains the hg19 → hg38 allele mapping with HRC-imputed `A1Freq`) and runs Python QC locally. Aligns SBayesRC's hg19 `A1Freq` onto the hg38 ALT allele (handling palindromic/strand-flipped SNPs via the Watson–Crick complement), applies the three filters, and writes the annotated QC-passed CSV.

4. **Sample 40k LD-reference cohort** (`sample_ld_cohort.sh` + `sample_ld_cohort.py`)
   Submits a Swiss Army Knife job that intersects `fit_pca_iids.txt` (Step 24 of `get_genotypes.sh`), `sex_covar.txt` (Step 27), and White British (UKBB field 22006 == 1), then random-samples `LD_COHORT_SIZE` individuals with a fixed seed for reproducibility. All sampled IIDs are WGS-sequenced (`fit_pca_iids.txt` already excludes imputed-only individuals). Outputs `ld_ref_40k_iids.txt` (two-column FID/IID) and `ld_ref_cohort_log.txt` to `$DX_OUTPUT_DIR/ld_reference/`.

5. **Derive hg38 LD-block boundaries** (`build_hg38_blocks.sh` + `build_hg38_blocks.py`)
   The `sbayesrc_liftover_results.csv` downloaded in Step 3 already contains a per-SNP `Block` id (SBayesRC's original 4 cM block partition) alongside each SNP's `pos_hg38`. For each `Block` id we compute `min_i = min(pos_hg38)` and `max_i = max(pos_hg38)`, then set `StartBP_i = min_i − 1` and `EndBP_i = max_i + 1`, each decremented/incremented further while the bound coincides with any SBayesRC panel SNP position (defensive against dense SNP regions where adjacent positions differ by 1 bp). This yields `StartBP_i < pos < EndBP_i` strictly for every SNP in block `i`, with no SNP ever landing on any boundary — regardless of which half-open convention SBayesRC's `LDstep1` uses internally. Rows with `pos_hg38 == -1` (1,771 unlifted SBayesRC SNPs) are skipped. No interval liftOver is performed; this trivially achieves 100 % coverage of every QC-passed SNP by some block. All 591 original blocks are retained, with 0 inter-block overlaps. Writes `data/ld_reference/ref4cM_hg38.pos` locally and uploads it to `$DX_OUTPUT_DIR/ld_reference/ref4cM_hg38.pos`.

6. **Build per-chromosome LD-reference bfiles** (`build_ld_ref_bfiles.sh`)
   Extracts the `variant_id` column from the QC-passed CSV and uploads it as `$DX_OUTPUT_DIR/ld_reference/ld_ref_snps.txt`, then submits 22 parallel Swiss Army Knife jobs. Each job runs `plink2 --pfile wgs_pfiles/chr{N} --keep ld_ref_40k_iids.txt --extract ld_ref_snps.txt --make-bed`, producing per-chromosome bfiles at `$DX_OUTPUT_DIR/ld_reference/bfiles/chr{1..22}.{bed,bim,fam,log}`. These bfiles, together with the hg38 block file, are the direct inputs to SBayesRC's `LDstep1`–`LDstep4`.

### What it produces

- **`data/freq_compare/wgs_vs_imputed_freq_qc_passed.csv`** — one row per QC-passed variant, with the original freq-comparison columns plus four derived columns: `alt_freq_wgs`, `alt_freq_topmed`, `abs_diff_wgs_vs_topmed`, and `alt_freq_hrc` (SBayesRC HRC-imputed ALT frequency in hg38). The `variant_id` column (RSIDs) is the SNP list for the LD reference.
- **`$DX_OUTPUT_DIR/ld_reference/ld_ref_40k_iids.txt`** — the 40k LD-reference cohort (FID/IID, sorted numerically).
- **`$DX_OUTPUT_DIR/ld_reference/ref4cM_hg38.pos`** — the hg38-lifted block-boundary file (also cached locally at `data/ld_reference/ref4cM_hg38.pos`).
- **`$DX_OUTPUT_DIR/ld_reference/bfiles/chr{1..22}.{bed,bim,fam,log}`** — per-chromosome WGS bfiles filtered to the 40k cohort × QC-passed SNPs.

### Variant counts (default thresholds)

| Stage | SNPs | Loss |
|---|---:|---:|
| Original SBayesRC panel (hg19 → hg38 liftover) | 7,356,518 | — |
| After WGS + TopMed extraction + merge (`get_genotypes.sh` outputs) | 7,349,366 | 7,152 |
| After QC filtering (`generate_ld.sh` output) | 7,346,329 | 3,037 |
| **Total lost from the original SBayesRC panel** | | **10,189** |

Of the 3,037 variants dropped by QC: 333 by the MAF filter, 1,766 by the WGS-vs-TopMed freq-diff filter, 1,462 by the WGS-vs-HRC freq-diff filter (with 447 variants failing both freq-diff filters jointly and 75 failing both the MAF and WGS-vs-TopMed filters jointly).

### Running it

```bash
source ~/venvs/dnanexus/bin/activate   # if not already active
bash generate_ld.sh
```

DNAnexus writes happen in Steps 1, 4, 5 (upload only), and 6. Steps 2, 3, and the local portion of Step 5 are local only. A full fresh run on the default thresholds takes ~30–45 minutes end-to-end (dominated by the Step 1 frequency-comparison job and the Step 6 per-chromosome bfile jobs); subsequent re-runs skip every step.

## Running a GWAS

After the pipeline completes, the `run_continuous_regenie_gwas.sh` script provides a command-line tool for running continuous-trait GWAS via the [REGENIE](https://rgcgithub.github.io/regenie/) app on DNAnexus. The pipeline produces all the required genotype files (Steps 12 and 14) and includes a worked example for height (Steps 28-29), but the tool is designed to be used repeatedly for any phenotype of interest.

### Usage

```bash
bash run_continuous_regenie_gwas.sh <input_name> <output_name> [OPTIONS]
```

| Argument | Description |
|---|---|
| `input_name` | Directory name under `$DX_OUTPUT_DIR/regenie_input/` containing the three required input files |
| `output_name` | Directory name for `$DX_OUTPUT_DIR/regenie_output/` where REGENIE results will be written |

| Option | Default | Description |
|---|---|---|
| `--apply-rint` | enabled | Apply rank inverse normal transformation to the phenotype |
| `--no-apply-rint` | | Disable RINT |
| `--step1-block-size N` | 1000 | Number of variants per block in REGENIE Step 1 (whole-genome ridge regression) |
| `--step2-block-size N` | 200 | Number of variants per block in REGENIE Step 2 (association testing) |
| `--priority LEVEL` | normal | DNAnexus job priority: `low`, `normal`, or `high` |

Example:

```bash
# Run the height example produced by the pipeline
bash run_continuous_regenie_gwas.sh height_example height_v1

# Same phenotype, different settings
bash run_continuous_regenie_gwas.sh height_example height_v2 --no-apply-rint --priority high
```

### What the script does

1. **Validates inputs.** Checks that the input directory exists on DNAnexus and contains the three required files (`phen.txt`, `covar.txt`, `training_iids.txt`). Also checks that the Step 1 bfiles (`$DX_OUTPUT_DIR/direct_bfile/chr1_22_merged.*`) and Step 2 BGENs (`$DX_OUTPUT_DIR/merged_bgens/chr{1..22}.*`) exist.

2. **Overwrite protection.** Checks that the output directory does **not** already exist on DNAnexus. If it does, the script exits with an error to prevent overwriting existing GWAS results. To re-run a GWAS, choose a different output name or remove the existing output directory first.

3. **Submits the REGENIE app.** Launches the DNAnexus REGENIE app with the validated inputs, configured options, and the following fixed settings:
   - Quantitative traits (`-iquant_traits=true`)
   - PRS mode disabled (`-iprs_mode=false`)
   - Additive test (`-itest_type=additive`)
   - First allele as reference for both steps (`-istep1_ref_first=true`, `-istep2_ref_first=true`)
   - `training_iids.txt` is used as the keep file for both Step 1 and Step 2, restricting the analysis to the specified sample

4. **Reports the job ID** and monitoring commands.

### REGENIE input files

To run a GWAS, you must prepare three input files in a directory under `$DX_OUTPUT_DIR/regenie_input/` on DNAnexus. The pipeline's height example (Step 28) demonstrates how to produce these files; you can follow the same pattern for any phenotype.

#### `training_iids.txt` — sample inclusion list

A two-column, space-separated file with no header. Each row is an FID/IID pair identifying an individual to include in the GWAS. This file is passed to REGENIE's `--keep` flag for both Step 1 and Step 2, so only individuals listed here are analyzed.

```
1 1
2 2
3 3
```

This file should be the **intersection** of two pipeline outputs:

- **`$DX_OUTPUT_DIR/train_test/final_train_iids.txt`** (Step 23) — the European training cohort. This is the definitive set of European-ancestry individuals selected for GWAS, with kinship separation from the test cohort enforced.

- **`$DX_OUTPUT_DIR/genetic_sex/sex_covar.txt`** (Step 27) — the genetic sex covariate file, which **excludes individuals with sex chromosome aneuploidy**. Aneuploidy individuals are removed during Step 27 because sex chromosome aneuploidy can confound genetic analyses. However, the training cohort from Step 23 is built before aneuploidy exclusion and may still contain these individuals.

The intersection ensures that training IIDs are both (a) members of the European training cohort and (b) free of sex chromosome aneuploidy. Any individual present in `final_train_iids.txt` but absent from `sex_covar.txt` has aneuploidy and is excluded from the GWAS. In the height example, this intersection is performed automatically by `setup_height_gwas.py`.

#### `phen.txt` — phenotype file

A tab-separated file with a header row. The first two columns must be `FID` and `IID`. Subsequent columns are phenotype values (one column per trait). Only individuals with non-missing phenotype values should be included.

> **⚠️ Disclaimer:** All example values shown below are entirely synthetic and do not correspond to any UK Biobank participant. These are random illustrative values only.

```
FID	IID	height
1	1	170.0
2	2	165.5
```

In the height example, `height` is standing height in centimeters (UK Biobank field 50). When multiple assessment instances are available for an individual, measurements below 140 cm are excluded as likely errors, and the median of the remaining measurements is used.

#### `covar.txt` — covariate file

A tab-separated file with a header row. The first two columns must be `FID` and `IID`. Subsequent columns are covariates to adjust for in the association model.

```
FID	IID	age_c	sex_c	age_c_sex_c_inter	PC1_AVG	PC2_AVG	PC3_AVG	...	PC10_AVG
1	1	-5.16	-0.5	2.58	0.001	-0.005	0.002	...	0.001
2	2	3.84	0.5	1.92	-0.002	0.007	-0.001	...	-0.001
```

In the height example, the covariates are:

| Column | Description |
|---|---|
| `age_c` | Age at assessment, **mean-centered** (age minus the sample mean age). Centering ensures the intercept is interpretable and reduces collinearity with the interaction term. |
| `sex_c` | Genetic sex, **centered at 0.5** (female = -0.5, male = 0.5). Centering at 0.5 rather than using raw 0/1 coding makes the main effects of age and sex interpretable as effects at the average level of the other variable, rather than at the reference category. |
| `age_c_sex_c_inter` | Interaction term: `age_c * sex_c`. Captures sex-dependent age effects on the phenotype. |
| `PC1_AVG` through `PC10_AVG` | The first 10 principal components from Step 26, projected onto all samples. These control for population stratification. |

### Preparing your own GWAS input files

To run a GWAS for a different phenotype, create a new directory under `$DX_OUTPUT_DIR/regenie_input/` on DNAnexus (e.g., `regenie_input/bmi_example/`) containing the three files described above. You can use `setup_height_gwas.py` as a template — the key steps are:

1. Intersect `final_train_iids.txt` with `sex_covar.txt` to get eligible training IIDs
2. Extract and clean your phenotype of interest from UK Biobank fields
3. Build covariates with centering and PC scores
4. Upload the three files to your input directory on DNAnexus

Then run:

```bash
bash run_continuous_regenie_gwas.sh your_input_name your_output_name
```

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
| `DX_PRIORITY` | `normal` | DNAnexus job priority (`low`, `normal`, or `high`). Set to `normal` to save money at the cost of longer queue times. |

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
| `setup_height_gwas.sh` | Submits SAK job to prepare height GWAS input files |
| `setup_height_gwas.py` | Builds phenotype, covariate, and training sample files for height GWAS |
| `run_continuous_regenie_gwas.sh` | Standalone CLI tool to launch continuous-trait GWAS via REGENIE on DNAnexus |
| `generate_ld.sh` | Optional orchestrator — QCs SNPs for downstream LD-matrix generation (runs after `get_genotypes.sh`) |
| `compute_freq_compare.sh` | Submits SAK job to compute WGS vs TopMed allele-frequency comparison on EUR samples |
| `download_freq_compare.sh` | Downloads the freq-comparison CSV from DNAnexus to `data/freq_compare/` |
| `qc_snps.sh` | Wrapper that caches the SBayesRC liftover CSV and invokes `qc_snps.py` |
| `qc_snps.py` | Applies MAF + WGS-vs-TopMed + WGS-vs-HRC freq-diff filters to produce the QC-passed SNP list |
| `sample_ld_cohort.sh` | Submits SAK job to sample the 40k LD-reference cohort |
| `sample_ld_cohort.py` | Intersects `fit_pca_iids`, `sex_covar`, White British (UKBB 22006); random-samples with fixed seed |
| `build_hg38_blocks.sh` | Wrapper that invokes `build_hg38_blocks.py` and uploads the resulting `.pos` file to DNAnexus |
| `build_hg38_blocks.py` | Derives hg38 block boundaries from the per-SNP `Block` column of `sbayesrc_liftover_results.csv` |
| `build_ld_ref_bfiles.sh` | Uploads QC-passed SNP list, submits 22 parallel SAK jobs to build per-chrom WGS bfiles for the LD reference |
