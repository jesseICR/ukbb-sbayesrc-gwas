# UK Biobank ~500K QC-Pass Directly Genotyped SNPs

## Overview

`ukbb_500k_qc_pass_direct_snps.txt` contains 501,477 quality-controlled, directly genotyped SNP rsIDs suitable for use as variant-level covariates in whole-genome regression models such as REGENIE Step 1. This set replicates the approach described in the REGENIE publication (Mbatchou et al., *Nature Genetics*, 2021; [doi:10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7)), which used ~472K QC-filtered genotyped SNPs from UK Biobank for ridge regression in the first stage of their two-step association framework.

## File Format

Plain text, one rsID per line. No header.

```
rs2765710
rs1018436
rs2735313
...
```

## Methods: SNP Selection and Quality Control

### 1. Starting Universe of Directly Genotyped SNPs

We began with the set of directly genotyped (non-imputed) SNPs available in the UK Biobank, genotyped on the Affymetrix UK BiLEVE Axiom and UK Biobank Axiom arrays. Only SNPs that passed basic genotype rate and MAF QC thresholds were considered.

### 2. Cross-Panel Intersection

To ensure broad utility and compatibility across downstream analyses, the candidate SNP set was intersected with several independent genomic resources:

- **SBayesRC 7M panel**: SNPs were required to be present in the SBayesRC 7-million-variant reference panel (derived from UK Biobank LD reference data), ensuring compatibility with Bayesian effect-size estimation workflows.
- **UK Biobank Dragen WGS**: SNPs were required to be called in the UK Biobank whole-genome sequencing (WGS) data (Dragen pipeline), confirming their presence in sequencing-based genotype calls.
- **1000 Genomes Phase 3 WGS**: SNPs were required to be present in the 1000 Genomes Project Phase 3 whole-genome sequencing data, ensuring these variants are observable across major reference panels and facilitating downstream LD score regression and cross-cohort analyses.

### 3. Missingness Filtering

- **WGS missingness**: Variants with a missing call rate > 0.2% in the UK Biobank WGS data were removed.
- **Imputed panel missingness**: Variants with a missing call rate > 1.5% in the UK Biobank TOPMed-imputed data were removed.

### 4. Allele Frequency Filtering

- **WGS minor allele frequency (MAF)**: Variants with MAF < 0.8% in the UK Biobank WGS data were excluded.
- **Imputed panel MAF**: Variants with MAF < 0.7% in the UK Biobank TOPMed-imputed data were confirmed to already be absent from the candidate set (i.e., no additional removals needed at this threshold).

These thresholds are slightly below the conventional 1% MAF cutoff to retain variants near the boundary that pass other QC metrics, while still excluding rare variants unsuitable for common-variant regression.

### 5. Inter-Chromosomal LD Filtering

Variants involved in inter-chromosomal linkage disequilibrium (ICLD) were removed using the list from the REGENIE publication (Supplementary Table 19 of Mbatchou et al., 2021). These ICLD signals arise from probe cross-hybridization artifacts on the Affymetrix genotyping arrays and produce spurious cross-chromosome correlations that can bias ridge regression.

### 6. Allele Frequency Concordance with SBayesRC LD Reference

For each remaining SNP, the allele frequency observed in the UK Biobank WGS data was compared to the allele frequency reported in the SBayesRC 7M LD reference panel. Variants with an absolute allele frequency difference exceeding 3% were removed, as such discrepancies may indicate allele coding errors, strand ambiguity issues, or population stratification artifacts between the two resources.

### 7. Hardy-Weinberg Equilibrium Filtering

Rather than applying a fixed HWE p-value threshold (which becomes overly stringent at UK Biobank sample sizes of ~500K), we used a relative heterozygosity deviation approach. For each variant, the relative deviation was computed as |O(Het) - E(Het)| / E(Het), where O and E are the observed and expected heterozygote counts under HWE. Variants exceeding a 15% relative deviation were removed. This approach is better calibrated for large sample sizes than raw mid-p HWE tests, which reject nearly all variants at N > 100K.

### 8. LD Pruning

The remaining variants were LD-pruned using PLINK2 with the following parameters, matching the REGENIE publication:

```
plink2 --indep-pairwise 1000 100 0.9
```

- Window size: 1,000 variants
- Step size: 100 variants
- R-squared threshold: 0.9

This retains one representative variant from each group of variants in high LD, reducing redundancy while preserving genome-wide coverage.

### 9. Low-Complexity Region Removal

Variants falling within low-complexity regions (LCRs) were excluded using the LCR-hs38 BED file ([Heng Li, varcmp](https://github.com/lh3/varcmp)), which defines genomic intervals with repetitive or low-complexity sequence where genotyping and variant calling are unreliable.

### 10. Final Output

After all filtering steps, **501,477 SNPs** were retained. These are written as a simple list of rsIDs in `ukbb_500k_qc_pass_direct_snps.txt`.

## Intended Use

This SNP set is designed for:

- **REGENIE Step 1**: Whole-genome ridge regression on directly genotyped variants to build polygenic predictors, prior to single-variant association testing in Step 2.
- **General-purpose covariate construction**: Any analysis requiring a high-quality, LD-pruned, QC-filtered set of common directly genotyped variants from UK Biobank.

## Population

Quality control was performed on unrelated individuals of European ancestry in UK Biobank, consistent with the REGENIE publication. Allele frequencies and HWE statistics reflect this population.

## Reference

Mbatchou, J., Barnard, L., Backman, J. et al. Computationally efficient whole-genome regression for quantitative and binary traits. *Nat Genet* **53**, 1097-1103 (2021). [https://doi.org/10.1038/s41588-021-00870-7](https://doi.org/10.1038/s41588-021-00870-7)
