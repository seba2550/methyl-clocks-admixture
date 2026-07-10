# Methylation Clocks Applied to Genetically Diverse and Admixed Populations

This repository holds the code, configuration, and manuscript documentation for the study:  
**"Methylation Clocks Fail to Generalize Across Genetically Admixed Individuals"**

Epigenetic aging clocks based on DNA methylation patterns across the genome have emerged as a potential biomarker for risk of age-related diseases, like Alzheimer’s disease (AD), and environmental and social stressors. However, methylation clocks have not been comprehensively validated in genetically diverse individuals. Here we evaluate a set of first-, second-, and third-generation methylation clocks in 621 AD patients and matched controls from African American, Hispanic, and White cohorts. 

### Citation & Publication Reference
Our work is published in *eLife*:  
**"Methylation Clocks Fail to Generalize Across Genetically Admixed Individuals"**  
Article DOI/Link: [https://doi.org/10.7554/eLife.105343.2](https://doi.org/10.7554/eLife.105343.2)

---

## Directory Structure

The repository is organized as follows:

```
methyl-clocks-admixture/
├── README.md                  # Master repository documentation
├── environment.yml            # Conda environment package specifications
├── data/                      # Reference BED coordinates, clock coefficients, and small metadata
│   ├── DunedinPACE_from_pkg.csv
│   ├── EN_from_pkg.csv
│   ├── Hannum_from_pkg.csv
│   ├── Horvath_from_pkg.csv
│   ├── PhenoAge_from_pkg.csv
│   ├── DunedinPACE_hg38.bed
│   ├── Hannum_hg38.bed
│   ├── Horvath_hg38.bed
│   ├── PhenoAge_hg38.bed
│   ├── Zhang_hg38.bed
│   └── bio_age_estimates_magenta_age_diff_metadata.csv
├── scripts/                   # Renumbered, step-by-step analysis and plotting scripts
│   ├── 01_age_distribution.R
│   ├── 02_magenta_accuracy.R
│   ├── 03_replication_accuracy.R
│   ├── 04_magenta_age_acceleration.R
│   ├── 05_diff_methylation.R
│   ├── 06_gnomad_variants.R
│   ├── 07_meqtl_overlap_analysis.R
│   ├── 08_variant_disruption.R
│   ├── 09_meqtl_frequencies.R
│   ├── 10_pc_clocks.R
│   └── 11_assemble_figures.R
└── manuscript/                # LaTeX source document parts, library, and final figure assets
    ├── *.tex                  # Main text and section source LaTeX files
    ├── *.bib                  # BibTeX library
    ├── main_figs/             # Main text Figures 1 to 6 (PDF formats)
    └── supplement/            # Supplementary text and PDFs (S1 to S3)
```

---

## Installation & Setup

We recommend managing the project dependencies using [Conda](https://docs.conda.io/en/latest/). You can install R, Python, and all required library dependencies into a dedicated environment using `environment.yml`:

```bash
# Clone the repository
git clone https://github.com/seba2550/methyl-clocks-admixture.git
cd methyl-clocks-admixture

# Create and activate the conda environment
conda env create -f environment.yml
conda activate methyl-clocks-admixture
```

---

## Data Acquisition & Preparation

### 1. Methylation Datasets
The whole-blood DNA methylation datasets used in the analyses are publicly accessible via the NCBI Gene Expression Omnibus (GEO) portal:
- **MAGENTA Cohort:** [GSE338167](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE338167)
- **Grady Trauma Project:** [GSE72680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72680)
- **GENOA Study:** [GSE210254](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210254)
- **Swedish Cohort:** [GSE87571](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571)

### 2. Genotyping Ancestry Proportions
Global ancestry proportions calculated from genotyping data for the MAGENTA cohorts are provided in this repository under:
- `data/genotyping_data/AA_ancestry_proportions.txt`
- `data/genotyping_data/HISPANIC_ancestry_proportions.txt`

### 3. meQTL and gnomAD Variant Data
- Population-specific variants and allele frequencies can be downloaded from the [gnomAD Browser (v3.0 / v4.1)](https://gnomad.broadinstitute.org/).
- meQTL data mapping genetic associations with CpG methylation levels can be retrieved from their respective repositories (such as the GENOA meQTL summary stats and SCREEN portals).

---

## Step-by-Step Replication Guide

To reproduce the findings and figures in the manuscript, run the scripts in numerical order:

1. **`scripts/01_age_distribution.R`**
   - Extracts replication cohort age metadata (from GSE72680, GSE210254, GSE87571) and compares age distribution boxplots (generates **Supplementary Figure 1**).
2. **`scripts/02_magenta_accuracy.R`**
   - Performs clock accuracy analysis in MAGENTA controls, maps global ancestry fractions, and saves RDS figures (generates data for **Figure 2**).
3. **`scripts/03_replication_accuracy.R`**
   - Assesses clock accuracy inside replication cohorts (GENOA, Grady, Swedish) and outputs absolute and relative clock accuracy grids (generates **Figure 3** and **Figure 4**).
4. **`scripts/04_magenta_age_acceleration.R`**
   - Models residual epigenetic age acceleration by case/control status across all admixed groups (generates **Figure 5A**).
5. **`scripts/05_diff_methylation.R`**
   - Runs differential methylation mapping to find CpG sites whose methylation levels associate with prediction error (generates **Figure 6B / Panel B**).
6. **`scripts/06_gnomad_variants.R`**
   - Summarizes genomic variants disrupting clock CpGs in gnomAD (generates **Supplementary Figure 9 & 10**).
7. **`scripts/07_meqtl_overlap_analysis.R`**
   - Identifies clock CpG sites showing differential methylation in African vs European ancestries and plots meQTL overlaps (generates **Figure 6A, 6F, 6G / Supplementary Figure 13**).
8. **`scripts/08_variant_disruption.R`**
   - Models variant disruption frequency patterns across cohorts (generates **Figure 6C / Panel C**).
9. **`scripts/09_meqtl_frequencies.R`**
   - Analyzes global and local ancestry meQTL frequencies across gnomAD populations (generates **Figure 6D, 6E / Supplementary Figure 12**).
10. **`scripts/10_pc_clocks.R`**
    - Runs principal component version of clocks and validates accuracy/generalization improvements across cohorts (generates **Supplementary Figures 4, 5, & 6**).
11. **`scripts/11_assemble_figures.R`**
    - Combines generated panel RDS files into ready-to-publish figure layouts (assembling **Figure 2, 3, & 6**).
