# GP2 XWAS: Notebook Pipeline for ChrX Analyses

## Overview
GP2 XWAS is a collection of notebooks for preparing and running X-chromosome wide association analyses: ingesting and annotating sample information, performing pre- and post-imputation QC, staging GWAS inputs (prinicipal components and phenotype files), running association tests, and summarizing results with chromosome-focused visualizations.

## Installation
### Clone the repository:

````
git clone https://github.com/nvk23/GP2_XWAS.git

cd GP2_XWAS
````

### [Optional] Create a Conda environment:

````
conda create -n "gp2_xwas" python=3.11

conda activate gp2_xwas
````

### Install the required Python packages:

````
pip install -r requirements.txt
````

### [Optional] Install the R packages:

````
Rscript -e 'install.packages(grep("^[^#]", readLines("requirements_R.txt"), value = TRUE))'
````
````
Rscript -e 'remotes::install_github("boxiangliu/locuscomparer")'
````

## Running the Pipeline

- `00_pre_imputation.ipynb`: prepare inputs and sample annotations ahead of imputation.
- `01_post_imputation.ipynb`: converts imputed VCFs to PLINK2 files and perform post-imputation checks and filtering.
- `02_gwas_prep.ipynb`: stage files for association testing.
- `03_run_gwas.ipynb`: execute the XWAS run.
- `04_gwas_results.ipynb`: review summary statistics and visualizations.
- `05_run_METAL.ipynb`: meta-analyze the per-ancestry results with METAL and annotate the output with LocusZoom.

The notebooks in `extra/` are optional follow-ups, each independent of the others:

- `COJO.ipynb`: conditional and joint analysis of the PAR1 signal with GCTA.
- `forest_plots.ipynb`: per-study forest plots for variants of interest.
- `sex_het.ipynb`: test male vs female effect-size heterogeneity.
- `wgs_concordance.ipynb`: NBA vs WGS genotype concordance on chrX.

## Plotting
`plots/plotting.py` generates chrX Miami or Manhattan plots directly from PLINK2 `.glm`
output or METAL meta-analyses, either as a library or from the command line.

The `.R` scripts in `plots/` are standalone figure scripts, and need the R packages
listed in `requirements_R.txt`.

- `manhattan.R`: chrX Manhattan plot of the sex-combined results, with nearest-gene labels.
- `miami.R`: the same, as a male vs female Miami plot.
- `locusCompare_setAxes.R`: LocusCompare scatter plus per-sex LocusZoom panels for one locus.
  Needs an [LDlink API token](https://ldlink.nih.gov/?tab=apiaccess) exported as `LDLINK_TOKEN`.

`plots/ref/nearest_genes.csv` maps each labelled rsID to its nearest gene and whether the
signal is novel or previously identified; edit it to add or correct an entry.

## Project Structure
```
GP2_XWAS/
├── 00_pre_imputation.ipynb
├── 01_post_imputation.ipynb
├── 02_gwas_prep.ipynb
├── 03_run_gwas.ipynb
├── 04_gwas_results.ipynb
├── 05_run_METAL.ipynb
├── extra/
│   ├── COJO.ipynb
│   ├── forest_plots.ipynb
│   ├── sex_het.ipynb
│   ├── wgs_concordance.ipynb
├── plots/
│   ├── locusCompare_setAxes.R
│   ├── manhattan.R
│   ├── miami.R
│   ├── plotting.py
│   └── ref/
│       └── nearest_genes.csv
├── src/
│   └── run_metal_plink2.py
├── requirements.txt
├── requirements_R.txt
└── README.md
```

## Software
|               Software              |      Version(s)     |                       Resource URL                       |       RRID      |                                               Notes                                               |
|:-----------------------------------:|:-------------------:|:--------------------------------------------------------:|:---------------:|:-------------------------------------------------------------------------------------------------:|
|               Python Programming Language              |      3.11     |        http://www.python.org/        | RRID:SCR_008394 |                Notebooks and plotting script run in Python; install dependencies via requirements.txt                |
|               R Project for Statistical Computing              |      4.x     |        https://www.r-project.org/        | RRID:SCR_001905 |                Figure scripts in plots/ run in R; install dependencies via requirements_R.txt                |
|               PLINK              |      1.9     |        https://www.cog-genomics.org/plink/1.9/        |        RRID:SCR_001757        |                Used for merging functionality            |
|               PLINK2              |      2.0     |        https://www.cog-genomics.org/plink/2.0/        |        RRID:SCR_001757        |                Used for core chromosome preprocessing and association testing              |
|               imputationbot              |      2.0.0     |        https://github.com/lukfor/imputationbot        |        -        |                CLI helper for submitting and accessing TOPMed Imputation Server jobs                |
|               METAL (random-metal build)              |      0.1.0     |        https://github.com/explodecomputer/random-metal        |        RRID:SCR_002013        |                Runs the meta-analysis in a GP2 container that builds the random-metal GitHub fork. This repo adds random-effects modeling to the fixed-effects analysis in standalone METAL (https://csg.sph.umich.edu/abecasis/metal/) |
|               LocusZoom              |      -     |        https://my.locuszoom.org/        |        RRID:SCR_021374        |                Web service used for regional association plots and provides annotated exports with rsIDs for GWAS results                |
|               GCTA-COJO              |      1.94.3     |        https://yanglab.westlake.edu.cn/software/gcta/        |        -        |                Used for conditional and joint analysis of PAR region                |
|               LDlink              |      -     |        https://ldlink.nih.gov/        |        -        |                Queried by plots/locusCompare_setAxes.R through the LDlinkR package for r^2 proxies around the lead variant; needs a free API token exported as LDLINK_TOKEN. Its 1000 Genomes panel is autosome-only, so the script falls back to the Ensembl REST LD endpoint for chrX                |
