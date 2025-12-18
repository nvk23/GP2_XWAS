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

### Install the required packages:

````
pip install -r requirements.txt
````

## Running the Pipeline
Open the root notebooks and run them in order to walk through the workflow:

- `00_pre_imputation.ipynb`: prepare inputs and sample annotations ahead of imputation.
- `01_post_imputation.ipynb`: converts imputed VCFs to PLINK2 files and perform post-imputation checks and filtering.
- `02_gwas_prep.ipynb`: stage files for association testing.
- `03_run_gwas.ipynb`: execute the GWAS/XWAS run.
- `04_gwas_results.ipynb`: review summary statistics and visualization-ready outputs.

Use modules in `src/plotting.py` to generate chrX Miami or Manhattan plots directly from PLINK2 `.glm` outputs.

## Project Structure
```
GP2_XWAS/
├── 00_pre_imputation.ipynb
├── 01_post_imputation.ipynb
├── 02_gwas_prep.ipynb
├── 03_run_gwas.ipynb
├── 04_gwas_results.ipynb
├── src/
│   └── plotting.py
├── requirements.txt
└── README.md
```

## Software
|               Software              |      Version(s)     |                       Resource URL                       |       RRID      |                                               Notes                                               |
|:-----------------------------------:|:-------------------:|:--------------------------------------------------------:|:---------------:|:-------------------------------------------------------------------------------------------------:|
|               Python Programming Language              |      3.11     |        http://www.python.org/        | RRID:SCR_008394 |                Notebooks and plotting script run in Python; install dependencies via requirements.txt                |
|               PLINK              |      1.9     |        https://www.cog-genomics.org/plink/1.9/        |        RRID:SCR_001757        |                Used for merging functionality            |
|               PLINK2              |      2.0     |        https://www.cog-genomics.org/plink/2.0/        |        RRID:SCR_001757        |                Used for core chromosome preprocessing and association testing              |
|               imputationbot              |      2.0.0     |        https://github.com/lukfor/imputationbot        |        -        |                CLI helper for submitting and accessing TOPMed Imputation Server jobs                |
