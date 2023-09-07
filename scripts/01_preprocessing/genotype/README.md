# Processing Genotypes

## Table of Contents

- [Step 1 VCF LiftOver Pipeline](#vcf-liftover-pipeline)
- [Step 2 PLINK2 VCF Processing Script Instructions](#plink2-vcf-processing-script-instructions)
- [Step 3 Running Mishap Test with PLINK](#running-mishap-test-with-plink)
- [Step 4 Pre-processing for ROSMAP QTL Analysis](#pre-processing-for-rosmap-qtl-analysis)
- [Step 5 Running PLINK and PLINK2 for data formatting](#running-plink-and-plink2-for-data-formatting)
- [Step 6 Principal Component Analysis (PCA) using Eigensoft SmartPCA](#principal-component-analysis-pca-using-eigensoft-smartpca)
- [Step 7 Formatting and Analyzing Eigenstrat Output](#formatting-and-analyzing-eigenstrat-output)
- [Step 8 Removing Population Outliers Using PLINK2](#removing-population-outliers-using-plink2)

## <ins>Step 1:</ins>
## VCF LiftOver Pipeline

This pipeline provides a step-by-step process to lift over VCF files from `hg19` to `hg38` using GATK's LiftoverVcf tool.

### Data Source

The data was sourced from the [ROMSAP harmonization whole genome sequencing as VCF files](https://www.synapse.org/#!Synapse:syn11707420). To download these files:

1. Navigate to the provided link.

2. Click on "Download Options" to find detailed instructions.

Note: Accessing this data requires a signed Data Use Agreement.

### Prerequisites

1. **Modules**: Before using the script, ensure the following modules are installed and available on your system:
    - `gatk/4.0.10.0`
    - `bcftools`

2. **Files & Directories**:
    - Source VCF files situated in the designated paths.
    - Reference fasta files for `hg19` and `hg38` must be correctly pointed in the script.
    - Chain file for liftOver, moving from `hg19` to `hg38`.

### Usage

1. Clone this repository to your local machine:

```
git clone [repository-url]
cd [repository-dir]
```

### Run the script:

```
./liftover_script.sh
```

### Workflow

- **Concatenation**: Initially, chromosomes from the ROMSAP VCFs are concatenated. The code prepends "chr" in front of each chromosome to match the annotations in the LiftOver chain file.
- **Normalization**: The VCF undergoes normalization using the `bcftools norm` function. This step uses the `hg19 GATK resource bundle` for accuracy, retains multiallelic sites together (avoiding issues later in the pipeline), and swaps ref/alt alleles if incorrectly placed.
- **LiftOver**: This step employs GATK/4.0.10.0's `LiftoverVcf` tool to transform the VCF from `hg19` to `hg38`. The tool uses the `hg19tohg38` chain and references the `hg38` fasta genome file located in the SCG reference folder. It's critical to adjust the VCF header to include "chr" before chromosome definitions to ensure compatibility. The liftOver operation occasionally faces "out of memory" issues. To mitigate this, set `MAX_RECORDS_IN_RAM` to 500. If the problem persists, rerunning the script typically resolves it.
- **Final Normalization**: A final round of normalization is performed post LiftOver.

### Parameters

Inside the script, specific parameters are set at the beginning. These include paths to the required VCF files, reference genomes, and the chain file for LiftOver. To ensure the script works correctly in your environment, adjust these parameters as necessary.

### Additional Details

- The main code for this pipeline is found in `01_liftover_vcf.sh`.
- For reference sequences, the pipeline utilizes `/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa`.
- If you encounter memory issues, rerunning the script is usually a quick solution, though it's not guaranteed to always work.

## <ins>Step 2:</ins>
## PLINK2 VCF Processing Script Instructions

This script provides utilities for converting VCF files into PLINK formats (both PGEN and BED) and extracting a specific sample from a VCF using PLINK2.

### Prerequisites:
1. **PLINK2**: Ensure you have PLINK2 installed. If you're using a system with module environments (like many HPC clusters), the script assumes you have a `plink2` module available.

2. **VCF File**: Have your VCF file ready. This script is designed for VCF files that have been normalized and can handle multiple chromosomes.

### Step-by-step Instructions:

1. **Clone the Repository & Navigate to the Script**:

```
git clone <repository_url>
cd <repository_directory>
```

2. **Script Parameters**: Open the script using your preferred text editor (e.g., `nano`, `vim`, `gedit`).

```
nano plink2_processing.sh
```

Adjust the following parameters at the top of the script:

- `INPUT_VCF_PATH`: Replace `../path_to_your_data/input_vcf_file.vcf` with the path to your input VCF file.

- `OUTPUT_PREFIX`: Replace `../path_to_output_dir/output_prefix` with the desired output directory path and prefix for your output files.

- `SAMPLES_KEEP_PATH`: Replace `../path_to_samples_dir/samples_to_keep.tsv` with the path to your samples file. This file should list the samples you wish to keep, one per line.

- `ONE_SAMPLE_PATH`: Replace `../path_to_samples_dir/one_sample.tsv` with the path to a file containing a single sample name you wish to extract from the VCF.

- Filters (`MIND_FILTER`, `HWE_FILTER`, `GENO_FILTER`, `MAF_FILTER`): Adjust these values as needed. They correspond to PLINK's `--mind`, `--hwe

## <ins>Step 3:</ins>
## Running Mishap Test with PLINK

### Overview
This script runs the mishap test on genomic data using PLINK. The mishap test is designed to detect Mendelian inconsistencies—errors or inconsistencies in the inheritance patterns among related individuals. It's especially useful in studies that involve family data. The script specifically focuses on the dataset 'rosmap_wgs_hg38_all_chroms' and serves as Step 3 in our data analysis pipeline, tailored for research reproducibility.

### Pre-requisites
PLINK must be installed and accessible either through a module system or installed globally.
The dataset 'rosmap_wgs_hg38_all_chroms' should be present at the specified path.

### Usage
To execute the script, run the following command in your terminal, replacing <chromosome_number> with the number of the chromosome you wish to analyze:

```
./03_mishap_test.sh <chromosome_number>
```

### Script Behavior and Output
- The script will exit immediately if any command returns a non-zero status to ensure that all steps are successfully completed.
- If the chromosome number is not specified as an argument, the script will display usage instructions and exit.
- The output of the mishap test will be stored in the directory specified by OUTPUT_DIR.

## <ins>Step 4:</ins>
## Pre-processing for ROSMAP QTL Analysis

### Overview

This script accomplishes various tasks to prepare the ROSMAP QTL dataset for further analysis. The tasks include:

1. Processing missing haplotype test results generated in Step 3.
2. Identifying significant variants based on a predefined xQTL study threshold.
3. Preparing lists of files to be merged and filtered for data alignment.
4. Matching and filtering samples to align with methylation data.

### Pre-requisites

- Make sure the R libraries `data.table` and `tidyverse` are installed.
- Ensure that you have access to various data files related to ROSMAP QTL, as specified in the script.
  
### Dependencies

- **R Libraries**: `data.table`, `tidyverse`
  
### Input Files

- Missing haplotype results file: `all_chroms_hg38.missing.hap`
- Metadata for whole genome sequencing: `ROSMAP_assay_wholeGenomeSeq_metadata.csv`
- Matched clinical metadata: `matched_clincal_meta.rds`

### Output Files

- List of significant variants to exclude: `exclude_variants_hg38.txt`
- List of BIM files for merging: `bim_files_hg38.txt`
- List of refiltered BIM files: `bim_files_hg38_refiltered.txt`
- List of samples to keep: `keep_samples.txt`
- List of samples to remove: `remove_samples.txt`

### Usage

Execute the R script in your R environment.

```
Rscript 04_filter_mishap.R
```


### Important Variables

- `thresh`: Threshold for significance in xQTL study. Default is \(1 \times 10^{-9}\).
- `savedir`: Output directory where the processed datasets are saved.


## <ins>Step 5:</ins>
## Running PLINK and PLINK2 for data formatting

### Overview
This script automates several steps in managing genomic data, particularly by using PLINK and PLINK2 to filter, merge, and export genomic data files. It is designed to run on the SCG computational environment and has several key functionalities:

1. Identifies multiallelic sites
2. Removes multiallelic sites from individual chromosome files
3. Merges cleaned chromosome files
4. Creates counts of reference alleles
5. Exports a PED file for the merged dataset

### Pre-requisites
Make sure PLINK and PLINK2 are loaded in your computational environment. The script automatically does this by running module load plink and module load plink2.

### Usage
Navigate to the directory where the script is stored and execute the script by running:

```
./05_merge_bims.sh
```

### Script Contents
The script executes the following steps:

- Step 1: Identifies multiallelic sites and saves them to a BED file.
- Step 2: Iterates through chromosomes 1 to 22, removing the identified multiallelic sites.
- Step 3: Merges the cleaned individual chromosome files.
- Step 4: Counts reference alleles and saves the counts to a file.
- Step 5: Creates a PED file for the merged dataset.

For detailed script logic, refer to the comments in the script itself.

## <ins>Step 6:</ins>
## Principal Component Analysis (PCA) using Eigensoft SmartPCA

### Description

This step utilizes the SmartPCA tool from the Eigensoft package to perform a Principal Component Analysis (PCA) on genotypic data. This is a Perl script that acts as a wrapper around the SmartPCA tool, tailored specifically for the ROSMAP QTL analysis.

### Prerequisites

- **Perl**: Make sure you have Perl installed on your system.
- **Eigensoft SmartPCA**: Ensure that you've installed Eigensoft and that the `smartpca.perl` executable is in your system's PATH.

### Script Details

The script is located in the repository as `06_eigenstrt_pca.pl`.

#### Parameters

- Input geno file in PACKEDPED format
- Input snp file in PACKEDPED format
- Input indiv file
- Number of Principal Components (PCs) to output (default = 10)
- Number of maximum iterations (default = 5)
- Number of PCs to remove outliers during each iteration (default = 10)
- Sigma, the number of standard deviations for outlier removal (default = 6.0)

#### Output Files

- PCA results
- Plot file
- Eigenvalues file
- Log file

### Usage

To run the script, navigate to the directory where the script is located and execute:

```bash
perl 06_eigenstrt_pca.pl
```

### Troubleshooting
If you encounter an error message saying that `smartpca.perl` is not found, make sure that the Eigensoft SmartPCA executable is correctly placed in your system's PATH.

### Note
Before running the script, you might need to adjust the hardcoded directories and file paths according to your specific setup.


## <ins>Step 7:</ins>
## Formatting and Analyzing Eigenstrat Output

### Overview

This R script is designed to read, format, and analyze Eigenstrat output files generated from a QTL (Quantitative Trait Loci) analysis. Specifically, the script performs the following tasks:

1. Reads and processes Eigenvalues from PCA (Principal Component Analysis) to calculate the variance explained by each PC.
2. Reads and formats Principal Components (PCs) for each individual in the dataset.
3. Identifies and handles samples that have been removed in previous steps.
4. Correlates PCs with various traits to identify significant associations.

### Prerequisites

- R (Version 3.6.0 or above)
- ggplot2 library
- tidyverse library

You can install the required R packages with the following command:

```R
install.packages(c("ggplot2", "tidyverse"))
```

### Directory Structure

Before running the script, make sure you have the following directory structure:

- `output/`
  - `05_qtl_analysis/`
    - `07_geno_pca/`
      - `eigenstrat_smartpca_hg38.pca.evec`
    - `05_filter_mishap/`
      - `keep_samples.txt`
  - `01_preprocessing/`
    - `01_subset_samples/`
      - `matched_clincal_meta.rds`

### How to Run

Navigate to the directory where the script is saved and execute the R script by running the following command in your R console:

```R
source("07_format_eigenstrat.R")
```

### Outputs

The script generates the following output files:

1. `pca_pct_var_hg38.tsv` - Contains the percentage of variance explained by each PC.
2. `genotype_pcs_hg38.rds` - Contains formatted PCs for each individual.
3. `population_outliers_hg38.txt` - List of samples that were removed.
4. `genotype_pc_trait_lm_pvals_hg38.rds` - Contains p-values of correlations between PCs and traits.
5. `genotype_pc_trait_pvals_heatmap.png` - Heatmap of the p-values.

The output files will be saved in the `08_format_eigenstrat/` directory under `output/05_qtl_analysis/`.



## <ins>Step 8:</ins>
## Removing Population Outliers Using PLINK2

### Overview

This section describes a Bash script that uses PLINK2 to remove population outliers from the dataset. The script creates new PLINK binary bed files without the outliers, making it easier to proceed with downstream genetic analyses.

### Requirements

- PLINK2: Make sure that the PLINK2 module is available on your system. 

### Configuration

The script uses the following directory variables, which you can modify as per your project's directory structure:

- `INPUT_DIR`: The directory containing the input .pfile (default is `../../output/05_qtl_analysis/06_merge_bims/`)
- `OUTPUT_DIR`: The directory where the output files will be saved (default is `../../output/05_qtl_analysis/09_remove_outliers/`)
- `OUTLIER_FILE`: The file containing the list of population outliers to be removed (default is `../../output/05_qtl_analysis/08_format_eigenstrat/population_outliers_hg38.txt`)

### Usage

To execute the script, navigate to its location and run:

```bash
./08_remove_outliers.sh
```

### Script Options

- `--make-bed`: Creates PLINK binary bed files.
- `--output-chr chrM`: Formats mitochondrial chromosome names as 'chrM'.
- `--pfile`: Specifies the input .pfile location.
- `--remove`: Specifies the list of samples to be removed.
- `--out`: Specifies the location for the output files.

### Output

The script will generate new PLINK binary bed files without the outliers and save them in the directory specified by the `OUTPUT_DIR` variable.

