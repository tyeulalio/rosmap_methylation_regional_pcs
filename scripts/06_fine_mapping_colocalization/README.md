# Fine-mapping / Colocalization

This pipeline contains a set of functions designed to perform statistical fine-mapping after QTL analysis and colocalization.

## Table of Contents

- [01_check_data](#dir-01)
    - [01_check_data.R](#script-01)
- [02_gwas_prep](#dir-02)
    - [01_liftover_gwas.sh](#script-02)
    - [02_check_liftover.R](#script-03)
    - [03_create_gwas_vcfs_for_ld_matrix.sh](#script-04)
    - [04_finemap_gwas.sh](#script-05)
    - [05_check_results.R](#script-06)
    - [06_format_finemapped_gwas.sh](#script-07)
- [03_qtl_prep](#dir-03)
    - [01_create_torus_input.sh](#script-08)
    - [02_run_torus.sh](#script-09)
    - [03.0_create_sbams_genotype_files.sh](#script-10)
    - [03.1_create_qtl_dapg_input_parallel.py](#script-11)
    - [04_run_qtl_dapg.sh](#script-12)
    - [05_make_fastqtl_qtl_annotations.sh](#script-13)
    - [06_format_fastqtl_qtl_input.py](#script-14)
- [04_fastenloc](#dir-04)
    - [01_run_fastenloc.sh](#script-15)
    
## <a name="dir-01">`01_check_data/`</a>

## Script 01                                    
## 01_check_data.R

This R script is designed to validate the alignment of genomic annotations between Genome-Wide Association Study (GWAS) data and Quantitative Trait Loci (QTL) data. The script's focus is on ensuring consistency between different genomic builds and preparing the data for further analysis or liftover.

## Prerequisites

Before running this script, ensure you have the following libraries installed:
- `GenomicRanges`: Provides data structures to represent and manipulate genomic intervals.
- `tidyverse`: A collection of R packages designed for data science, including data manipulation and visualization tools.

## Functionality

This script performs the following tasks:

1. **Load GWAS Data**:
   - Reads a GWAS dataset from a compressed text file. The data contains key information on genetic variants, including chromosome, position, and allele information.
   
2. **Format for Liftover**:
   - Transforms the GWAS data into a format suitable for liftover, including appropriate columns for chromosome, position, reference allele, and alternative allele.
   - Ensures consistent naming conventions and cleans the data by removing invalid or non-standard alleles.
   - Changes the `ID` field to a unique identifier combining chromosome, position, and alleles.

3. **Save as VCF**:
   - Writes the formatted GWAS data to a VCF (Variant Call Format) file, a widely used format for storing genomic variants.
   - Adds the appropriate VCF header and column names to ensure compatibility with downstream tools.

## Execution Steps

1. **Load the GWAS Data**:
   - Load the dataset from the specified file path, extracting the key columns required for liftover and further analysis.
   
2. **Format the Data**:
   - Transform the GWAS data to ensure compatibility with liftover and other genomic operations.
   - Create a unique `ID` field and adjust column names as necessary.
   
3. **Save the Formatted Data**:
   - Write the formatted GWAS data to a VCF file, ensuring it is in a standard format suitable for subsequent analysis.

## Example Usage

```r
# Load GWAS data
datafile <- "/path/to/data/ad_gwas/wightman_etal_2021_results.txt.gz"
gwas <- read_tsv(datafile)

# Format and save as VCF
formatted_gwas <- format_gwas(gwas)
save_vcf(formatted_gwas, "/path/to/output/formatted_gwas.vcf")
```

Ensure the paths to input and output files are correct before running the script. This script is designed to be modular, allowing for integration into larger genomic pipelines or standalone use to check and format GWAS data.

---

## <a name="dir-02">`02_gwas_prep/`</a>

## Script 02                                    
## 01_liftover_gwas.sh

This bash script performs a liftover operation to convert variant call format (VCF) files from one genomic build to another, specifically from the hg37 to hg38 build. The script uses GATK's LiftoverVcf utility to remap genomic coordinates based on a predefined chain file, ensuring that variant annotations are consistent with the reference genome.

## Prerequisites

Before executing this script, ensure the following:

- The GATK module is loaded and available in your environment. This script specifies GATK version 4.0.10.0, but you can adjust as needed.
- A liftover chain file that defines the mapping between the hg19 and hg38 genomic builds.
- A reference genome in FASTA format corresponding to hg38.

## Script Overview

The script performs the following tasks:

1. **Load GATK**:
   - Ensures that the GATK module is loaded to access the LiftoverVcf utility.

2. **Set Input and Output Files**:
   - Defines the input VCF file to be lifted over.
   - Specifies the output VCF file after liftover and the file for rejected records.

3. **Perform Liftover**:
   - Executes the GATK LiftoverVcf command to remap the coordinates in the VCF file from hg37 to hg38.
   - Uses the provided chain file and reference sequence for accurate conversion.
   - Captures rejected records for further investigation.

## Execution Steps

1. **Load GATK**:
   - Ensure GATK and related modules are loaded into your environment.

2. **Set Data Paths**:
   - Define the input data file, output file, and path to the liftover chain.
   - Ensure paths are accurate and point to existing files.

3. **Run Liftover**:
   - Execute the GATK LiftoverVcf command with the specified parameters.

## Example Usage

```bash
# Activate the GATK module
module load gatk/4.0.10.0

# Perform the liftover operation
gatk LiftoverVcf \
    --INPUT="/path/to/ad_gwas/wightman_etal_2021_results.vcf" \
    --OUTPUT="/path/to/output/01_liftover_gwas/wightman_etal_2021_results_hg38.vcf" \
    --VERBOSITY=DEBUG \
    --CHAIN="/path/to/liftOver/chains/hg19ToHg38.over.chain.gz" \
    --REFERENCE_SEQUENCE="/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa" \
    --REJECT="/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/01_liftover_gwas/rejected_records.vcf" \
    --MAX_RECORDS_IN_RAM 10000 \
    --RECOVER_SWAPPED_REF_ALT
```

## Additional Notes

- The `VERBOSITY=DEBUG` option is used for detailed output during liftover. Adjust to `ERROR` or `INFO` for less detailed output.
- The `--RECOVER_SWAPPED_REF_ALT` flag is set to recover variants with swapped reference and alternative alleles.
- Make sure to have sufficient memory and compute resources, as liftover operations can be resource-intensive.
  
Check for errors or rejections after running the script to ensure the liftover was successful.

---

## Script 03                                    
## 02_check_liftover.R

This R script is designed to integrate and check the consistency of GWAS (Genome-Wide Association Studies) and QTL (Quantitative Trait Loci) SNP data. It includes functions to examine liftover consistency, join GWAS and QTL data, check allele swapping, and create LD (Linkage Disequilibrium) annotations. The script also creates Z-scores for GWAS data and performs various quality checks.

## Prerequisites

Ensure that you have installed the following R packages:

- `GenomicRanges`
- `liftOver`
- `ggplot2`
- `tidyverse`
- `data.table`

Additionally, ensure you have access to the appropriate GWAS data files and chain files for genomic liftover operations.

## Script Structure

This script consists of several key functions, each responsible for a specific task:

- **check_failed()**: Checks variants that failed to map during liftover and attempts to correct them.
- **load_qtl()**: Loads QTL SNP data from specified files.
- **join_gwas_qtls()**: Joins GWAS and QTL data, checks for allele swapping, and identifies mismatches.
- **load_gwas_scores()**: Loads GWAS scores, computes Z-scores, and addresses potential infinite values.
- **check_zscores()**: Conducts quality checks on computed Z-scores to ensure consistency.
- **create_ld_annotations()**: Creates LD matrix annotations for SNPs based on predefined LD blocks.
- **load_gwas_snps()**: Loads GWAS SNP data with position and allele information.
- **main()**: Orchestrates the script's flow, integrating the above functions to perform GWAS-QTL matching and create annotations.

## Execution Flow

1. **Load GWAS Data**:
   - Loads the specified GWAS dataset with genomic coordinates and allele information.

2. **Load QTL Data**:
   - Loads the QTL SNPs from a specified file for later comparison with GWAS data.

3. **Match GWAS and QTL**:
   - Joins GWAS and QTL data on chromosome and position.
   - Identifies and addresses variants with swapped alleles.
   - Removes mismatched variants.
   - Saves the final matched variants to a specified output file.

4. **Create LD Annotations**:
   - Uses LD blocks to annotate the matched GWAS and QTL SNPs.
   - Generates the corresponding annotations and saves them to a specified output directory.

5. **Generate Z-Scores**:
   - Computes Z-scores for the GWAS data and saves the updated scores.
   - Ensures proper handling of infinite values and selects unique SNPs based on LD blocks.

## Additional Notes

- Ensure that paths to input data, output directories, and reference files are correctly set.
- Regularly check for missing values, swapped alleles, or other inconsistencies during data processing.
- This script is flexible and allows for customization based on different GWAS and QTL datasets, so make sure to adjust parameters and paths accordingly.

---

## Script 04                                    
## 03_create_gwas_vcfs_for_ld_matrix.sh

This bash script creates Variant Call Format (VCF) files for the ROSMAP dataset using Plink2. It generates files suitable for constructing LD (Linkage Disequilibrium) matrices for further analyses, such as fine-mapping and colocalization. The script requires a conda environment with the necessary Plink and other genomic tools activated.

## Prerequisites

Ensure that you have installed the following software and packages:

- `Plink2`
- `Conda/Mamba`
- ROSMAP WGS (Whole Genome Sequencing) data in VCF format.

Additionally, verify that the specified conda environment is properly configured and activated with Plink2 installed.

## Execution Steps

1. **Environment Activation**:
   - Activate the conda environment where Plink2 is installed.
   - Ensure that all necessary dependencies are available.

2. **Output Directories**:
   - Create output directories for storing the generated files, specifically for Plink and LD matrix files.
   - Ensure that these directories exist before running the script.

3. **Plink Operations**:
   - **Create BED Files**: Converts the input VCF file to Plink's BED format, allowing for further manipulation. This operation includes:
     - Setting chromosome inclusion (`--chr 1-22`).
     - Allowing extra chromosome formats (`--allow-extra-chr`).
     - Handling missingness (`--mind 0.05`, `--geno 0.05`).
     - Limiting alleles to two (`--max-alleles 2`).
   - **Handle Duplicate SNP IDs**:
     - Checks for SNPs with duplicate IDs and creates a VCF containing only these SNPs.
     - Removes duplicates from the main BED file.

4. **Quality Checks**:
   - Ensure all Plink operations complete without errors by reviewing standard output and error logs.
   - Validate that duplicate SNPs are handled appropriately.

## Output Files

- **Plink2 BED Files**: Stores the primary genomic data for ROSMAP in Plink's binary format.
- **VCF Files**: Stores VCF files with corrected SNP IDs and without duplicates.
- **Plink Logs**: Contains stdout and stderr logs from Plink operations for debugging and quality checks.

## Important Notes

- This script should be executed only once, as it creates foundational files for subsequent analyses.
- Ensure that paths to data files, output directories, and other dependencies are correctly set.
- Adjust the script parameters according to your dataset's structure and desired filtering criteria.

---
    
## Script 05                                    
## 04_finemap_gwas.sh

This script fine-maps Genome-Wide Association Studies (GWAS) using DAP-G, a Bayesian approach to identify credible sets of causal variants. It integrates with a [QTL pipeline](
https://github.com/tyeulalio/QTL_pipeline) and requires specific input files and conda environments.

## Prerequisites

Ensure the following software and packages are installed and activated:

- `Conda/Mamba`
- `Python3`
- `DAP-G`
- `Plink2`

Additionally, verify that the specified conda environment is correctly set up with all dependencies installed.

## Execution Parameters

- `TMPDIR`: Temporary directory for intermediate files during fine-mapping.
- `REGIONNUM`: The region number to fine-map.

## Execution Steps

1. **Environment Activation**:
   - Activate the conda environment with DAP-G and other necessary packages.
   - Ensure the proper environment activation scripts are sourced.

2. **Output Directories**:
   - Create the output directory for storing fine-mapping results.
   - Ensure the directory exists before executing the script.

3. **Input Files**:
   - **SNP Map File**: Contains the mapping from SNP positions to variant names used in QTL analysis.
   - **Z-Score File**: Contains GWAS summary statistics, including Z-scores for fine-mapping.
   - **Plink File**: Plink-formatted genotype data for the GWAS analysis.
   - Verify that these files are accessible and correctly formatted.

4. **Fine-Mapping with DAP-G**:
   - The script uses a Python-based pipeline to fine-map GWAS using DAP-G.
   - Key parameters include:
     - `-m`: Path to the SNP map file.
     - `-r`: Region number for fine-mapping.
     - `-z`: Path to the Z-score file.
     - `-o`: Output directory.
     - `-p`: Path to the Plink genotype data.
     - `-th`: Number of threads for parallel processing.
     - `-t`: Temporary directory for intermediate files.
   - Adjust parameters based on your environment and dataset specifics.

5. **Script Execution**:
   - Execute the script with the correct parameters.
   - Monitor stdout/stderr for errors and ensure the expected output is generated.

## Output Files

- **Fine-Mapped GWAS Data**: Contains the results from DAP-G fine-mapping, with credible sets of causal variants.
- **Logs**: Standard output and error logs from the fine-mapping process for debugging and quality checks.

## Important Notes

- This script requires significant computational resources, depending on the dataset size and fine-mapping complexity.
- Ensure that paths to all input files, output directories, and other dependencies are correctly set before executing the script.
- Adjust the script parameters according to your dataset structure and desired fine-mapping resolution.

---

## Script 06                                    
## 05_check_results.R

This script is designed to check for missing fine-mapped files in the DAP-G output for a specific Genome-Wide Association Study (GWAS) type.

## Prerequisites

Ensure the following libraries are installed and loaded:

- `tidyverse`

## Execution Steps

1. **Define Parameters**:
   - `gwas_type`: Specifies the type of GWAS analysis to check.
   - `datadir`: Directory path where fine-mapped files from DAP-G are stored.

2. **Retrieve Existing Files**:
   - List all files in the specified directory and count the total number.
   - Display a sample of the existing files to ensure correct retrieval.

3. **Identify LD Blocks**:
   - Extract the LD (Linkage Disequilibrium) block identifiers from the file names.
   - Remove additional text and keep only the LD block numbers for comparison.

4. **Identify Missing Files**:
   - Create a range of expected LD block numbers.
   - Use `setdiff` to identify which LD blocks are missing from the expected range.

## Output

The script returns the LD block numbers that are missing from the expected range. This information helps identify which fine-mapped files might need to be re-generated or investigated further.

## Important Notes

- The expected range of LD blocks might vary based on the dataset or analysis type. Adjust the range accordingly.
- Ensure the directory paths are correctly set to retrieve the fine-mapped files.
- The script assumes a specific naming convention for the fine-mapped files. Modify as needed if the file structure differs.

---

## Script 07                                    
## 06_format_finemapped_gwas.sh

This script formats the fine-mapped Genome-Wide Association Study (GWAS) output to be compatible with FastEnloc for colocalization analysis. It integrates with the [QTL pipeline](
https://github.com/tyeulalio/QTL_pipeline).

## Prerequisites

Ensure the following environment and tools are set up:

- Micromamba or Conda is installed and activated.
- `FastEnloc` and its dependencies are available in the conda environment.
- The output directory structure is prepared.

## Execution Steps

1. **Activate Conda Environment**:
   - Ensure the correct environment is activated for the colocalization analysis.

2. **Define Parameters**:
   - `GWASTYPE`: Specifies the type of GWAS analysis for which the fine-mapped output will be formatted.
   - `OUTPUT_DIR`: Directory path for storing the formatted fine-mapped output.
   - `DAP_DIR`: Directory path where the fine-mapped GWAS results from DAP-G are stored.

3. **Create Output Directory**:
   - Ensure the specified `OUTPUT_DIR` exists or create it if not.

4. **Run Formatting Script**:
   - Use a Python script to format the DAP-G fine-mapped GWAS output for FastEnloc.
   - Provide the output directory and fine-mapped DAP-G directory as parameters to the script.

## Notes

- Ensure that the `DAP_DIR` contains the expected fine-mapped files from the previous DAP-G analysis.
- The `OUTPUT_DIR` should have sufficient space to store the formatted outputs.
- The script requires appropriate permissions to read from `DAP_DIR` and write to `OUTPUT_DIR`. Adjust permissions as needed.

---

## <a name="dir-03">`03_qtl_prep/`</a>

## Script 08                                    
## 01_create_torus_input.sh

This script generates input data for Torus from the output of TensorQTL analysis, preparing files for further statistical analysis and colocalization. It integrates with the [QTL pipeline](
https://github.com/tyeulalio/QTL_pipeline).

## Prerequisites

- Ensure Python 3 is installed and the script `01_create_torus_input.py` is accessible.
- Confirm TensorQTL output files are available for processing.
- Set the necessary command line arguments.

## Command Line Arguments

- `SUMMARY_TYPE`: Specifies the type of summary (e.g., averages, rPCs, CpGs).
- `CELL_TYPE`: Indicates the cell type used in the analysis.
- `REGION_TYPE`: Specifies the genomic region type.
- `PROPORTION_TYPE`: Refers to whether cell-type proportions were included in the analysis.

## Output Directory

- `OUTPUT_DIR`: The directory where Torus input files will be stored. This is set based on the provided `PROPORTION_TYPE`.

## Script Execution

1. **Set Up Output Directory**:
   - Ensure the specified output directory exists or create it if it doesn't.

2. **Define File Paths**:
   - `QTL_PREFIX`: A prefix for identifying the correct TensorQTL output files based on the given summary type, cell type, and region type.
   - `QTL_FILES`: A pattern that matches all relevant TensorQTL output files.

3. **Execute the Script**:
   - Run the Python script to create the Torus input, providing the necessary command line arguments for output directory, run name, and TensorQTL files.

## Example Command

```bash
./01_create_torus_input.sh avgs neurons full_gene with_proportions
```

Ensure to replace `HOME` and other placeholders with the actual paths on your system.

## Additional Notes

- This script creates Torus input files using the TensorQTL cis-QTL results, stored in Parquet format.
- Ensure appropriate permissions for reading and writing to the specified directories.
- Verify that the `QTL_FILES` pattern correctly matches the expected TensorQTL output files.

---

## Script 09                                    
## 02_run_torus.sh

# Running Torus to Obtain SNP PIPs

This script runs Torus to calculate Posterior Inclusion Probabilities (PIPs) for Single Nucleotide Polymorphisms (SNPs), providing a measure of the likelihood that each SNP is causally associated with the trait of interest. It integrates with the [QTL pipeline](
https://github.com/tyeulalio/QTL_pipeline).

## Command Line Arguments

- `SUMMARY_TYPE`: Type of data summarization used (e.g., averages, rPCs, CpGs).
- `CELL_TYPE`: Cell type being analyzed.
- `REGION_TYPE`: Genomic region type.
- `PROPORTION_TYPE`: Indicates whether cell-type proportions were included in the analysis.

## Output Directories

- `OUTPUT_DIR`: The primary output directory for storing Torus results, organized by proportion type and run name.
- `RUNNAME`: A unique identifier for the analysis run, derived from the summary type, cell type, and region type.

## Script Execution

1. **Set Up Output Directories**:
   - Create the main output directory if it doesn't exist.
   - Create sub-directories based on the `PROPORTION_TYPE` and `RUNNAME`.

2. **Define Input File**:
   - `INPUT_FILE`: Path to the file containing the data for Torus analysis, derived from the previous script that creates Torus input.

3. **Run Torus**:
   - Call the Torus script to calculate SNP PIPs, passing in the input file and output directory.

## Example Command

```bash
./02_run_torus.sh avgs neurons full_gene with_proportions
```

## Additional Notes

- Ensure the necessary Python environment is activated before running Torus.
- Confirm that the `INPUT_FILE` is correctly generated by the previous script and is accessible.
- Adjust permissions as needed to allow reading and writing to the specified directories.

--- 

## Script 10                                    
## 03.0_create_sbams_genotype_files.sh

This script creates a template file for the genotype counts for each genotype from an existing set of sbam files. The templates are used by a subsequent script to process genotype data.

## Directory Paths

- `SBAM_DIR`: Path to the directory containing sbam files from previous processing steps.
- `OUTPUT_DIR`: Directory to store the output genotype count templates. This is created if it doesn't exist.

## Script Execution

1. **Create Output Directory**:
   - If the output directory doesn't exist, create it to store the resulting genotype count templates.

2. **Process Sbam Files**:
   - Iterate over each file in the `SBAM_DIR`.
   - Extract the gene name from the sbam file name.
   - Create an output file for each gene's genotype counts.
   - Extract lines starting with `geno` from the sbam file and save them to the corresponding output file.

3. **Progress Monitoring**:
   - The script outputs progress every 1000 genes processed, indicating the current count of processed genes.

---

## Script 11                                   
## 03.1_create_qtl_dapg_input_parallel.py

This script generates input files for ENLOC analysis by processing genotype and phenotype data and creating prior probability files. The main() function orchestrates these tasks by extracting gene information, loading genotype data, and processing phenotype and variable data.

## Command Line Arguments

- `summary_type`: Index representing the summary type (e.g., averages, principal components, or CpGs).
- `cell_type`: Index for the cell type (e.g., bulk, astrocytes, endothelium, neurons, oligodendrocytes).
- `region_type`: Index for the region type (e.g., full gene, promoters, gene body, preTSS).
- `job_idx`: Index for the current job, used for dividing work across multiple jobs.
- `savedir`: Directory where output files will be saved.
- `num_jobs`: Total number of jobs to be processed.
- `priors_dir`: Directory containing prior probability files for ENLOC.
- `pheno_dir`: Directory containing phenotype data files.
- `genotype_dir`: Directory with genotype data files.

## Script Tasks

1. **Directory Setup**:
   - Create output directories for saving results, priors, and sbam files.

2. **Gene List Extraction**:
   - Retrieve a list of genes to process from the phenotype file.
   - Save the list of genes for reference.

3. **Process Phenotype Data**:
   - Read the phenotype data and format it into an appropriate structure for sbam files.
   - Check that the order of samples matches between phenotype and genotype data.

4. **Process Genotype Data**:
   - Use the gene name to identify the corresponding genotype counts file.
   - Format and save the genotype data to the sbam file.
   - Create and save prior probability files for ENLOC.

5. **Process Controlled Variables**:
   - Extract and format controlled variable data from the matched covariates file.
   - Save the controlled variable data to the sbam file.

6. **Progress Reporting**:
   - Output the number of genes processed and monitor progress during script execution.

## Execution and Usage

To run the script, execute it with the required command line arguments:

```bash
python3 /path/to/script.py <summary_type> <cell_type_number> <region_type_number> <job_idx> <savedir> <num_jobs> <priors_dir> <pheno_dir> <genotype_dir>
```

Ensure all necessary dependencies are installed and accessible, and appropriate permissions are set to read from and write to the specified directories.

---

## Script 12                                  
## 04_run_qtl_dapg.sh

This script runs the Deterministic Approximation of Posteriors (DAP-G) to obtain probabilistic QTL annotations for fastENLOC, following instructions from GTEx's DAP-G analysis.

## Command Line Arguments

- `i`: Gene index indicating which gene to process.
- `SUMMARY_TYPE`: Type of data summary (e.g., averages, principal components, or CpGs).
- `CELL_TYPE`: Cell type (e.g., bulk, astrocytes, neurons, etc.).
- `REGION_TYPE`: Genomic region type (e.g., full gene, promoters, gene body, preTSS).

## Execution Steps

1. **Setting Up Environment**:
   - Activate the necessary conda environment with colocalization tools.
   - Set up output directories for storing DAP-G results.

2. **Creating File Names and Paths**:
   - Construct the run string based on summary type, cell type, and region type.
   - Define the path to the filtered genes list and generate it if needed.
   - Define the output directory structure for DAP-G results.

3. **Running DAP-G**:
   - Retrieve the gene to process based on the gene index.
   - Format the gene name to ensure compatibility with DAP-G.
   - Check if the output file already exists to avoid redundant processing.
   - Execute DAP-G with the specified sbam and prior files, using parallel threads and an LD control threshold of 0.25.

4. **DAP-G Execution Command**:
   - The following command is used to run DAP-G:

     ```bash
     dap-g -d "${GENESDIR}/sbams/${GENE}_fastqtl_singl_snp_output.sbam" \
          -p "${GENESDIR}/priors/${GENE}_priors.txt" \
          -ld_control 0.25 \
          --all \
          -t 4 \
          -n ${GENE} \
          > $OUTFILE
     ```

## Running the Script

To run the script, execute it with the gene index and other required arguments:

```bash
bash script_name.sh <gene_index> <summary_type> <cell_type> <region_type>
```

Ensure that the output directories are properly set and the required files are available in their expected locations. If the output file already exists, the script will exit to avoid overwriting.

---

## Script 13                                    
## 05_make_fastqtl_qtl_annotations.sh

This script processes fine-mapped QTL data to create SNP annotations in VCF format for FastEnloc using the `summarize_dap2enloc.pl` script from FastEnloc's GitHub tutorial . It requires a directory with fine-mapped QTL results and a VCF file to annotate all SNP positions, which is created in previous scripts.

## Command Line Arguments

- `SUMMARYTYPE`: Type of summary data (e.g., averages, principal components, or CpGs).
- `CELLTYPE`: Type of cell (e.g., bulk, astro, neurons, etc.).
- `REGIONTYPE`: Genomic region type (e.g., full gene, gene body, promoters, preTSS).

## Execution Steps

1. **Set Up Environment**:
   - Create an output directory for the annotation results.
   - Define the DAP-G directory containing the fine-mapped QTL results.
   - Define the path to the VCF file for SNP annotation.

2. **Run FastQTL**:
   - Use `summarize_dap2enloc.pl` to create annotations from the fine-mapped QTL results.
   - The SNP annotation is output as a compressed VCF file for FastEnloc.

3. **FastQTL Command**:
   - The following command is used to create annotations for FastEnloc:

     ```bash
     /path/to/programs/fastenloc/src/summarize_dap2enloc.pl \
         -dir $DAPDIR \
         -vcf "/path/to/data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms.vcf.gz" \
         | gzip - > "${OUTDIR}/fastenloc.qtl.annotation.vcf.gz"
     ```

## Running the Script

To execute the script, use the following command with the required arguments:

```bash
bash script_name.sh <summary_type> <cell_type> <region_type>
```

Ensure that the DAP-G directory and the VCF file are correctly set. The script generates an output file containing annotations in VCF format for further processing with FastEnloc.

---

## Script 14                                    
## 06_format_fastqtl_qtl_input.py

This script formats the output from QTL DAP-G into a VCF (Variant Call Format) file for further analysis with FastQTL. It reads fine-mapped QTL annotations and generates a VCF file with a gene-based structure.

## Command Line Arguments

- `summary_type`: Type of summary data (e.g., averages, principal components, or CpGs).
- `cell_type`: Type of cell (e.g., bulk, astro, neuron, etc.).
- `region_type`: Genomic region type (e.g., full gene, gene body, promoters, preTSS).

## Execution Steps

1. **Set Up Output Directory**:
   - Define a base output directory for storing the formatted VCF file.
   - Create subdirectories for specific summary, cell, and region types.

2. **Format VCF Output**:
   - Read fine-mapped QTL annotations from the specified file.
   - Format each annotation to include the corresponding gene before the SNP identifier.
   - Save the formatted annotations to the output VCF file.

3. **Compress the Output File**:
   - After formatting the VCF file, compress it with gzip for efficient storage and use with FastQTL.

## Running the Script

To execute the script, use the following command with the required arguments:

```bash
python3 script_name.py <summary_type_index> <cell_type_index> <region_type_index>
```

Ensure that the appropriate index values are used for summary, cell, and region types. The script generates a formatted VCF file and compresses it for FastQTL use.

---

## <a name="dir-04">`04_fastenloc/`</a>

## Script 15                                   
## 01_run_fastenloc.sh

This script runs FastENLOC to perform colocalization analysis on QTL and GWAS data. It requires input VCF files with fine-mapped QTL annotations and GWAS data, and outputs the colocalization results.

## Command Line Arguments

- `SUMMARYTYPE`: The summary type (e.g., averages, principal components, or CpGs).
- `CELLTYPE`: The cell type (e.g., bulk, astro, neuron, etc.).
- `REGIONTYPE`: The genomic region type (e.g., full gene, gene body, promoters, preTSS).
- `GWASTYPE`: The type of GWAS data being used (e.g., Wightman, Kunkle).
- `NUM_GWAS_VARS`: The total number of GWAS variants tested.

## Output Directory Structure

The script creates the following directory structure to store FastENLOC outputs:

1. **Main Output Directory**:
   - A base directory for storing FastENLOC results.
   - The structure is organized by summary, cell, region types, and GWASTYPE.

2. **Subdirectories**:
   - Specific subdirectories for different summary, cell, and region types, and GWASTYPE.

## Execution Steps

1. **Create Output Directory**:
   - Create the main output directory for storing FastENLOC results.
   - Create specific subdirectories based on the summary, cell, and region types, and GWASTYPE.

2. **Run FastENLOC**:
   - Use the `-total_variants` flag to specify the total number of GWAS variants tested.
   - Set the number of threads with the `-thread` flag (here, it uses 1 thread).
   - Define the output prefix for FastENLOC results.
   - Specify the input QTL file and GWAS file for colocalization analysis.

## Running the Script

To execute the script, use the following command with the required arguments:

```bash
bash script_name.sh <SUMMARYTYPE> <CELLTYPE> <REGIONTYPE> <GWASTYPE> <NUM_GWAS_VARS>
```

Ensure the appropriate values are provided for summary, cell, region types, and GWASTYPE, along with the total number of GWAS variants tested. The script performs colocalization analysis using FastENLOC and saves the results in the specified output directory.