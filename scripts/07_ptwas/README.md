# PTWAS

This pipeline contains a set of functions designed to perform PTWAS analysis.

## Table of Contents


- [01_compute_ptwas_weights.sh](#script-01)
- [02_concat_weights.sh](#script-02)
- [03_run_make_gambit.sh](#script-03)
- [04_format_gwas.R](#script-04)
- [05_check_ld_panel.R](#script-05)
- [06_swap_ld_panel.sh](#script-06)
- [07_ptwas_scan.sh](#script-07)
- [08_format_scan_results.R](#script-08)
- [11_run_intact.R](#script-11)


## Script 01                                    
## 01_compute_ptwas_weights.sh
This script creates the required files to build Transcriptome-Wide Association Study (TWAS) weights using the PTWAS Builder tool. It takes input from fine-mapped DAP files and SBAM files and outputs the PTWAS weights for the specified gene.

## Command Line Arguments

- `DAPFILE`: The input DAP-G file containing fine-mapped QTL data.
- `SBAMFILE`: The input SBAM file with gene-level information.
- `GENE`: The name of the gene for which to build PTWAS weights.
- `OUTDIR`: The base output directory to store the PTWAS weights.

## Output Directory Structure

The script creates the following directory structure to store PTWAS weight files:

1. **Weights Output Directory**:
   - A subdirectory within the specified output directory to store the PTWAS weights.
   - The structure is organized by gene names.

2. **Output File**:
   - A text file containing the PTWAS weights for the specified gene.

## Execution Steps

1. **Create Output Directory**:
   - Create the main output directory for storing PTWAS weight files.
   - Create a specific subdirectory for the PTWAS weights based on the gene name.

2. **Run PTWAS Builder**:
   - Use the PTWAS Builder tool to create the weights for the specified gene.
   - Provide the DAP file, SBAM file, gene name, and output file to the PTWAS Builder command.

## Running the Script

To execute the script, use the following command with the required arguments:

```bash
bash script_name.sh <DAPFILE> <SBAMFILE> <GENE> <OUTDIR>
```

Ensure the appropriate input DAP file, SBAM file, gene name, and output directory are specified. The script runs PTWAS Builder to create the weights for the given gene and stores them in the specified output directory.

---

## Script 02                                   
## 02_concat_weights.sh

This script concatenates multiple PTWAS weight files into a single file and compresses the final output. It combines all individual gene-level PTWAS weights into one comprehensive file for easier analysis.

## Command Line Arguments

- `WEIGHTSDIR`: The directory containing individual PTWAS weight files.
- `OUTDIR`: The directory where the final concatenated PTWAS weights file will be saved.

## Execution Steps

1. **Gather PTWAS Files**:
   - Retrieve all files in the specified weights directory that contain PTWAS weights (files with "_ptwas_weights" in the name).

2. **Create Output File**:
   - Create a new output file to store all concatenated PTWAS weights.
   - If a previous output file exists, delete it to ensure a clean slate.

3. **Concatenate PTWAS Files**:
   - Loop through the individual PTWAS weight files and append their contents to the output file.

4. **Compress the Final Output**:
   - If a compressed output file exists, remove it to avoid conflicts.
   - Use `gzip` to compress the final concatenated file for efficient storage.

## Running the Script

To execute the script, use the following command with the specified arguments:

```bash
bash script_name.sh <WEIGHTSDIR> <OUTDIR>
```

Ensure the correct directory paths are provided for the weights and output directories. The script concatenates all individual PTWAS weight files into one comprehensive file and compresses it to reduce storage space.

---

## Script 03                                    
## 03_run_make_gambit.sh

This script creates a VCF file from PTWAS weights and formats the output to be compatible with GAMBIT for further processing. It also compresses and indexes the VCF file using `bgzip` and `tabix`.

## Command Line Arguments

- `PTWAS_WEIGHTS_FILE`: The input file containing PTWAS weights to be converted into VCF format.
- `OUTDIR`: The output directory where the formatted VCF files will be stored.

## Execution Steps

1. **Conda Environment Activation**:
   - Activate the appropriate conda environment with R and other necessary tools.

2. **Run R Script**:
   - Use `Rscript` to convert PTWAS weights into a VCF file using an external script, `make_GAMBIT.DB.R`.
   - The script outputs the VCF to the specified location.

3. **Format Output by Chromosome**:
   - Loop through each chromosome (1-22) and create a separate VCF file for each.
   - Add a VCF header and process each chromosome's data by removing "chr" from chromosome names and IDs.
   - Append the formatted lines to the appropriate output file for each chromosome.

4. **Compress and Index Output**:
   - Use `bgzip` to compress the VCF files.
   - Use `tabix` to index the compressed VCF files, allowing for quick access to data by chromosome.

---

## Script 04                                    
## 04_format_gwas.R

This R script processes GWAS data and formats it for use with PTWAS and FastENLOC. It includes various tasks such as data transformation, alignment, and compression.

## Functionality Overview

The script performs the following tasks:

1. **Load and Join GWAS Data**:
   - Load GWAS data and map it to variant names using specified files.
   - Retrieve original GWAS data to join with additional information.

2. **Estimate Effect Sizes**:
   - Compute standard errors (SE) and betas from z-scores and p-values.
   - Plot relationships between z-scores and betas to ensure accurate estimations.

3. **Add Minor Allele Frequencies (MAFs)**:
   - Join minor allele frequencies with GWAS data for further computations.

4. **Prepare Final GWAS Data**:
   - Separate and format columns for chromosome, position, reference, and alternative alleles.
   - Reorder by chromosome and position.
   - Compress and index the VCF file using `bgzip` and `tabix`.

5. **Check Overlaps with QTL**:
   - Optionally, check that variants overlap with QTL data for validation.

## Script Execution

### Running the Script

To execute the script, use an R environment with necessary libraries installed. Ensure that all required data paths and directories are properly set.

### Output Files

The script outputs formatted GWAS data, which can be used for subsequent analysis with PTWAS and FastENLOC. It creates VCF files for further processing and compresses and indexes them for easy storage and access.

### Notes

- Ensure all input file paths are correct before running the script.
- The script includes optional checks for variant overlaps with QTL, which can be useful for validation.
- Be cautious with large datasets, as operations on GWAS data can be resource-intensive.

---

## Script 05                                    
## 05_check_ld_panel.R

This R script processes LiftOver to check if there are matches between the GWAS data and LD (Linkage Disequilibrium) panel after the LiftOver from hg19 to hg38. It also generates outputs required for subsequent analyses.

## Functionality Overview

The script carries out the following tasks:

1. **Load and LiftOver LD Panel**:
   - Import the LD panel in VCF format.
   - Use LiftOver to map coordinates from hg19 to hg38.
   - Handle duplicates and mismatches created by LiftOver.

2. **Filter and Sort the Lifted LD Panel**:
   - Remove duplicate records that arise from LiftOver.
   - Filter and order by chromosome and position.
   - Ensure consistent formatting.

3. **Match LD Panel with GWAS**:
   - Join LiftOver data with GWAS to determine overlapping SNPs.
   - Identify and handle swapped reference and alternate alleles.
   - Output relevant data for further analysis and diagnostics.

4. **Create Files for PLINK**:
   - Generate files to correct swapped reference/alternate alleles using PLINK.
   - Format the output to align with PLINK conventions.

## Script Execution

### Parameters

- `chromosome`: The chromosome to process.
- `gwas_type`: The type of GWAS to process.
- `savedir`: The directory where the outputs will be saved.

### Running the Script

Ensure the required R libraries are installed. The script reads specific data files, so verify that the paths are correct. Execute the script in an R environment with the necessary permissions to read and write files.

### Output Files

- Lifted LD panel in VCF format.
- Files for swapping reference and alternate alleles for PLINK.
- Outputs showing matched, mismatched, and swapped variants.

### Notes

- The script includes detailed comments and intermediary checks to ensure accurate processing.
- Be cautious when handling large datasets, especially when performing LiftOver operations.
- Proper error handling and resource management are crucial for seamless execution.


---

## Script 06                                    
## 06_swap_ld_panel.sh

This script uses PLINK to swap reference and alternate alleles for a subset of variants in the LD (Linkage Disequilibrium) panel. This operation aligns the LD panel with other genomic data, such as GWAS and QTL annotations.

## Functionality Overview

The script accomplishes the following:

1. **Swap Reference/Alternate Alleles**:
   - Use PLINK to update alleles in the LD panel according to a specified swap file.
   - Ensure updated variants match with GWAS and QTL annotations.

2. **Format and Compress Updated LD Panel**:
   - Convert phased genotypes from "x/y" format to "x|y".
   - Compress the updated VCF files using `bgzip`.
   - Create an index for the compressed files using `tabix`.

3. **Create Output Directories and Files**:
   - Generate output directories for storing updated LD panels.
   - Save the updated LD panel in compressed VCF format.

## Script Execution

### Parameters

- `CHROM`: Chromosome numbers to process.
- `GWAS_TYPE`: Type of GWAS data.
- `LD_VCF`: Path to the original LD panel in VCF format.
- `SWAPFILE`: Path to the file containing variants that need to be updated.
- `OUTDIR`: Output directory for the updated LD panel.
- `OUTFILE`: Base name for the output files.
- `TMPFILE`: Temporary file used for intermediate operations.

### Running the Script

Ensure PLINK and `bgzip` are installed and available in your environment. The script should be executed in a Unix-like environment with proper permissions for reading and writing files.

### Output Files

- Updated LD panel in compressed VCF format.
- Index file created by `tabix`.
- Intermediary files (temporary text file).

### Notes

- This script processes chromosomes in a loop, updating alleles and reformatting the LD panel.
- Proper error handling and resource management are crucial for seamless execution.
- Ensure that the paths and directory structures are correct before running the script.

---

## Script 07                                    
## 07_ptwas_scan.sh

This script runs a PTWAS (Probabilistic Transcriptome-Wide Association Study) scan using the GAMBIT software. It combines GWAS data, PTWAS weight files, and LD panel files to perform a PTWAS analysis.

## Functionality Overview

The script performs the following operations:

1. **PTWAS Scan**:
   - Runs the PTWAS scan using GAMBIT, with provided GWAS, PTWAS weight, and LD panel files.
   - Outputs the results to a specified location.

2. **Configuration of Input and Output Files**:
   - Configures the file paths for the GWAS data, PTWAS weights, and LD panels.
   - Sets the output prefix for the results.

## Script Execution

### Parameters

- `GWAS_FILE`: Path to the GWAS input file.
- `PTWAS_WEIGHT_FILE`: Path to the PTWAS weight file.
- `LD_PANEL_FILES`: Path to the LD panel input files.
- `PREFIX`: Prefix for the output files.

### Output

- PTWAS scan results with the specified prefix.
- Intermediate output indicating the files used and the script's operations.

### Running the Script

Ensure that GAMBIT is installed and accessible in your environment. To execute the script:

- Set the correct file paths for the input data.
- Define the output directory and prefix for the results.
- Run the script in a Unix-like environment with proper permissions.

### Notes

- Before running the script, verify the correctness of input file paths.
- Check that the necessary dependencies and tools, such as GAMBIT, are installed.
- Ensure you have sufficient compute resources, as PTWAS scans can be resource-intensive.

---

## Script 08                                    
## 08_format_scan_results.R

This script checks the results of PTWAS (Probabilistic Transcriptome-Wide Association Study) scans for various runs, summarizing and stratifying the PTWAS output.

## Overview

- Processes multiple PTWAS runs, specified by combinations of summary type, region type, and cell type.
- Reads the summary and stratified output from PTWAS scans.
- Combines and summarizes the results across multiple runs and chromosomes.
- Outputs the results to summary and stratified files.

## Script Execution

### Input Data

- `summary_types`: Types of summaries for analysis.
- `region_types`: Regions of interest in the analysis.
- `cell_types`: Cell types to consider.
- `proportion_type`: Specifies the proportion used for PTWAS analysis.

### Configuration

- Creates data frame with run details, including region, cell, and summary types.
- Sets the GWAS type to "wightman".

### Processing Runs

- Defines a process for reading summary and stratified PTWAS results for each chromosome.
- Combines the results across all chromosomes for each run.
- Provides an option to check for missing results or genes.

### Output

- Combined summary and stratified PTWAS results.
- Optional report on missing results or genes.

### Execution Details

1. **Process Runs**:
   - Iterates over the list of runs, processing each chromosome individually.
   - Reads the summary and stratified output files.
   - Combines results across chromosomes.

2. **Output Results**:
   - Saves the combined summary and stratified results to output files.
   - Optionally checks for missing genes in the results.

### Notes

- Before running the script, ensure that the PTWAS output files are correctly generated and located.
- This script is designed for a specific GWAS type ("wightman") but can be modified to accommodate other types.
- Ensure that the directory structure for saving results exists or is created as needed.


---

## Script 11                                    
## 11_run_intact.R

We skip scripts 9 & 10 because they weren't used in the final manuscript.

This script runs INTACT (Integration of Transcriptome-Wide Association Study and Colocalization) on PTWAS (Probabilistic Transcriptome-Wide Association Study) and colocalization output. It processes and visualizes the results, providing insights into potential causal genes.

## Overview

- Loads PTWAS and colocalization data for specified runs.
- Joins PTWAS and colocalization data.
- Runs the INTACT pipeline to calculate posterior inclusion probabilities (PIPs) for potential causal genes.
- Outputs significant results and plots them.
- Checks and adds gene annotations.
- Provides plot functions for visualization of significant results.

## Script Execution

### Input Data

- `summary_types`: Types of summaries for analysis.
- `region_types`: Regions of interest in the analysis.
- `cell_types`: Cell types to consider.
- `gwas_type`: Type of genome-wide association study (GWAS) data to use.
- `savedir`: Directory to save output results.

### Functions

#### `load_ptwas(runname, gwas_type)`

- Loads PTWAS data for a specific run.
- Returns a data frame with PTWAS results.

#### `load_fastenloc(runname, gwas_type)`

- Loads colocalization results from FastEnloc.
- Returns a data frame with colocalization results.

#### `run_intact(ptwas, fastenloc, runname)`

- Joins PTWAS and FastEnloc data.
- Runs INTACT to calculate PIPs.
- Outputs significant results and saves them to a file.
- Checks for gene annotations and attaches them to the results.
- Returns the formatted results.

#### `plot_results(res_df, region_type)`

- Plots significant results with a specified PIP threshold.
- Annotates known and novel genes.
- Saves the plot to an output directory.

### Execution Details

1. **Process Runs**:
   - Iterates over a list of runs, loading PTWAS and colocalization data.
   - Runs INTACT to calculate PIPs for each run.
   - Saves the significant results to a file.

2. **Generate Plots**:
   - Plots significant results with an adjustable PIP threshold.
   - Annotates known Alzheimer's disease (AD) genes.
   - Saves the plot to an output directory.

### Notes

- Ensure that all necessary input files are accessible and correctly formatted.
- Adjust the `summary_types`, `region_types`, `cell_types`, and `gwas_type` as needed.
- Check the directory structure to ensure it exists or create it as needed.
