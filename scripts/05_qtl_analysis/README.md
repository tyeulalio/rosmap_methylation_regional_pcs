# Quantitative Trait Loci (QTL) Analysis

This pipeline contains a set of functions designed to perform QTL analysis and analyze the results.

## Table of Contents

- [01_liftover_vcf.sh](#script-01)
- [02_run_plink.sh](#script-02)
- [03_format_pc_beds.R](#script-03)
    - [load_methylation](#load_methylation)
    - [load_annotations](#load_annotations)
    - [convert_ids](#convert_ids)
    - [create_bedfile](#create_bedfile)
    - [save_bed](#save_bed)
    - [main](#main1)
- [04_mishap_test.sh](#script-04)
- [05_filter_mishap.R](#script-05)
- [06_merge_bims.sh](#script-06)
- [07_eigenstrat_pca.sh](#script-07)
- [08_format_eigenstrat.R](#script-08)
- [09_remove_outliers.sh](#script-09)
- [10_compute_global_pcs.R](#script-10)
- [11_format_qtl_input.R](#script-11)
- [12_map_qtls.sh](#script-12)


## Script 01                                    
## 01_liftover_vcf.sh

This bash script is intended to be used for the liftover of ROSMAP Whole Genome Sequencing (WGS) VCFs from the hg17 to hg19 genome build. It employs GATK's LiftoverVCF utility to remap the coordinates of variants contained within VCF files according to the specified chain file, ensuring that variant annotations are accurate with respect to the updated reference genome.

### Prerequisites

Before running this script, make sure to load the necessary modules:
- GATK/4.0.10.0 for the liftover tool
- bcftools for VCF file manipulations

### Temporary Files and Directories

The script utilizes a temporary directory to store intermediate files during processing:

- `tmp_files`: Directory to store temporary files.

The following temporary files are used:

- `tmp_concat_chroms.vcf`: Stores concatenated chromosome files.
- `tmp.vcf`: Temporary storage for VCF file with 'chr' added to chromosome numbers.
- `tmp_norm.vcf`: Normalized VCF file before liftover.

### Execution Steps

1. **Concatenating Chromosome VCFs**: Concatenate VCFs of all chromosomes into a single file, as some variants may map to different chromosomes after the liftover.

2. **Storing Temporary VCF Data**: If selected, store the modified VCF with 'chr' prefixing chromosome numbers to align with the reference fasta file.

3. **Normalization**: If selected, the VCF file is normalized using bcftools norm to improve match rates during the liftover.

4. **Liftover**: If selected, GATK LiftoverVcf is executed to map the variants from hg17 to hg19.

### Script Flags

The following flags can be toggled to enable or disable certain steps of the script:

- `STORE_DATA`: If set to true, stores the modified VCF with 'chr' prefixing chromosome numbers.
- `NORMALIZE`: If set to true, normalizes the VCF file before liftover.
- `LIFTOVER`: If set to true, performs the liftover from hg17 to hg19.

### Post-Liftover Normalization

After the liftover, the script performs a final normalization step on the hg19 VCF file using bcftools.

## Script Execution

To execute the script, run it from the command line as follows:

```bash
bash liftover_vcf.sh
```

Make sure all the paths to the input files, chain files, and reference genomes are correctly set within the script before execution. It's also important to ensure that the necessary compute resources are available, as genome-wide liftover operations are resource-intensive.

---

## Script 02                                  
## 02_run_plink.sh

# Plink2 Data Conversion Script Documentation

This bash script uses Plink2 to convert normalized WGS VCF files into Plink's native file formats, including .pgen (Plink 2 binary genotype table), .bed (binary genotype/phenotype table), and .vcf (Variant Call Format).

### Prerequisites

Ensure Plink2 is loaded in your environment:

```bash
module load plink2
```

### Execution Steps

1. **Pgen File Creation**: Generates Plink2-specific .pgen files from the VCF.

2. **Bed File Creation**: Converts VCF files to .bed files, commonly used in genomic analyses for their compact size and computational efficiency.

3. **VCF File Creation for Individual Samples**: Exports VCF files for a subset of samples, typically for quality control or specific analyses.

### Script Flags and Parameters

- `--chr`: Specifies the chromosome numbers to include.
- `--output-chr chrM`: Sets the output format for chromosome names.
- `--max-alleles 2`: Limits alleles per variant to two, focusing on biallelic variants.
- `--mind`, `--hwe`, `--geno`, `--maf`: Set various thresholds for sample missingness, Hardy-Weinberg equilibrium, genotype missingness, and minor allele frequency, respectively.
- `--set-all-var-ids`: Formats variant IDs with a custom string pattern.
- `--new-id-max-allele-len`: Sets the maximum length for allele representation in variant IDs.
- `--vcf-half-call`: Determines how to handle half-calls in the VCF.

### Output Files

- `.pgen`, `.pvar`, and `.psam`: Set of files comprising the Plink2 binary genotype table.
- `.bed`, `.bim`, and `.fam`: Binary files representing genotype/phenotype data for Plink.
- `.vcf`: Variant call format files for individual samples.

### Quality Control

- Samples are filtered based on inclusion criteria defined in `keep_wgs_samples.tsv`.
- Variants are filtered based on missingness, Hardy-Weinberg equilibrium, genotype missingness, and minor allele frequency thresholds.

## Script Execution

To run the script, execute it from the command line:

```bash
bash plink2_data_conversion.sh
```

Ensure the paths to the VCF files and the sample keep-lists are correct. The script assumes a specific directory structure and naming conventions which should be adjusted to fit the project's setup.


---

## Script 03                                    
## 03_format_pc_beds.R

This R script is tailored for processing methylation data to prepare BED files for principal component-based analyses across various cell types and genomic regions. It incorporates several specialized functions for effective data manipulation and analysis, particularly focusing on creating formatted BED files that integrate methylation data and annotations. The key stages of this script include:

1. **Data Loading**: Utilizes functions like `load_methylation` to fetch methylation data, and `load_annotations` to retrieve region-specific genomic annotations, aligning the data with the necessary context for further processing.

2. **ID Conversion**: The function `convert_ids` adjusts sample IDs to ensure consistency across different data types, facilitating accurate data integration and comparison.

3. **Data Formatting**: Converts methylation data into a BED file format through `create_bedfile`, adapting it for genomic analyses that require specific data structuring.

4. **File Management**:
   - `save_bed` function is responsible for saving the formatted data into BED files, managing file paths, and compressing the final outputs for efficient storage and handling.

5. **Main Function**: Directs the execution of the pipeline by orchestrating the function calls in a sequential manner based on the predefined parameter sets. It ensures that each combination of cell type, region type, and summary type is processed systematically.

6. **Execution Control**: Implements a flexible execution environment that allows for batch processing through command line arguments. This design aids in the scalable execution across a computing cluster, making the script adaptable to various computational and data contexts.

7. **Utility Functions**:
   - Functions like `head()` and `dim()` are frequently used for data inspection and debugging throughout the script to ensure data integrity and correct processing steps.

The script is built with an emphasis on modularity and practical utility, making it highly reusable for different datasets and adaptable to evolving research needs. It is especially useful in high-throughput settings where efficient data processing and summarization are critical.

Below are detailed descriptions and usage guidelines for each of the functions employed in this script.

---

### <a name="load_methylation">`load_methylation`</a>

#### Description

The `load_methylation` function is designed to load and preprocess methylation data for specific cell types, regions, and summary methods. It adjusts the loading process based on whether methylation data needs to be linked to CpG sites or summarized by principal components (PCs) or averages. The function is structured to handle different types of methylation data formats, providing flexibility in handling various methylation summary types such as CpGs, averages, and PCs.

#### Parameters

- `cell_type`: A character string specifying the cell type for which methylation data is being loaded.
- `region_type`: A character string indicating the genomic region of interest (e.g., `full_gene`, `promoters`).
- `summary_type`: A character string defining the summary method used (`cpgs`, `avgs`, `pcs`).
- `match_cpgs`: A boolean indicating whether to match CpG sites to their respective regions. This is particularly relevant when summary_type is 'cpgs'.

#### Returns

- `meth`: A dataframe of methylation data prepared for further analysis. The structure and content of the dataframe vary based on the summary type specified.

#### Details

- The function begins by determining the correct file path and loading the methylation data from an RDS file based on input parameters.
- If `summary_type` is 'avgs', unnecessary columns such as `duration_secs` are removed.
- For `summary_type` 'pcs', it ensures that only principal component related columns are retained unless specific region types require the removal of failed genes or unnecessary metrics.
- When `match_cpgs` is true and `summary_type` is 'cpgs', it performs additional filtering to align CpGs with their respective gene regions based on a predefined mapping file.

#### Example Usage

```R
# Load methylation data for neuronal cells summarized by CpGs in gene bodies
neuron_meth_data <- load_methylation("neuron", "gene_body", "cpgs", TRUE)
```

---

### <a name="load_annotations">`load_annotations`</a>

#### Description

The `load_annotations` function is designed to load and format genomic region annotations specific to different types of analysis. It processes annotations depending on the region type (e.g., `full_gene`, `promoters`). This function filters and formats these annotations to suit various downstream genomic analyses, such as setting transcription start site (TSS) positions based on gene orientation and filtering out unnecessary chromosome data.

#### Parameters

- `region_type`: A character string indicating the type of genomic region for which annotations are being loaded.

#### Returns

- `tss_annots`: A dataframe with formatted annotations including TSS positions, which are essential for defining QTL analysis windows.

#### Details

- The function begins by constructing a file path based on the `region_type` and loading the relevant annotation data.
- For `full_gene`, it formats the annotations to exclude chromosomes X and Y, selects regions from the promoter to the 5' UTR, and calculates the TSS based on gene strand.
- For `promoters`, it calculates the center position of each region.
- If the `region_type` is not specifically handled, it defaults to loading TSS positions prepared for the `full_gene` category.
- Each annotation type undergoes specific transformations to prepare it for downstream applications, such as QTL window definitions.

#### Example Usage

```R
# Load and format gene annotations for full_gene region types
full_gene_annotations <- load_annotations("full_gene")
```

---


### <a name="convert_ids">`convert_ids`</a>

#### Description

The `convert_ids` function is designed to synchronize sample identifiers in methylation datasets with those in a genotype database. It ensures that methylation data samples match the identifiers used in genotype files, facilitating combined genetic and epigenetic analyses.

#### Parameters

- `meth`: A data frame or list containing methylation data where column names represent sample IDs.

#### Returns

- A list containing:
  - `meth`: The methylation data with updated sample IDs.
  - `clinical`: The ordered clinical metadata corresponding to the new sample IDs in the methylation dataset.

#### Details

- This function starts by extracting sample IDs from the methylation data.
- It then loads matched clinical metadata from a predefined path.
- The clinical data is reordered to match the order of sample IDs in the methylation dataset.
- Sample IDs in the methylation data are replaced with corresponding specimen IDs from the genotype data to ensure consistency across datasets.
- It checks and confirms that the order of the clinical data matches the methylation sample IDs after conversion.

#### Example Usage

```R
# Load methylation data and convert sample IDs
meth_data <- load_meth(...)  # Assume meth_data is already loaded
converted_data <- convert_ids(meth_data)
meth_updated <- converted_data$meth
clinical_info <- converted_data$clinical
```

---

### <a name="create_bedfile">`create_bedfile`</a>

#### Description

The `create_bedfile` function is tailored for converting methylation data into BED file format suitable for genomic analyses. This function integrates annotations with methylation data and handles different summarization types (`avgs`, `pcs`, `cpgs`). It ensures that each methylation measurement is appropriately aligned with genomic coordinates and annotations for downstream processing.

#### Parameters

- `converted_dat`: A list containing normalized and converted methylation data (`meth`) and potentially other associated data.
- `annots`: Annotation data that includes genomic coordinates and other relevant genomic features.
- `summary_type`: A character string indicating the summarization method used (`avgs` for averages, `pcs` for principal components, `cpgs` for CpG sites).

#### Returns

- A data frame formatted as a BED file with columns for chromosome (`chrom`), start position (`start`), end position (`end`), and additional columns for each sample's phenotype data.

#### Details

- Initially, the function extracts the methylation data and removes any rows with missing values.
- Depending on the `summary_type`:
  - **Avgs**: Processes average methylation values and merges them with annotations.
  - **Pcs**: Splits principal component identifiers from methylation data, merges with annotations, and retains unique identifiers for each principal component.
  - **Cpgs**: Processes CpG-specific data, aligns CpG identifiers with genomic coordinates, and adjusts the end positions to represent the length of the CpG.
- All types involve normalization of methylation data, but the example omits the details of normalization for simplicity.
- The function concludes by formatting the merged data into a standard BED format.

#### Example Usage

```R
# Assuming 'converted_data' and 'annotations' are pre-loaded
bed_data <- create_bedfile(converted_data, annotations, 'avgs')
write.table(bed_data, "output_path.bed", quote=FALSE, sep="\t", row.names=FALSE)
```

---

### <a name="save_bed">`save_bed`</a>

#### Description

The `save_bed` function is designed to save the formatted BED file to disk, followed by compressing the file using gzip. It constructs the filename based on the cell type, region type, summary type, and additional parameters, then saves the file in a specified directory. This function ensures data is stored efficiently and is readily accessible for further genomic analysis or archiving.

#### Parameters

- `formatted_bed`: A data frame containing the formatted BED data ready for saving.
- `cell_type`: A character string indicating the cell type to be included in the filename.
- `region_type`: A character string indicating the region type to be included in the filename.
- `summary_type`: A character string indicating the summary type to be included in the filename.

#### Returns

- No explicit return value. The function's primary purpose is to save the data to disk.

#### Details

- The function constructs the file path using the specified `cell_type`, `region_type`, `summary_type`, and a directory path stored in `savedir`.
- It uses `write_tsv` from the `readr` package to save the BED file as a tab-separated values file.
- After saving, the function checks for any existing gzip version of the file and removes it if found.
- Finally, it compresses the new BED file using `gzip`, optimizing storage space.

#### Example Usage

```R
# Assuming 'formatted_bed' is pre-loaded and cell_type, region_type, summary_type are defined
save_bed(formatted_bed, "neuron", "promoters", "avgs")
```

---


### `main`

#### Description

The `main` function runs a comprehensive pipeline to process methylation data for different cell types, region types, and summary types. It sequentially executes the processes of loading methylation data, converting sample IDs, formatting data into BED files, and saving them. The function is designed to handle both specific regions and a genome-wide approach for CpGs, ensuring a versatile and robust analysis flow.

#### Execution Steps

1. **Parameter Initialization**: Defines the combinations of cell types, region types, and summary types for which the analysis will be performed.
2. **Data Processing Loop**:
   - For each combination of parameters, the following steps are executed:
     - Load methylation data.
     - Load region-specific annotations.
     - Convert sample IDs to match external genotype data.
     - Format the data into BED file format.
     - Save the formatted data to disk.
3. **Genome-wide CpGs Processing**:
   - Additionally, for genome-wide CpGs, a separate process run is executed for each cell type to handle comprehensive genomic analysis.

#### Parameters

- No direct parameters are passed to the function. All necessary parameters are defined within the function body.

#### Details

- The function uses nested processing functions to handle each step, ensuring modular and traceable execution.
- Output is primarily in the form of files, specifically BED files for subsequent genomic analyses.
- Logging statements are included to provide real-time feedback during processing.

#### Example Usage

```R
# Execute the main function
main()
```

This function represents the core of the methylation analysis pipeline, integrating multiple data handling steps into a seamless workflow, ideal for high-throughput genomic data processing.

---

## Script 04                                    
## 04_run_mishap_test.sh

This bash script is designed to perform mishap testing using PLINK on genomic data. The script focuses on analyzing data from Whole Genome Sequencing (WGS) that has been harmonized and formatted into PLINK's binary file format (.bed, .bim, .fam). It utilizes the mishap test function of PLINK to evaluate potential mishaps or errors in the genetic data. Key components and features of this script include:

### Prerequisites

Before executing this script, ensure that the PLINK module is loaded in your environment:

```bash
module load plink
```

### Execution Steps

1. **Data Specification**: Specifies the input binary file (`--bfile`) located at `/path/to/data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms` which contains the genotype data prepared in PLINK binary format.

2. **Mishap Testing**: Executes the mishap test (`--test-mishap`) to check for inconsistencies or errors in genotype data across all chromosomes.

### Output

- The output is directed to a designated path `../../output/05_qtl_analysis/04_mishap_test/all_chroms_hg38`, which will contain the results of the mishap tests.

### Script Execution

To run the script, use the following command from the terminal:

```bash
bash run_mishap_test.sh
```

Ensure that the script is executed in an environment where the necessary computational resources and file access permissions are available. Adjust the paths and file names as necessary to match your project's directory structure and naming conventions.

### Considerations

- Verify that the input data files are correctly formatted and located at the specified paths.
- The mishap test is computationally intensive; ensure that sufficient computational resources are available.

This script is designed to integrate seamlessly into a larger genomic analysis workflow, providing essential checks for data quality and consistency.

---

## Script 05                                    
## 05_filter_mishap.R

This R script is utilized to filter out genetic variants and manage sample inclusion based on results from mishap tests and sample metadata. The script operates on data formatted in previous steps and outputs files for further genomic analysis, emphasizing data integrity and relevance.

### Libraries and Dependencies

The script uses the `data.table` and `tidyverse` libraries for data manipulation and I/O operations:

```R
library(data.table)
library(tidyverse)
```

### Data Processing and Variant Filtering

1. **Load Mishap Test Results**: Loads haplotype test results from a previous mishap test to identify variants that may be problematic.

   ```R
   hap <- read_table("../../output/05_qtl_analysis/04_mishap_test/all_chroms_hg38.missing.hap")
   ```

2. **Filter Significant Variants**: Filters out variants based on a significance threshold determined from xQTL studies.

   ```R
   sig <- hap %>% filter(P < 1e-9)
   ```

3. **Save Excluded Variants**: Exports a list of excluded variant IDs to a file for further use in genomic analyses.

   ```R
   write_tsv(exclude_df, "../../output/05_qtl_analysis/05_filter_mishap/exclude_variants_hg38.txt", col_names=FALSE)
   ```

### File Management and Sample Inclusion

1. **Create File Lists for Merging**: Generates text files listing .bim files for merging in subsequent analyses.

2. **Sample Filtering**: Matches sample IDs to those in methylation data, ensuring consistency across genetic and epigenetic data sets.

3. **Output Filtered Sample Lists**: Outputs lists of samples to include and exclude from further analysis, aiding in the coordination of sample data across multiple studies or batches.

### Execution

Run the script from the command line or within an R environment, ensuring all file paths and library dependencies are correctly configured:

```R
source("05_filter_variants_and_samples.R")
```

### Output Files

- `exclude_variants_hg38.txt`: List of variant IDs to be excluded.
- `bim_files_hg38.txt`: List of .bim files for merging.
- `keep_samples.txt`: List of samples to keep.
- `remove_samples.txt`: List of samples to remove.

### Considerations

- Verify all paths and file names correspond to the actual directory structure and naming conventions used in your project.
- Ensure that the threshold for filtering variants is set appropriately based on previous results and scientific rationale.

---

## Script 06                                    
## 06_merge_bims.sh

This bash script uses Plink and Plink2 to merge genetic data files and export them in various formats suitable for downstream genetic analyses. The script facilitates the management of genetic data by merging files and creating derived datasets, ensuring they contain only relevant variants and samples.

### Prerequisites

Before running this script, ensure that Plink and Plink2 are loaded in your environment:

```bash
module load plink
module load plink2
```

### Script Execution Steps

1. **Merging Files**: Merges genotype files while excluding certain variants and keeping specific samples.
   
   ```bash
   plink2 --make-pgen --pfile "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms" --exclude "../../output/05_qtl_analysis/05_filter_mishap/exclude_variants_hg38.txt" --keep "../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt" --out "../../output/05_qtl_analysis/06_merge_bims/all_chroms_hg38"
   ```

2. **Export to BED Format**: Converts the merged genotype data to BED format for compatibility with a wide range of bioinformatics tools.
   
   ```bash
   plink2 --make-bed --pfile "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms" --exclude "../../output/05_qtl_analysis/05_filter_mishap/exclude_variants_hg38.txt" --keep "../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt" --out "../../output/05_qtl_analysis/06_merge_bims/all_chroms_hg38"
   ```

3. **Creating Allele Counts**: Generates a file containing counts of the reference allele for each variant.
   
   ```bash
   plink2 --pfile "../../output/05_qtl_analysis/06_merge_bims/all_chroms_hg38" --export A-transpose --out "../../output/05_qtl_analysis/06_merge_bims/all_chroms_ref_counts_hg38"
   ```

4. **Export to PED Format**: Converts the merged data into PED format, a plain text format widely used for representing genotypic data alongside phenotypic data.

   ```bash
   plink2 --pfile "../../output/05_qtl_analysis/06_merge_bims/all_chroms_hg38" --export ped --out "../../output/05_qtl_analysis/06_merge_bims/all_chroms_ref_counts_hg38"
   ```

### Output Files

- `.pgen`, `.pvar`, and `.psam`: Files created in Plink2 binary genotype table format.
- `.bed`, `.bim`, and `.fam`: BED files representing genotype/phenotype data.
- `all_chroms_ref_counts_hg38`: Contains allele counts for each variant.
- `all_chroms_hg38.ped`: PED file containing genotype alongside phenotype data.

### Running the Script

Execute this script from the command line within a suitable compute environment:

```bash
bash merge_and_export_plink_files.sh
```

Ensure the paths to the input files, exclusion lists, and sample keep-lists are correct. Adjustments may be necessary to align with your project's specific directory structure and file naming conventions. This script assumes a particular workflow that should be validated against your data processing requirements.

---

## Script 07                                    
## 07_eigenstrat_pca.pl

This Perl script orchestrates the execution of the `smartpca` tool from the EIGENSOFT package to perform principal component analysis (PCA) on genetic data. It is tailored to handle large-scale genomic datasets and is particularly useful for identifying and visualizing population substructure.

### Prerequisites

Ensure that the `smartpca` executable from EIGENSOFT is correctly installed and included in the system PATH:

```bash
export PATH="/path/to/programs/EIG-7.2.1/bin:$PATH"
export PATH="/path/to/programs/EIG-7.2.1/src/eigensrc:$PATH"
```

### Script Execution Steps

1. **Set Environment Variables**: Adjust the PATH to include directories containing the `smartpca` executable and related scripts.

2. **Define Input and Output Files**: Specify paths for input files (genotype, SNP, and individual data) and output files (PCA results, plots, and logs).

3. **Execute SmartPCA**:
   - `-i`: Path to the BED file generated by Plink.
   - `-a`: Path to the MAP file corresponding to the BED file.
   - `-b`: Path to the PED file for individual IDs.
   - `-k`: Number of principal components to compute.
   - `-o`, `-p`, `-e`, `-l`: Paths for output PCA results, plots, eigenvalues, and log files, respectively.
   - `-m`, `-t`, `-s`: Options for maximum iterations, number of PCs for outlier removal, and standard deviations for outlier definition.

### Running the Script

Execute the script from the command line with appropriate paths and settings:

```bash
perl run_smartpca.pl
```

### Output Files

- `.pca`: File containing the PCA scores for each individual.
- `.plot`: File suitable for plotting PCA results.
- `.eigenvalues`: Contains the eigenvalues for each principal component, indicating the variance explained.
- `.log`: Log file containing details of the PCA execution process.

### Additional Considerations

- Make sure the file paths and directory structure are set correctly according to your project's organization.
- Review the PCA settings and adjust the number of principal components and iteration settings as necessary to match the specifics of your analysis.
- Validate the PCA results by inspecting the plots and eigenvalues to ensure that the analysis captures significant genetic variation and substructure.

This script is integral for genetic studies aiming to account for population structure, which is crucial for downstream genetic association studies.

---

## Script 08                                  
## 08_format_eigenstrat.R

This R script processes Eigenstrat output for QTL analysis by formatting principal components (PCs), computing variance, and identifying and removing outliers. It uses various tidyverse functions, particularly from the `ggplot2` and `dplyr` libraries, to manipulate and visualize data effectively. The script organizes and prepares data for downstream genetic and clinical trait correlation analyses.

### Processing Steps

1. **Variance Calculation**: Computes the total and percentage variance explained by the PCs, formatting the results into a structured summary table.
  
2. **PC Formatting**: Transforms the principal components data to match clinical metadata, ensuring consistency across datasets.

3. **Outlier Identification**: Detects and removes samples that are outliers in the PC analysis, which could bias further analyses.

4. **Metadata Matching**: Aligns processed PC data with metadata, filtering to include only relevant samples.

5. **Correlation Analysis**: Performs correlation tests between PCs and various clinical traits, adjusting for multiple comparisons to identify significant relationships.

### Key Functions

- `dir.create()`: Ensures the necessary directories exist.
- `read_table()`: Loads data from specified file paths.
- `write_tsv()`: Saves data frames to disk in TSV format.
- `saveRDS()`: Serializes R objects to disk.

### Output Files

- `pca_pct_var_hg38.tsv`: Contains the percentage of variance explained by each PC.
- `genotype_pcs_hg38.rds`: Serialized R data frame of formatted principal components.
- `population_outliers_hg38.txt`: List of sample IDs identified as outliers.
- `genotype_pc_trait_lm_pvals_hg38.rds`: Results of the correlation tests between PCs and clinical traits.
- Plots visualizing correlation p-values, highlighting significant findings.

### Example Usage

To run the script, ensure that all file paths and directories are correctly set up as specified in the script. Load the necessary libraries and execute the script within an R environment. This script assumes that preliminary PCA has been conducted and that the necessary input files are available in the specified directories.

```bash
Rscript format_eigenstrat.R
```

Ensure the script is executed in a directory structure that matches the expected input paths, or modify the paths in the script to match your environment. This script is part of a larger analysis pipeline and relies on outputs from previous steps, so it should be run sequentially after those steps have completed.

---

## Script 09                                  
## 09_remove_outliers.sh

This bash script uses Plink2 to refine genotype data by removing identified outliers from a dataset. 

### Prerequisites

Before running this script, ensure that Plink2 is available and loaded into your environment:

```bash
module load plink2
```

### Execution Steps

1. **Input Data**: The script processes genotype data stored in Plink2's `.pfile` format, which is specified along with its directory path.

2. **Remove Outliers**: It utilizes a previously generated list of outlier sample IDs to exclude these from the analysis, enhancing the dataset integrity.

3. **Output Data**: The script outputs a new set of Plink binary files (`.bed`, `.bim`, `.fam`) without the outliers, prepared for further genetic analysis.

### Script Flags and Parameters

- `--make-bed`: Command to create a binary fileset from Plink2 files.
- `--output-chr chrM`: Ensures that the output files maintain chromosome nomenclature consistent with mitochondrial DNA labeling.
- `--pfile`: Specifies the path to the input file in Plink2's `.pfile` format.
- `--remove`: Indicates the file containing the list of outliers to be excluded from the dataset.
- `--out`: Specifies the path and base name for the output files.

### Output Files

- `.bed`, `.bim`, `.fam`: These files represent the Plink binary dataset after outliers have been removed.

### Script Execution

To execute the script, navigate to the directory where the script is located and run it from the command line:

```bash
bash remove_outliers.sh
```

Ensure the paths to the input `.pfile` and outlier list are correctly set within the script before execution. Adjust the `CHROMDIR` variable as necessary to match the directory structure of your specific project. This script assumes that the necessary pre-processing steps, including outlier identification, have already been completed and that the results are available in the expected directories.

---

## Script 10                                  
## 10_compute_global_pcs.R

This R script computes the global methylation PCs to use as covariates in the QTL analysis model.

### Prerequisites

This script uses multiple R packages, which need to be installed and loaded before running the script:

- `RnBeads.hg19`, `RnBeads.hg38`, and `RnBeads` for DNA methylation data analysis.
- `GenomicRanges` and `liftOver` for handling genomic data transformations.
- `isva` and `biomaRt` for statistical analysis and annotation lookup, respectively.
- `RNOmni` and `PCAtools` for principal component analysis.
- `tidyverse` for data manipulation and visualization.

### Functions

1. **get_pcs**: Estimates significant principal components for gene methylation data, applying normalization and dimensionality reduction techniques like the Gavish-Donoho and Marchenko-Pastur methods.

2. **match_mvals**: Filters and matches methylation values to a subset of samples based on clinical metadata, ensuring consistency across datasets.

### Main Function

The `main` function orchestrates the overall process:

- It iterates over a set of predefined cell types.
- For each cell type, it loads normalized methylation data, subsets this data to relevant samples, computes global principal components, and saves these components for further analysis.

### Execution Steps

1. **Directory Setup**: Establishes the working and output directories.
2. **Data Loading and Subsetting**: Loads methylation data specific to each cell type and subsets this data based on predefined sample lists.
3. **Principal Component Analysis**: Computes principal components for the methylation data to capture the major variance and saves these results.

### Usage

To run this script, you need to invoke the `main` function with appropriate indices for the array, region, and cell type as parameters. This is done within an R environment, ensuring all paths and directories are correctly set up:

```R
# Run the main function for a specific setup
main(array_idx = 1, region_idx = 1, celltype_idx = 1)
```

### Output Files

The script outputs files containing the principal components for each cell type, stored as `.rds` files in the specified output directory. These files are crucial for downstream analyses, including correlation studies and association testing with phenotypic data.

### Notes

- Ensure that the directory paths and file names are correctly set up in your environment.
- Adjust the `savedir` and `datadir` paths according to your project's directory structure.
- This script is highly dependent on the structure and quality of input data and requires thorough quality control and preprocessing steps before execution.

---

## Script 11                                  
## 11_format_qtl_input.R

This R script is crafted to prepare data for quantitative trait loci (QTL) analysis, focusing on matching covariates and phenotypic outcomes across different summary types and cell types. 

### Prerequisites

Before running this script, ensure that the following R packages are loaded:

- `tidyverse`: For data manipulation and visualization.

### Main Function

The `main` function oversees the execution flow of the script, processing through combinations of cell types, region types, and summary methods:

- **Covariate Matching**: Aligns multiple covariate sources including genetic principal components, methylation principal components, and other phenotypic variables.
- **Outcome Matching**: Aligns the matched covariates with phenotypic outcomes, which are summarized methylation values across different genomic regions.

### Execution Steps

1. **Covariate Preparation**: Assembles all relevant covariates needed for QTL analysis, adjusting for factors like ancestry, global methylation, and additional phenotype variables.
   
2. **Phenotypic Data Alignment**: Matches these covariates with summarized methylation data to prepare for regression analysis in QTL studies.

### Function Descriptions

- **match_covariates**: Combines various covariates (genetic PCs, methylation PCs, and clinical data) into a single dataset for each cell type and region.
  
- **match_covariates_outcome**: Aligns the covariate dataset with outcome variables from methylation summaries to ensure that each sample's data aligns across all inputs and outputs.

### Usage

Invoke the `main` function to start the process. Ensure that directory paths (`datadir`, `savedir`) and variable settings (like `without_proportions`) are correctly configured:

```R
# Run the script
main()
```

### Output Files

Outputs include TSV files containing matched covariate and outcome data, formatted for use in downstream statistical analysis:

- Matched covariate files: `*_matched_covariates.tsv`
- Matched phenotypic outcome files: `*_matched_pheno.bed`

These files are compressed and stored in the specified output directory.

### Notes

- Adjust the `chromdir` variable to point to specific chromosome directories if needed.
- Ensure that the paths to data files are correct and accessible.
- The script includes parameters that can toggle the inclusion of certain covariates (`without_proportions`), which should be set based on the specifics of the QTL analysis plan.
- It's essential to verify that the data used in this script has passed quality control checks to avoid issues during the QTL analysis.

This script is pivotal for ensuring that the data inputs into a QTL analysis are correctly aligned and formatted, which is crucial for the accuracy of the analysis results.

---

## Script 12                                  
## 12_map_qtls.sh

This Bash script is designed to map quantitative trait loci (QTLs) using a specific QTL analysis pipeline, integrating genotype, phenotype, and covariate data. It utilizes the TensorQTL toolkit, a powerful resource for the integration of transcriptional profiling and genome-wide association data. It calls the scripts provided in the [QTL pipeline](https://github.com/tyeulalio/QTL_pipeline) GitHub.

### Prerequisites

Ensure the following are in place before running the script:

- Conda or Mamba environment with TensorQTL installed.
- Access to the required directories and files.
- Download scripts provided in the [QTL pipeline](https://github.com/tyeulalio/QTL_pipeline) GitHub.

### Configuration

Set variables to specify the summary type, region type, cell type, and other parameters relevant to the analysis:

- `SUMMARY_TYPE`, `REGION_TYPE`, `CELL_TYPE`: Define the scope of analysis based on genomic summaries and cell-specific data.
- `PROPORTION_TYPE`: Specifies whether the analysis includes proportion-based adjustments.

### Execution Steps

1. **Environment Setup**: Activates the necessary Conda environment to use TensorQTL.
   
2. **File Paths**: Configures paths for phenotype data, covariates, and genotype data, ensuring that they are correctly addressed based on the analysis setup.

3. **Output Directory**: Establishes and prepares the directory where the analysis results will be stored.

4. **Parameter Adjustment**: Sets the cis-window parameter, which is adjustable based on the summary type.

### Usage

Execute the script by passing the summary type, region type, and cell type as arguments. Make sure the Conda environment is properly configured with TensorQTL:

```bash
bash map_qtls.sh avgs gene_body neuron
```

### Key Variables

- `RUNNAME`: Combines cell type, region type, and summary type to define the specific run.
- `PHENO_FILE`: Path to the phenotype file.
- `COVS_FILE`: Path to the covariates file.
- `PLINK_PATH`: Path to the directory containing PLINK binary files.
- `OUTPUT_DIR`: Directory for storing output files.

### Output Files

The script generates QTL mapping results in the specified output directory, which includes various statistical and graphical representations of QTLs.

### Notes

- The `CIS_WINDOW` is particularly crucial for defining the genomic range around transcription start sites (TSS) considered in the QTL mapping.
- Ensure that the file paths are correct and that required data files (phenotype, covariates, genotype) are accessible.
- The script is part of a larger pipeline, and its proper function depends on the successful execution of previous steps in the pipeline, such as phenotype and covariate preparation.

This script facilitates the complex task of QTL mapping by automating the integration and analysis of various data types, aiming to elucidate the genetic mechanisms underlying trait variability.