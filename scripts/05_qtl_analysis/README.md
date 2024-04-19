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