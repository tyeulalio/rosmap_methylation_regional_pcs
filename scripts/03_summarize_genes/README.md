# Cell type deconvolution and Batch correction

This pipeline contains a set of functions designed to deconvolve cell type-specific methylation  and perform batch correction.

## Table of Contents

- [01_summarize_genes.R](#script-01)
    - [load_methylation](#load_methylation1)
    - [lift_methylation](#lift_methylation)
    - [get_annotations](#get_annotations)
    - [get_pcs](#get_pcs)
    - [compute_averages](#compute_averages)
    - [select_genes](#select_genes)
    - [summarise_genes](#summarise_genes)
    - [combine_results](#combine_results)
    - [compute_pcs](#compute_pcs)
    - [map_cpgs](#map_cpgs)
    - [save_results](#save_results)
    - [main](#main1)
                                    

## Script 01                                    
## 01_summarize_genes.R

This R script is designed as a comprehensive pipeline for analyzing DNA methylation data across different cell types and genomic regions. It integrates a series of functions to handle data loading, transformation, annotation, summarization, and result compilation. Specifically, the script involves:

1. **Data Loading**: Functions to load methylation data (`load_methylation`) and annotations (`get_annotations`), preparing them for analysis based on cell type specificity.

2. **Data Conversion**: A `lift_methylation` function to update the methylation data from the hg19 to hg38 reference genome using provided chain files and annotations.

3. **Annotation Mapping**: The `map_cpgs` function associates CpG sites with annotated genomic regions, creating a mapping that facilitates region-specific analysis.

4. **Data Summarization**: 
   - `compute_averages` calculates average methylation values across genes for each cell type.
   - `get_pcs` derives principal components from gene-specific methylation values to capture the major variance within the data, providing a reduced representation of methylation patterns.
   - `select_genes` identifies genes to be processed in a particular batch job, aiding in the distribution of workload for parallel processing.
   - `summarise_genes` uses the previous selections to calculate either averages or principal components as specified by the `summary_type` parameter.

5. **Combining Results**: The `combine_results` function consolidates the outputs from multiple batch jobs into a single dataset, facilitating downstream analysis.

6. **Utility Functions**: 
   - `save_results` for persisting the summarized data onto the disk.
   - `compute_pcs` to calculate and save principal components for specified genes.

7. **Main Function**: Orchestrates the workflow by invoking the aforementioned functions in sequence and managing the parallel processing logistics.

8. **Execution Control**: The script includes conditions to handle batch processing and arguments passed from the command line, allowing for flexible and scalable execution on a computing cluster.

The script is designed with modularity and reusability in mind, making it suitable for a wide range of methylation datasets. It emphasizes batch processing efficiency and is intended for use in a high-throughput computing environment.


Below are descriptions and usage guidelines for each of the functions for this script.

### <a name="load_methylation1">`load_methylation`</a>

#### Description

The load_methylation function is designed for importing processed methylation data specific to a given cell type. It reads in the pre-computed and cleaned methylation values stored in an RDS file format. The function assumes that the data has already undergone batch correction and SNP filtering as part of preprocessing. It returns the methylation beta values in a structured format suitable for further analysis steps such as PCA, differential methylation analysis, or correlation studies.

#### Parameters

- `cell_type`: A character string specifying the cell type for which methylation data is to be loaded. The name provided must correspond to the suffix of the RDS file containing the methylation data for that cell type after batch correction.

#### Returns

- `mvals`: Returns a data frame or matrix with CpG sites as rows and samples as columns, containing the beta values of methylation data for further analysis.

---
    
### <a name="lift_methylation">lift_methylation</a>

#### Description

The lift_methylation function transitions methylation data from the human genome build hg19 to hg38 using the provided chain file for remapping. It starts by importing the appropriate chain file for the conversion process. The function then cross-references the provided methylation values (mvals) with the hg19 CpG site annotations to create a new data frame that only includes the sites present in the input data. After executing the liftOver function, it creates a mapping between the original and the new hg38 CpG site positions and purges any CpGs lost during the transition process. The finalized, lifted methylation data is then ordered and saved for further use.

#### Parameters
- `mvals`: A matrix or data frame containing beta-values of methylation data with row names as CpG site IDs based on the hg19 genome build.
- `cell_type`: A string indicating the cell type corresponding to the methylation data, used to tailor the output file name.

### Returns
- `ordered_mvals`: A data frame of methylation beta values with CpG site IDs updated to the hg38 genome build, ordered to align with the liftOver results.


#### Details
The lifting over process may result in a loss of some CpG sites that do not map directly from hg19 to hg38, which will be excluded from the returned data. It's important to ensure that the CpG site IDs in the input mvals are compatible with the CpG annotations used by the liftOver tool.

---
    
### <a name="get_annotations">`get_annotations`</a>

#### Description

Retrieves genomic annotations for specified regions from the RnBeads hg38 dataset and, if required, combines them with gene biotype information from Ensembl via BioMart. It is capable of handling different types of genomic regions such as enhancers, gene bodies, preTSS (regions before transcription start sites), or custom-defined annotations based on the input `region_type`. The function can merge annotations with gene biotype data, creating a comprehensive dataset for downstream analysis.

#### Parameters

- `region_type`: A string specifying the type of genomic region for which to retrieve annotations. Accepted values include 'enhancers', 'gene_body', 'preTSS', or other region types corresponding to the available annotation files.

#### Returns

- `full_annots`: A data frame of annotations that includes genomic coordinates, gene IDs, and additional information such as gene biotype, if available. This data frame is ready for downstream analysis or for mapping CpG sites to these regions.

#### Details

Depending on the `region_type` specified, the function will either load pre-saved RDS files containing enhancer annotations, gene body annotations, or perform a filtering operation to retrieve only the relevant regions. For non-standard region types, it will attempt to load a compiled annotation file named after the `region_type`. If the gene biotype information is not already available, it will use BioMart to retrieve and join this data to the annotations.

#### Examples

```r
# To retrieve annotations for enhancers:
enhancer_annotations <- get_annotations("enhancers")

# To retrieve annotations for gene bodies excluding certain regions:
gene_body_annotations <- get_annotations("gene_body")
```

#### Note

- The `datadir` variable must be defined in your working environment and should point to the directory where the annotation RDS files are located.
- `savedir` needs to be specified to indicate the output directory for saving the extended annotations RDS file.

#### See Also

- `readRDS`, `saveRDS`: to read and write R objects in RDS format.
- `useEnsembl`, `getBM`: to retrieve gene biotype information from Ensembl BioMart.
- `write_csv`: to write CSV files in R.
- `dir.create`: to create directories within the file system.

Ensure to replace `"/path/to/annotate_regions/01_annotate_sites/"` with the actual directory path that contains your annotation RDS files before running the function.

---

### <a name="get_pcs">`get_pcs`</a>

#### Description

Performs principal component analysis (PCA) on gene-level methylation values (mvals) to identify significant principal components (PCs) using both the Gavish-Donoho (GD) and Marchenko-Pastur (MP) methods. It returns the significant PCs, the corresponding loadings, estimated dimensions for each method, and the computation time.

#### Parameters

- `gene_mvals`: A matrix of methylation values for a specific gene across multiple samples.

#### Returns

A list containing the following components:

- `sig_pcs`: A data frame of significant principal components for the given gene's methylation data.
- `loadings`: The loadings corresponding to the significant PCs.
- `gv_dim`: The estimated dimension using the Gavish-Donoho method.
- `mp_dim`: The estimated dimension using the Marchenko-Pastur law.
- `error`: NA placeholder, which can be used to store error messages in future extensions of the function.
- `duration_secs`: The time taken to perform the PCA and dimension estimation, in seconds.

#### Details

The function computes the principal components of the methylation data, estimates the number of significant components using the GD and MP methods, and then filters the PCs based on these estimates. It calculates the proportion of variance explained by each significant PC and formats the results into a user-friendly structure.

#### See Also

- `PCAtools::pca`: Function from the PCAtools package used for principal component analysis.
- `chooseGavishDonoho`, `chooseMarchenkoPastur`: Functions for estimating the optimal number of principal components.
- `Sys.time`, `as.numeric`: Functions used to calculate the duration of the PCA process.

#### Notes

- It is important to ensure `gene_mvals` is a correctly formatted matrix with genes as rows and samples as columns before using this function.
- The `error` component in the return list is currently not in use but provides a template for error handling in more complex workflows.
- The duration measured reflects the computational efficiency of the PCA process and may vary based on system architecture and dataset size.

---

### <a name="compute_averages">`compute_averages`</a>

#### Description

Calculates the average methylation values for genes within specified regions for a given cell type. It checks for an existing results file and skips computation if the file is found. Otherwise, it performs averaging across all CpGs associated with each gene and saves the results in a specified directory.

#### Parameters

- `mvals`: A matrix containing methylation beta values for various CpGs across samples.
- `gene_map`: A data frame mapping genes to their corresponding CpGs.
- `region_type`: A string representing the genomic region of interest (e.g., 'promoters', 'gene_body').
- `cell_type`: A string denoting the cell type for which the methylation data is summarized.
- `cleaned_str`: A string used to prefix saved file names, indicating whether the data is cleaned or raw.

#### Returns

Integer `1` if the operation completes successfully and the results are saved. `NA` is returned if the savefile already exists, indicating that the function skipped the operation to avoid overwriting or redundant computations.

#### Details

The function iterates over a list of unique genes derived from `gene_map` and applies a mean function to compute the average beta value across all CpGs for each gene. If only a single CpG is associated with a gene, that CpG's beta value is returned without averaging.

#### See Also

- `mclapply`: Function from the parallel processing package to apply operations in parallel.
- `saveRDS`, `readRDS`: Functions to save and load R objects to and from RDS files.
- `paste0`, `file.exists`: Functions used to concatenate strings and check for file existence.

#### Notes

- It is crucial to provide the correct `savedir` path where the output RDS files will be saved.
- The parallel processing with `mclapply` is set to use 16 cores; this value should be adjusted based on the available computational resources.
- The script is intended to be part of a larger pipeline, so it assumes the presence of other scripts and data structures.
- The returned value of `1` is a placeholder to indicate successful execution; error handling and messaging can be further developed to provide more detailed feedback.

#### Files Created

- A `.rds` file named following the pattern `<cleaned_str><region_type>_<cell_type>_avg_mvals.rds`, containing a data frame of averaged methylation values.

---

### <a name="select_genes">`select_genes`</a>

#### Description

Determines the subset of genes to be processed in a batch job based on the array index and total number of jobs. It calculates the start and end indices for the subset of genes to ensure an even distribution of workload across batch jobs.

#### Parameters

- `array_idx`: An integer representing the index of the current batch job array.
- `mvals`: A matrix containing methylation beta values for various CpGs across samples.
- `gene_map`: A data frame mapping genes to their corresponding CpGs.
- `region_type`: A string representing the genomic region of interest.
- `cell_type`: A string indicating the cell type for which the methylation data is summarized.
- `cleaned_str`: A string used to differentiate between cleaned and uncleaned data sets in file naming.

#### Results

A data frame (`job_map`) containing the filtered gene-to-CpG map for the specified job array index, ready for further processing such as summarizing genes.

#### Details

This function is part of a larger batch processing script intended to distribute the load of methylation data processing. It assumes that `total_jobs` is a predefined variable accessible within the function scope.

#### See Also

- `filter`: dplyr function to subset rows based on condition.
- `ceiling`: Math function to round up values to the nearest integer.
- `print`: Basic R function to output text in the console.

#### Notes

- The function prints the progress and selected gene range for the current job, which is useful for logging and debugging purposes.
- The end index is calculated to ensure it does not exceed the total number of genes available.
- The function returns `NA` if the start index exceeds the total number of genes, indicating completion of jobs.

#### Files Created

None directly, but `job_map` can be used by subsequent functions to create output based on the selected gene subset.

---

### <a name="summarise_genes">`summarise_genes`</a>

#### Description

The `summarise_genes` function compiles gene-specific data summaries using provided methylation values (mvals) and a mapping of CpG sites to genes (job_map). It supports summarization by averaging methylation values or by computing principal components (PCs). This function is flexible and can be applied to any subset or 'job' of genes, making it suitable for parallel processing where each job may represent a set of genes to be processed.

#### Parameters

- `mvals`: Methylation value matrix with CpG sites as rows and samples as columns.
- `job_map`: A data frame mapping CpG sites to genes, specifying which genes are included in the current job's processing.
- `summary_type`: A character string specifying the type of summary to compute, either 'avgs' for averages or 'pcs' for principal components.

#### Returns

This function returns a different structure based on the `summary_type` parameter:
- If `summary_type` is 'avgs', it returns a data frame of averaged methylation values for the genes processed.
- If `summary_type` is 'pcs', it returns a list containing principal component scores, loadings, and additional PCA information for each gene.

#### Side Effects

- Prints status updates and processing times to the console.
- Potentially performs computationally intensive operations to calculate PCs.

#### See Also

- `lapply()`: Used for applying functions over lists, commonly used in the context of parallel processing.
- `PCAtools::pca()`: PCA computation used when `summary_type` is 'pcs'.
- `t()`: Transpose function, used here to structure data frames correctly before returning.

#### Notes

- For principal components analysis, dimensionality is estimated using both Gavish-Donoho and Marchenko-Pastur methods.
- It is important that the `job_map` accurately reflects the CpG-to-gene mapping for the subset of genes processed to ensure valid summaries.
- The function is designed to be run as part of a larger pipeline that includes pre-processing steps and subsequent analysis.

--- 

### <a name="combine_results">`combine_results`</a>

#### Description

The `combine_results` function orchestrates the collation of processed methylation data across different cell types and genomic regions. It compiles principal component scores, loadings, and related summary information into consolidated datasets for each specified cell type and region type. This function is particularly useful after parallel processing of large datasets where results are generated in a batch-wise manner. The consolidation is essential for subsequent analyses that require aggregated data across various partitions.

#### Parameters

This function does not take parameters directly. Instead, it operates on predefined lists of `region_types` and `cell_types` within its scope and uses the `total_jobs` variable to determine the scope of the consolidation.

#### Returns

Nothing is explicitly returned. However, this function saves multiple files to disk:
- A combined principal components (PCs) dataset for each cell and region type.
- A dataset of PC loadings for each cell and region type.
- A summary information dataset containing dimensions and other metrics for each gene analyzed in each cell and region type.

#### Side Effects

- Prints status updates to the console regarding the processing and saving of results.
- Saves combined results to the specified directory as RDS files.
- Calls `lapply` to iterate over an index range, applying the `load_pcs` and `load_avgs` functions to accumulate results from individual job files.

#### See Also

- `lapply()`, `do.call()`, `rbind()`: Functions used to apply operations over list structures and combine results into data frames.
- `readRDS()`, `saveRDS()`: Functions to read from and write to RDS files, R's native format for serializing data objects.

#### Notes

- It's critical to ensure that the directory structure and file naming conventions are consistent with those expected by the script.
- Ensure that the `savedir` path is correctly set to the location where the individual job results are stored before running the function.
- The function's performance and efficiency are contingent on the organization and size of the processed files and the number of jobs to combine.

---

### <a name="compute_pcs">`compute_pcs`</a>

#### Description

The `compute_pcs` function performs Principal Component Analysis (PCA) on methylation data for specified genes within a given cell type and region. It is part of a larger pipeline that processes methylation data, where PCA is used for dimensionality reduction to identify the major patterns of variation within the data. The function calculates significant principal components and their corresponding loadings, along with the explained variance. It then saves the significant components and related data as a summarized result for subsequent analysis.

#### Parameters

- `mvals`: A matrix of methylation values for further PCA analysis.
- `gene_map`: A mapping of CpG sites to gene regions for the subset of data being analyzed.
- `region_type`: The type of genomic region being summarized (e.g., 'promoters', 'exons').
- `array_idx`: An index used to identify and process a subset of genes in a larger batch processing context.
- `cell_type`: The cell type for which the PCA is being computed.
- `cleaned_str`: A string indicating the type of data cleaning applied, used in naming the output file.

#### Returns

The function does not return a value but saves PCA results to a file.

#### Side Effects

- Prints progress messages to the console, indicating the job index and the range of genes being processed.
- Generates an RDS file containing the PCA results for the specified subset of genes.
- Creates necessary directories for storing the results, if they do not already exist.

#### See Also

- `pca`: Function from the `PCAtools` package for performing PCA analysis.
- `get_pcs`: A function to determine the number of significant principal components using the Gavish-Donoho method.
- `readRDS`, `saveRDS`: Functions to read from and write to RDS files.

#### Notes

- Ensure that the `datadir` and `savedir` paths are correctly set to the location where the input data is stored and where the output results are to be saved, respectively.
- The function's efficiency in processing and saving data can be affected by the size of the methylation matrix and the computing resources available.
- The actual computation of PCA is contingent upon successful data preparation steps executed prior to this function.

---

### <a name="map_cpgs">`map_cpgs`</a>

#### Description

The `map_cpgs` function associates CpG sites with their corresponding gene regions based on genomic annotations. It utilizes the GenomicRanges data structure to find overlaps between CpG sites and annotated genomic regions, creating a mapping that is used in downstream analysis for summarizing data by these regions. The function handles the mapping for one specific cell type and genomic region at a time and saves the mapping results to an RDS file for persistence.

#### Parameters

- `lifted_mvals`: A matrix of methylation values where row names are CpG site identifiers lifted to the hg38 assembly.
- `annots`: A data frame or list containing genomic annotations with which to align the CpG sites.
- `region_type`: A string representing the type of genomic region being mapped, such as 'enhancers' or 'gene_body'.
- `cell_type`: A string specifying the cell type for which the mapping is being performed.

#### Returns

The function returns a gene map as a data frame with 'gene_id' and 'cpg_id' columns indicating the association between genes and CpG sites.

#### Side Effects

- Reads from or writes to RDS files located at paths constructed using `savedir`, `region_type`, and `cell_type`.
- Creates a directory for saving annotated regions if it doesn't already exist.
- May print to the console the contents of loaded or saved data frames for verification purposes.

#### See Also

- `makeGRangesFromDataFrame`, `findOverlaps`: Functions from the `GenomicRanges` package to manipulate genomic intervals and find overlaps.
- `readRDS`, `saveRDS`: Functions to read from and write to RDS files.

#### Notes

- It is important to ensure that `savedir` is correctly set to the desired output directory before running the function.
- The function checks if a saved gene map already exists before computing the mapping to prevent redundant calculations.
- The `remove=FALSE` parameter in the `separate` function call ensures the original CpG identifiers are retained for reference.
- The mapping takes into account the possibility of multiple CpG sites overlapping with a single gene, which is common in genomic datasets.

---

### <a name="save_results">`save_results`</a>

#### Description

The `save_results` function is responsible for writing the summarized gene data to disk. It creates a directory structure based on the type of summary, the genomic region, and the cell type involved, then saves the summary information in an RDS file for each batch job indexed by `array_idx`. This modular approach facilitates the organization of results from potentially large and computationally intensive analyses.

#### Parameters

- `gene_summaries`: A data frame or list containing the summary statistics for genes. It can be the output of a function that computes gene averages or principal components.
- `summary_type`: A string denoting the type of summary performed, such as 'avgs' for averages or 'pcs' for principal components.
- `region_type`: A string indicating the genomic region associated with the summary, e.g., 'promoters' or 'enhancers'.
- `cell_type`: The cell type for which the summary was computed.
- `array_idx`: An index used to differentiate between summary files from different batch jobs.

#### Returns

The function returns `1` as a placeholder for successful execution and side-effect-driven nature (saving files).

#### Side Effects

- If the directory specified by `summarydir` does not exist, it is created.
- Saves the summary results to an RDS file at a specified path, potentially overwriting existing files without warning.
- Prints a message to the console indicating the action being taken.

#### See Also

- `dir.create`, `saveRDS`: Functions used to create directories and save RDS files.

#### Notes

- The function is designed to be used within a batch processing system where multiple jobs may be run in parallel, necessitating a unique identifier (`array_idx`) for each job's output.
- It's important to check the logic for overwriting files and create necessary backups or checks to prevent data loss.
- The print statement at the beginning helps users track the progress and results of their jobs when running the pipeline.

---

### <a name="main">`main`</a>

#### Description

The `main` function orchestrates the overall process of methylation data analysis for different cell types across specified genomic regions. It performs a sequence of steps, including loading methylation data, lifting CpG positions from the hg19 to hg38 assembly, obtaining relevant annotations, mapping CpGs to those annotations, summarizing genes, and saving the results. The function is designed to work in a batch processing environment, allowing for parallel processing of data across an array of cell types and genomic regions.

#### Parameters

- `array_idx`: The index of the current batch job, used for processing a subset of the data.
- `region_idx`: The index that specifies which genomic region type is to be processed.
- `celltype_idx`: The index corresponding to the cell type to be processed.
- `total_jobs`: The total number of jobs that will be run, used for dividing the workload.

#### Returns

The function is used for its side effects and does not return any object; it outputs to the console and saves files to the disk.

#### Side Effects

- Calls several other functions that read, process, and write data, producing outputs in the form of RDS files.
- Generates console output to track the progress of data processing.
- Alters the state by saving the processed data to files and potentially overwriting existing ones.

#### Notes

- The `main` function is the entry point of the script and governs the flow of execution.
- It assumes the existence of global variables such as `savedir`, `precleaned_str`, and `cleaned_str` which should be defined before calling the function.
- It uses `commandArgs` to allow passing arguments from the command line, making it suitable for batch processing.
- The function ends by calling `quit()` to exit the R session, which is typical behavior for a script designed to run batch jobs on a cluster.

#### See Also

- `readRDS`, `writeRDS`, `lapply`, `print`, `Sys.time`: Other functions and constructs used within `main`.
