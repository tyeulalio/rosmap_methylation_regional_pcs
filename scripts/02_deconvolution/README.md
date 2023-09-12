# Cell type deconvolution and Batch correction

This pipeline contains a set of functions designed to deconvolve cell type-specific methylation  and perform batch correction.

## Table of Contents

- [01_deconvolve_data.R](#script-01)
    - [load_methylation](#load_methylation1)
    - [load_clinical](#load_clinical1)
    - [format_data](#format_data1)
    - [get_cell_proportions](#get_cell_proportions)
    - [deconvolve_data](#deconvolve_data)
    - [main](#main1)
- [02_batch_correction.R](#script-02)
    - [load_methylation](#load_methylation2)
    - [correlate_traits](#correlate_traits)
    - [load_clinical](#load_clinical2)
    - [format_data](#format_data2)
    - [process_betas](#process_betas)
    - [run_pca](#run_pca)
    - [correlate_pcs](#correlate_pcs)
    - [limma_remove_batch](#limma_remove_batch)
    - [main](#main2)
- [03_celltype_umap.R](#script-03)
    - [load_meth](#load_methylation3)
    - [format_meth](#format_meth)
    - [run_pca](#run_pca3)
    - [plot_umap](#plot_umap)
    - [cluster_data](#cluster_data)
    - [main](#main3)
                                    

## Script 01                                    
## 01_deconvolve_data.R

This R script serves as an end-to-end pipeline for DNA methylation analysis, with a particular focus on data deconvolution and cell-type proportion estimates. After loading and formatting raw methylation and clinical data, the script utilizes EpiSCORE for estimating the proportions of various cell types present in each sample. Subsequently, it applies Tensor Composition Analysis (TCA) to deconvolve the methylation data, thereby isolating cell-type-specific methylation patterns. The output includes formatted clinical data, estimated cell-type proportions, and deconvolved methylation datasets, all saved as RDS files for further analysis.

Below are descriptions and usage guidelines for each of the functions for this script.

### <a name="load_methylation1">`load_methylation`</a>

#### Description

This function is responsible for loading methylation data based on the given threshold (`thresh`). The function loads  raw data after basic preprocessing and SNP filtering.

#### Parameters

- `thresh`: A threshold value for SNP filtering in the methylation data.

#### Returns

- `betas`: A data frame or matrix containing the beta-values of methylation data.

---
    
### <a name="load_clinical1">`load_clinical`</a>

#### Description

This function is responsible for loading clinical and phenotype data, mapping column names to specimen IDs, and then connecting these mapped specimen IDs to individual IDs. The function achieves the following:
    
1. Reads the phenotype data from an RDS file.
2. Reads a metadata manifest (in TSV format) which contains the mapping from file formats to specimen IDs.
3. Creates a unique mapping from column names to specimen IDs by filtering and transforming the metadata manifest.
4. Maps specimen IDs to individual IDs based on the loaded phenotype data.
5. Joins the mapped specimen IDs with the phenotype data to create a fully mapped clinical dataset.
6. Saves the mapped clinical data as an RDS file for future use.

#### Parameters

None.

#### Returns

- `mapped_clinical`: A data frame containing the clinical data, including individual and specimen IDs. This information is connected to methylation column names.

#### Side Effects

- Reads an RDS file containing phenotype data and a TSV file containing metadata manifest.
- Saves the mapped clinical data into an RDS file specified by `savefile`.
- Updates the global state by loading `mapped_clinical` into the workspace.

#### Files Used

- `"../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds"`: File containing the phenotype data.
- `"../../data/ROSMAP_data/Epigenetics/Epigenetics (DNA methylation array)/SYNAPSE_METADATA_MANIFEST.tsv"`: File containing the metadata manifest for mapping.

#### Files Created

- An RDS file specified by `savefile` containing the mapped clinical data.

---
    
### <a name="format_data1">`format_data`</a>

#### Description

This function formats the methylation beta-values and clinical data to ensure that they are compatible for downstream analysis. The function performs the following operations:

1. Identifies the common column names between the beta-values and clinical data.
2. Orders both datasets based on these common names.
3. Checks that the ordering of both datasets matches; if not, the function will halt with an error message.
4. Sets row and column names appropriately and formats categorical and numerical variables.
5. Saves the formatted beta-values as an RDS file.
6. Returns a list containing both formatted beta-values and clinical data.

#### Parameters

- `betas`: A data frame or matrix containing methylation beta-values. Columns should represent samples and rows should represent genomic features (CpG sites).
- `mapped_clinical`: A data frame containing clinical data with mapping to individual and specimen IDs.

#### Returns

- `formatted_dat`: A list containing two elements:
    - `formatted_betas`: A data frame or matrix of formatted beta-values.
    - `formatted_clinical`: A data frame of formatted clinical data.

#### Side Effects

- Reads global variable `savedir` to determine where to save the RDS file.
- Saves the formatted beta-values into an RDS file specified by `savefile`.

#### Files Used

None, the function relies on data passed as parameters.

#### Files Created

- An RDS file specified by `savefile` containing the formatted beta-values.

---

### `get_cell_proportions`

#### Description

This function estimates the cell-type proportions in DNA methylation data using the EpiSCORE algorithm. The function performs various operations, such as:

1. Loading noob-preprocessed methylation data and matching it with the clinical dataset.
2. Loading the brain reference matrix from EpiSCORE.
3. Filtering methylation data based on variance thresholds.
4. Calculating average methylation beta-values in specific genomic regions (200bp upstream of TSS).
5. Running the EpiSCORE algorithm to get cell-type proportions.

#### Parameters

- `formatted_betas`: A data frame or matrix containing methylation beta-values, already formatted.
- `formatted_clinical`: A data frame containing clinical data, already formatted.
- `thresh`: Threshold for SNP filtering applied in previous steps.

#### Returns

- `celltype_proportions$estF`: A data frame containing estimated cell-type proportions for each sample.

#### Side Effects

- Reads global variables `savedir` and `precleaned_str` to determine where to save RDS files.
- Saves estimated cell-type proportions and reference panel into separate RDS files.

#### Files Used

- Noob-preprocessed methylation data from a specified file path.
- Brain reference matrix from EpiSCORE (loaded via `data(BrainRef)`).

#### Files Created

- An RDS file containing the estimated cell-type proportions (`snpThresh_episcore_estimated_celltypes.rds`).
- An RDS file containing the reference panel used (`snpThresh_episcore_reference_panel.rds`).

#### Functions Called

- `load_noob()`: Internal function to load noob-preprocessed methylation data.
- `constAvBetaTSS()`: Function to compute average methylation beta-values upstream of TSS.
- `wRPC()`: EpiSCORE function to estimate cell-type proportions.

#### Notes

Make sure the reference matrix (`BrainRef`) is loaded into the environment before running this function. This is typically from the EpiSCORE package and can be loaded using `data(BrainRef)`.

---

### `deconvolve_data`

#### Description

This function performs data deconvolution using Tensor Composition Analysis (TCA) on DNA methylation data. The function specifically aims to:

1. Summarize and modify cell-type proportions.
2. Create a design matrix for covariates that do not directly affect methylation at the cell-type level.
3. Select the relevant subsets of betas and cell-type proportions based on individuals in the clinical dataset.
4. Run TCA to deconvolve the DNA methylation matrix into cell-type-specific components.
5. Save the TCA outputs and generate a cell-type-specific tensor.

#### Parameters

- `formatted_betas`: A data frame or matrix of the methylation beta-values, formatted.
- `proportions_df`: A data frame containing cell-type proportions.
- `formatted_clinical`: A data frame containing the formatted clinical data.

#### Returns

The function doesn't explicitly return anything but saves multiple RDS files containing:

1. Formatted cell-type proportions.
2. TCA outputs.
3. Cell-type-specific tensor.

#### Side Effects

- Reads from and writes to global variables `savedir` and `precleaned_str` to determine the directory for saving output files.
  
#### Files Used

None explicitly, although the function is dependent on the input parameters, which are assumed to be pre-processed and loaded beforehand.

#### Files Created

1. An RDS file containing the formatted cell-type proportions (`formatted_celltype_proportions.rds`).
2. An RDS file containing the TCA output (`tca_out.rds`).
3. An RDS file containing the cell-type-specific tensor (`tca_tensor.rds`).

#### Functions Called

- `tca()`: Function to run Tensor Composition Analysis.
- `tensor()`: Function to convert the betas matrix and TCA output into a cell-type-specific tensor.

#### Notes

1. The function assumes that the EpiSCORE-derived cell-type proportions are already available.
2. The function makes specific decisions about combining or omitting cell types based on the summary statistics of cell-type proportions. Therefore, this function is specific to the data and may need adjustments for different datasets.
3. The function assumes that the TCA package or equivalent is loaded and available for use.
4. Batch correction and missing data imputation for PMI (Post Mortem Interval) are performed within the function.
5. Parallel processing is commented out (`#parallel = TRUE, num_cores = 16`), but can be enabled for more efficient computations.

---

### <a name="main1">`main`</a>

#### Description

This function serves as the primary workflow for the deconvolution pipeline. It orchestrates the following steps:

1. Loading the methylation data.
2. Loading clinical phenotype data.
3. Formatting both methylation and clinical data.
4. Saving and re-loading the formatted clinical data.
5. Estimating cell-type proportions.
6. Performing data deconvolution using the estimated cell-type proportions and the methylation data.

#### Parameters

The function doesn't take any parameters; instead, it assumes that necessary functions and variables are loaded into the environment.

#### Returns

The function doesn't explicitly return anything but orchestrates the entire workflow of loading, formatting, estimating cell-type proportions, and deconvolving the methylation data.

#### Side Effects

- Reads from and writes to global variables `savedir` and possibly others (not explicitly mentioned) to determine the directory for saving output files.
- Calls various other functions (`load_methylation`, `load_clinical`, `format_data`, `get_cell_proportions`, and `deconvolve_data`) that read and write files.
  
#### Files Used

1. Reads a pre-saved RDS file named `formatted_clinical.rds` from a directory specified by `savedir`.

#### Files Created

1. Saves an RDS file containing the formatted clinical data (`formatted_clinical.rds`).

#### Functions Called

- `load_methylation(thresh)`: Loads methylation data with a specific threshold.
- `load_clinical()`: Loads clinical phenotype data.
- `format_data(betas, mapped_clinical)`: Formats methylation and clinical data.
- `get_cell_proportions(formatted_betas, formatted_clinical, thresh)`: Estimates cell-type proportions based on methylation data.
- `deconvolve_data(formatted_betas, proportions_df, formatted_clinical)`: Performs data deconvolution using Tensor Composition Analysis (TCA).

#### Notes

1. Assumes that all helper functions and necessary libraries are loaded and available for use.
2. The function makes use of the system time to print time stamps, which can be helpful for debugging and performance tracking.

---

## Script 02
## 02_batch_correction.R

This script performs batch correction on the methylation data. It correlates phenotypic traits with the top methylation PCs to find traits that have high correlation with the data.

### <a name="load_methylation2">`load_methylation`</a>

#### Description

This function is responsible for loading methylation data specific to a given cell type. The function can either load preprocessed 'bulk' methylation data or deconvolved methylation data for specific cell types as obtained through Tensor Composition Analysis (TCA).

#### Parameters

- `cell_type`: Specifies the type of cell for which the methylation data should be loaded. This can be 'bulk' for preprocessed methylation data or any other specific cell type for deconvolved data.

#### Returns

- `betas`: A data frame or matrix containing the beta-values of methylation data either for the 'bulk' sample or for a specific cell type as specified.

#### Notes

1. The function reads from global variables like `datadir` to determine the directory for data files.
2. Assumes that necessary directories and files exist for reading.
3. If you're loading cell-type specific data, ensure that the cell type specified is available in the TCA output.
  
---

### <a name="correlate_traits2">`correlate_traits`</a>

#### Description

This function is designed to analyze the correlation between different phenotypic traits in a given clinical dataset. The function not only calculates pairwise correlations between various traits but also performs statistical tests to determine the significance of these correlations. The correlation methods vary depending on the nature of the traits (numeric or categorical). While the main script does not call this function, it can be run interactively to examine how various traits are related.

#### Parameters

- `clinical_cts`: A data frame containing clinical and phenotypic information for each individual in the study.

#### Returns

- `formatted_corrs`: A data frame containing the pairwise correlations between the traits, p-values, and Bonferroni-adjusted p-values.

#### Side Effects

- Reads from the global variable `savedir` to determine the directory for saving output files.
- Creates a heatmap visualization of trait correlations and saves it in the specified directory.

#### Files Used

No external files are read in this function.

#### Files Created

1. Saves an RDS file containing pairwise trait correlations (`trait_trait_correlations.rds`).
2. Saves a PNG file of a heatmap visualizing the correlations between traits (`trait_correlation_heatmap.png`).

#### Functions Called

- `get_correlation(row)`: Internal helper function that calculates the correlation between a pair of traits based on their data type.

#### Notes

1. Assumes that the `tidyverse` and `ggplot2` packages are loaded and available for use.
2. The function also prints intermediate data frames and variable types for debugging purposes.
3. The function makes use of various statistical tests like Spearman's rank correlation, Chi-squared test, and Kruskal-Wallis test based on the data types of the traits being compared.

---

### <a name="load_clinical2">`load_clinical`</a>

#### Description

The `load_clinical` function serves to read in and preprocess clinical data, enhancing it with cell type proportion estimates. The function primarily reads two RDS files: one for clinical metadata and another for Tensor Composition Analysis (TCA) output. It then joins these datasets and creates new clinical groups based on predefined conditions for Alzheimer's disease (AD), control, and other categories. Additionally, the function provides the option to run trait correlation analyses through a call to `correlate_traits`.

#### Parameters

This function does not take any parameters; it reads files based on hardcoded paths.

#### Returns

- `clinical_cts`: A data frame containing merged clinical data and cell-type proportions. It also includes newly created clinical categories like "AD," "Control," and "Other."

#### Side Effects

- Reads from the global variable `datadir` to determine the directory for fetching TCA output.
- Modifies the global `clinical_cts` data frame with new clinical categories and cell type proportions.

#### Files Used

1. Reads a pre-saved RDS file containing clinical metadata, the path for which is hardcoded (`../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds`).
2. Reads another pre-saved RDS file containing TCA output (`tca_out.rds`), the path for which is determined by `datadir`.

#### Files Created

No new files are created.

#### Functions Called

- Optionally calls `correlate_traits(clinical_cts)`: To perform trait-trait correlation analyses if needed.

#### Notes

1. Assumes that all necessary R packages (`tidyverse` etc.) are loaded and available for use.
2. The function prints intermediate data frames and variable names, useful for debugging and ensuring data integrity.
3. New clinical categories ("AD," "Control," "Other") are created based on predefined conditions relating to cognitive diagnosis, Braak staging, and CERAD score.
  
---

### <a name="format_data2">`format_data`</a>

#### Description

The `format_data` function takes beta values, clinical data, and cell type information as input and returns a formatted version of this data. It performs several key preprocessing steps such as:

- Ensuring that both the clinical and beta data have the same ordering of individuals.
- Converting specific variables to either factors or numeric types.
- Creating new categorical variables based on existing ones (e.g., `cat_braaksc` based on Braak score).

#### Parameters

- `betas`: A matrix or data frame containing the beta-values of methylation data.
- `mapped_clinical`: A data frame containing mapped clinical data.
- `cell_type`: The cell type information. (Note: The function code currently does not use this parameter, so its purpose is unclear.)

#### Returns

- `formatted_dat`: A list containing two elements:
  - `formatted_betas`: A data frame or matrix with the beta-values, reordered to match `mapped_clinical`.
  - `formatted_clinical`: A data frame containing the clinical data, reordered and formatted.

#### Error Checks

- The function will stop execution and produce an error message if the ordering of columns in `formatted_betas` and `individualID` in `formatted_clinical` do not match.

#### Side Effects

- Performs in-place modifications to convert specific columns in the clinical data to either factor or numeric types.

#### Functions Called

No external functions are called.

#### Notes

1. The function assumes that the necessary R packages (`dplyr`, etc.) are loaded and available for use.
2. Intermediate data frames and variable names are printed for debugging and data integrity purposes.
3. The function performs error checking to ensure that the ordering of individuals is consistent between the beta and clinical data.
4. Factor and continuous covariates are specifically listed in `factor_covariates` and `cont_covariates` variables, and the function performs type conversions based on these.
  
---

### `process_betas`

#### Description

The `process_betas` function takes formatted beta values and cell type information to process methylation data. It performs several key steps:

- Removal of low-variance CpG sites.
- Conversion of methylation values based on a specified method: either Inverse Normal Transformation (INT) or M-values (mvals).
- Scaling and centering of CpG sites.

#### Parameters

- `formatted_betas`: A matrix or data frame containing the formatted beta-values of methylation data.
- `cell_type`: The cell type information. (Note: The function code currently does not use this parameter, so its purpose is unclear.)

#### Returns

Depending on the method chosen for methylation conversion (`meth_conversion`), it returns:

- `int_meth`: A data frame containing the Inverse Normal Transformed beta-values.
- `mvals`: A data frame containing the M-values of the beta-values.

#### Error Checks

- Checks and corrects the range of beta values to be between 0 and 1.

#### Side Effects

- The function prints intermediate diagnostic information, such as dimensions and head of the data frames involved, and the method of methylation conversion chosen.

#### Functions Called

- `RankNorm`: Presumably a function for performing rank-based inverse normal transformation (not explicitly defined in the code).
- `beta2m`: Presumably a function for converting beta values to M-values (not explicitly defined in the code).

#### Notes

1. The function assumes that the necessary R packages (`dplyr`, etc.) are loaded and available for use.
2. `meth_conversion` is a global variable that determines the method of methylation value conversion. It's not passed as a parameter, so it's assumed to be defined elsewhere in the code.
3. It's worth noting that the number of CpGs lost to low variance is calculated but not used further.
4. For `mvals` conversion, it scales the beta values to be between 0 and 1 before the transformation.

#### Important Issues

1. The code references a global variable `meth_conversion` that determines the method for methylation conversion but is not passed as a parameter.
2. Functions `RankNorm` and `beta2m` are assumed to be loaded or defined elsewhere in the code; their exact functionalities are not provided.

---

### `run_pca`

#### Description

The `run_pca` function performs Principal Component Analysis (PCA) on methylation values (or other types of input data, as defined by the parameter `datatype`). The function also estimates the dimensionality of the input data using Gavish-Donoho and Marchenko-Pastur methods. The significant Principal Components (PCs) are then saved to an RDS file.

#### Parameters

- `mvals`: A matrix or data frame containing the methylation values (or other data) to perform PCA on.
- `cell_type`: Cell type information. (Note: This parameter is used in naming the saved RDS file, but is not used in the computations.)
- `cv_savedir`: Directory path where the PCA results will be saved.
- `datatype`: Optional, a string specifying the type of data being processed. Used in naming the saved RDS file.

#### Returns

The function doesn't return anything but saves the significant PCs and their percent variance explained to an RDS file.

#### Error Checks

None explicitly performed within the function.

#### Side Effects

- The function prints various diagnostic and status messages to the console.
- An RDS file is saved to disk, containing the significant PCs.

#### Functions Called

- `PCAtools::pca`: An R function for performing PCA.
- `chooseGavishDonoho`: Function to estimate dimensionality using the Gavish-Donoho method.
- `chooseMarchenkoPastur`: Function to estimate dimensionality using the Marchenko-Pastur method.

#### Notes

1. The function assumes that the required R packages (`PCAtools`, etc.) are loaded and available for use.
2. The dimensionality of the data is estimated using two different methods: Gavish-Donoho and Marchenko-Pastur. This is done to potentially cross-verify the results.

#### Important Issues

1. Although `cell_type` is passed as a parameter, it is only used for naming the saved file and not for any computations within the function.
2. Functions `chooseGavishDonoho` and `chooseMarchenkoPastur` are presumed to be defined or loaded elsewhere.
3. There are no explicit error checks in the code. For instance, it doesn't check if the save directory (`cv_savedir`) exists or is writable.
4. The function performs PCA and dimensionality selection but does not return these results to the user except as saved files. This may limit its reusability in a pipeline that requires these results for subsequent steps.

---

### `correlate_pcs`

#### Description

The `correlate_pcs` function correlates Principal Components (PCs) with various traits from clinical data. Traits could be both categorical and numerical. The function performs different statistical tests based on the type of the trait. The results are then saved as a .RDS file.

#### Parameters

- `formatted_clinical`: A data frame containing formatted clinical data.
- `cell_type`: Cell type information, used for file naming.
- `cv_savedir`: The directory where the output will be saved.
- `datatype`: Optional, a string specifying the type of data being processed, used for file naming.

#### Returns

The function does not return anything but saves a .RDS file containing the p-values from the correlations with each trait.

#### Error Checks

No explicit error checks are performed within the function.

#### Side Effects

- Prints various diagnostic and status messages to the console.
- Saves a .RDS file containing p-values of the correlations.

#### Functions Called

- `readRDS`: Reads an RDS file.
- `kruskal.test`: Kruskal-Wallis test for comparing multiple groups (used for categorical traits).
- `cor.test`: Spearman correlation test (used for numerical traits).
  
#### Notes

1. The function assumes that the required R packages and variables (`precleaned_str`, `clean_global`, `cleaned_str`, etc.) are loaded or defined elsewhere.
2. Assumes that the PCs have been previously computed and saved as an RDS file.
  
#### Important Issues

1. The function does not check if the directory `cv_savedir` exists or if it's writable.
2. No explicit error handling is provided, e.g., what should happen if the RDS file does not exist or is corrupt.
3. The function assumes that the traits to be tested are hard-coded into the function, limiting its reusability.
4. While the function performs corrections for multiple comparisons using Bonferroni, this is hard-coded and not configurable by the user.
5. The function performs different statistical tests depending on the trait type but does not verify whether the chosen test is appropriate for the given data.
6. The variables `clean_global`, `precleaned_str`, and `cleaned_str` appear to be global but are not passed as arguments, which could create dependency issues.

---

### `limma_remove_batch`

#### Description

The `limma_remove_batch` function removes batch effects from methylation values (m-values) using the `removeBatchEffect` function from the limma package. Additionally, it handles covariates that are included in the model to remove the batch effect. The cleaned m-values are then saved as an RDS file.

#### Parameters

- `formatted_mvals`: A matrix/data frame containing formatted m-values.
- `formatted_clinical`: A data frame containing clinical data.
- `cell_type`: The cell type for which the analysis is being done.
- `cv_savedir`: Directory where the cleaned m-values will be saved.

#### Returns

Returns the cleaned m-values and saves them as a .RDS file.

#### Error Checks

- `stopifnot(identical(rownames(ordered_batch), colnames(formatted_mvals)))`: Checks if the row names of `ordered_batch` match the column names of `formatted_mvals`.

#### Side Effects

- Prints status updates and diagnostic information to the console.
- Saves the cleaned m-values to an RDS file.

#### Functions Called

- `readRDS`: Reads an RDS file.
- `removeBatchEffect`: Removes batch effects from the data.

#### Notes

1. Assumes that required R packages and global variables (`cleaned_str`, `savedir`, etc.) are loaded or defined elsewhere.
2. Assumes that global PCs have been previously computed and saved as an RDS file.

#### Important Issues

1. The function uses `readRDS` to read an RDS file but does not perform any error handling in case the file is not found or corrupt.
2. Variables to be removed (`remove_vars`) are hard-coded in the function, limiting its flexibility.
3. Uses `mean(pmi, na.rm=TRUE)` to fill NA values in `pmi` without notifying the user.
4. The directory `cv_savedir` is used for saving the output, but there's no check to see if the directory exists or is writable.
5. The function mutates the `formatted_clinical` data frame within the function but does not return it, potentially leading to side effects.
6. There's no check to ensure that `formatted_mvals` and `formatted_clinical` have the same samples/individuals or that they're in the same order before starting the batch correction.

Overall, while the function accomplishes its goal of batch correction, there are several hard-coded elements and assumptions that limit its flexibility and robustness.

---

### <a name="main2">`main`</a>

#### Description

The `main` function serves as a driver to perform a series of tasks, including loading methylation data, formatting clinical data, running Principal Component Analysis (PCA), and optionally removing batch effects. It is designed to process multiple cell types in a single run.

#### Parameters

- None

#### Returns

- Returns 1 upon successful completion of all steps for each cell type.

#### Error Checks

- The function itself doesn't include error checks; however, it calls other functions that should include error checks.

#### Side Effects

- Prints status updates and diagnostic information to the console.
- Creates a new directory for saving the output.
- Saves processed methylation data and clinical data to RDS files.

#### Functions Called

- `load_methylation`: Loads the methylation data.
- `load_clinical`: Loads the clinical phenotype data.
- `format_data`: Formats methylation and clinical data.
- `process_betas`: Converts beta values to m-values or other formats.
- `saveRDS`: Saves an R object to an RDS file.
- `run_pca`: Runs Principal Component Analysis on the m-values.
- `correlate_pcs`: Correlates principal components with traits.
- `limma_remove_batch`: Optionally removes batch effects from the data using the limma package.

#### Notes

1. Assumes that required R packages and global variables (`savedir`, `precleaned_str`, `cleaning`, etc.) are loaded or defined elsewhere.
2. Processes data for multiple cell types defined in `cell_types`.
3. Conditional `cleaning` flag controls whether or not batch effects are removed from the data.

#### Important Issues

1. Global variables (`savedir`, `precleaned_str`, etc.) are used without initial definitions within the function.
2. Assumes that all needed functions (`load_methylation`, `load_clinical`, etc.) are loaded or defined elsewhere but doesn't check for their existence.
3. It uses `dir.create(cv_savedir)` to create a directory but doesn't check if the directory already exists or if it was successfully created.
4. The function doesn't handle errors or issues that might occur in the functions it calls (`load_methylation`, `format_data`, etc.).
5. Variable `datatype` is declared but initially set to an empty string, which might be confusing.
6. The function uses `lapply` to iterate through `cell_types` but also calls `process_celltype(cell_type)` afterward, which might result in redundant processing for the last cell type.
7. Lack of return value or summary of what was accomplished makes it hard to programmatically check if all steps were successful.

Overall, the function provides a comprehensive pipeline for the tasks it aims to complete but lacks in robustness and error handling. It would be better to make the function more modular and include appropriate error checks and validation steps.

---

## Script 03
## 03_celltype_umap.R

This script loads methylation data for specified cell types, optionally incorporates schizophrenia-related data, and then runs Principal Component Analysis (PCA) to identify significant PCs. The PCA results are then visualized using Uniform Manifold Approximation and Projection (UMAP). Various helper functions like `load_meth`, `format_meth`, `run_pca`, and `plot_umap` modularize the tasks, and the `main` function serves as the orchestrator of these operations. This script also performs k-means clustering for cell type deconvolution validation.


### <a name="load_methylation3">`load_meth`</a>

#### Description

This function serves to load DNA methylation data for a specified cell type from a pre-defined directory. The function performs the following key actions:

1. Constructs the data file name based on the cell type specified.
2. Reads the methylation data stored in RDS format.
3. Optionally, downsamples the dataset for testing or debugging purposes.

#### Parameters

- `cell_type`: Specifies the type of cell for which methylation data is to be loaded.
- `downsample`: A Boolean flag that allows for testing on smaller datasets by downsampling.

#### Returns

The function returns a data frame or matrix containing the methylation values. This can be a downsampled subset based on the `downsample` parameter.

#### Side Effects

- Reads from the global variable `datadir` to determine the directory for loading methylation data.
  
#### Files Used

The function reads an RDS file containing the DNA methylation data, named in the format `[cell_type]_formatted_meth.rds`.

#### Files Created

No additional files are created by this function.

#### Functions Called

None explicitly, although the function heavily relies on the `readRDS` function to read methylation data from RDS files.

#### Note

Ensure that the directory path (`datadir`) from where the methylation data will be loaded is correctly specified in the global environment before running this function.

---

### `format_meth`

#### Description

The `format_meth` function is designed to format and combine methylation data across different brain regions for a given cell type. The function performs the following key steps:

1. Identifies the common CpGs (Cytosine-Phosphate-Guanine sites) that are present across all the brain regions.
2. Filters the methylation data in each brain region based on these common CpGs.
3. Combines the filtered methylation data from all brain regions into a single data matrix.

#### Parameters

- `ct_meth`: A list of data frames or matrices, each containing the methylation data for a specific brain region for the given cell type.

#### Returns

The function returns a combined data frame or matrix containing methylation data of common CpGs across all brain regions for a given cell type.

#### Side Effects

None explicitly mentioned.

#### Files Used

None explicitly, although the function relies on the `ct_meth` parameter which should be pre-loaded.

#### Files Created

No additional files are created by this function.

#### Functions Called

- `Reduce()`: Used to identify the intersection of common CpGs across all brain regions.
- `lapply()`: Applied to filter the methylation data in each brain region based on common CpGs.
- `do.call()`: Used to combine the filtered methylation data from all brain regions.

#### Note

Ensure that the list `ct_meth` passed to the function contains methylation data matrices that have row names corresponding to CpGs. The function will use these to identify common CpGs across all included brain regions.

---

### <a name="run_pca3">`run_pca`</a>

#### Description

The `run_pca` function performs Principal Component Analysis (PCA) on formatted methylation data and optionally on cleaned methylation data (`cleaned_mvals`). The key steps include:

1. Estimating the dimensions and identifying significant principal components (PCs).
2. Removing CpGs with zero variance from the methylation matrix.
3. Applying PCA to the filtered matrix to get principal components.
4. Employing dimensionality estimation methods, Gavish-Donoho and Marchenko-Pastur, to identify significant dimensions.
5. Subsetting PCs based on estimated dimensions.
6. Calculating and appending the percentage of variance explained by each principal component.
7. Saving the formatted PCs to a specified RDS file.

#### Parameters

- `formatted_meth`: A data frame or matrix containing the formatted methylation data.
- `cleaned_mvals`: An optional data frame or matrix containing cleaned methylation data (if applicable).

#### Returns

The function doesn't explicitly return anything but saves an RDS file containing formatted PCs.

#### Side Effects

The function reads from and writes to a global variable `savedir` to determine the directory for saving output files. Additionally, the function may use the global variable `include_schizo` to modify the output file name.

#### Files Used

None explicitly, although the function depends on the input parameters, which are assumed to be pre-processed and loaded beforehand.

#### Files Created

1. An RDS file containing the formatted principal components (`int_celltype_pcs.rds` or `schizo_int_celltype_pcs.rds`).

#### Functions Called

- `PCAtools::pca()`: Function to perform the principal component analysis.
- `chooseGavishDonoho()`: Function to estimate dimensions using the Gavish-Donoho method.
- `chooseMarchenkoPastur()`: Function to estimate dimensions using the Marchenko-Pastur method.

#### Note

Ensure that the `formatted_meth` and optional `cleaned_mvals` data frames or matrices passed to the function have proper formatting. The function assumes that the rows correspond to CpGs and the columns to samples.

--- 

### `plot_umap`

#### Description

The `plot_umap` function visualizes methylation data in a 2D space using UMAP (Uniform Manifold Approximation and Projection). The function:

1. Reads in principal components (PCs) from saved files.
2. Optionally modifies the input source based on whether the data includes schizophrenia data (`include_schizo`) and whether the data has been cleaned (`cleaned_mvals`).
3. Loads the relevant clinical data to format the UMAP plot.
4. Executes UMAP based on the input data and set parameters, with special considerations for schizophrenia data.
5. Formats UMAP results for plotting, including separating cell types from the sample names.
6. Matches UMAP data with clinical data to add relevant contextual information to the plot.
7. Finally, it visualizes the UMAP results, highlighting different cell types, and saves the resulting plot to a specified location.

#### Parameters

- `cleaned_mvals`: A Boolean indicator specifying whether the methylation values provided have been cleaned.

#### Returns

The function doesn't explicitly return a value but saves:

1. The UMAP results in an RDS file.
2. The UMAP plot in a PNG file.
3. The UMAP plot object in an RDS file.

#### Side Effects

The function reads from and writes to global variables, notably `savedir`, `include_schizo`, and `cell_type_colors`, to determine save locations, handle optional data sources, and use pre-specified color schemes.

#### Files Used

1. Principal components stored in RDS files like `celltype_pcs.rds`, `int_cleaned_celltype_pcs.rds`, and `schizo_int_celltype_pcs.rds`.
2. Clinical data from the specified RDS file path.

#### Files Created

1. An RDS file containing the UMAP results (`celltype_umap_results.rds`).
2. A PNG file for the UMAP plot (`celltype_umap_plot.png` or `schizo_celltype_umap_plot.png`).
3. An RDS file containing the UMAP plot object (`cleaned_celltype_umap_plot.rds`).

#### Functions Called

- `umap()`: The main function used for the Uniform Manifold Approximation and Projection.
- Various `tidyverse` functions for data wrangling, such as `separate()`, `mutate()`, `left_join()`, etc.
- `ggplot()`, `geom_point()`, and associated functions for visualizing UMAP results.

#### Note

This function assumes that relevant packages (like `umap`, `tidyverse`, and `ggplot2`) and global variables (`savedir`, `include_schizo`, and `cell_type_colors`) are available in the R environment. Properly formatted methylation data and clinical data are also prerequisites.

--- 

### `cluster_data`

#### Description
The `cluster_data` function in the script performs k-means clustering on Principal Component (PC) scores obtained from UMAP analysis to examine if schizophrenia cell types cluster with deconvolved cell types. This function allows the user to interactively adjust which cell types to include in the analysis and specifies the number of clusters to be formed. Additionally, it computes the silhouette score for each cluster, providing a measure of how similar an object is to its own cluster compared to other clusters. 

#### Parameters
- `formatted_meth`: The formatted methylation data.
- `cleaned_mvals`: Boolean flag to indicate if m-values should be cleaned or not.

#### Output
- k-means clustering results saved in `km_out`.
- Silhouette scores for the clusters.
- UMAP plots with clusters indicated.

#### Key Functions Used
- `kmeans()`: Performs k-means clustering.
- `silhouette()`: Calculates silhouette scores for clusters.
- `ggplot()`: Generates UMAP plots annotated with k-means cluster and cell type information.

#### Notes
1. The function allows for manual adjustments for the number of clusters and which cell types to include/exclude.
2. The function can be run interactively to adjust these parameters based on each run.
3. It produces UMAP plots saved in the specified directory to visualize the k-means clusters.
4. Two separate clusterings can be created by holding out one of the sorted cell types.

By running `cluster_data`, you'll obtain a comprehensive view of how various cell types are clustering based on their methylation profiles. This will be crucial for any downstream analysis or interpretation related to cell type-specific effects.

---

### <a name="main4">`main`</a>

#### Description

The `main` function serves as the entry point for the workflow, orchestrating the various data loading, preprocessing, and analysis steps for the methylation data. The function:

1. Defines the cell types of interest.
2. Loads cell-type-specific methylation data, with an option for downsampling.
3. Combines the methylation data across different cell types.
4. Optionally includes schizophrenia data in the analysis.
5. Runs Principal Component Analysis (PCA) on the combined methylation data.
6. Visualizes the PCA results using Uniform Manifold Approximation and Projection (UMAP).

#### Parameters

- `downsample`: A Boolean flag that indicates whether to downsample the data. Default is `FALSE`.

#### Returns

The function does not return any value but triggers side effects like printing the status, running PCA, and generating UMAP plots.

#### Side Effects

- Reads from and writes to global variables such as `cleaned_mvals`, `include_schizo`, and `savedir`.
- Outputs various diagnostic information, including the first few rows of data frames, their dimensions, etc.
- Saves PCA and UMAP results in specified directories.

#### Files Used

1. Methylation data for various cell types, loaded using the `load_meth` function.
2. Optional schizophrenia data, loaded using the `load_schizo` function.

#### Files Created

1. PCA results saved as `.rds` files, e.g., `celltype_pcs.rds`, `int_cleaned_celltype_pcs.rds`, or `schizo_int_celltype_pcs.rds`.
2. UMAP plot saved as `.png` file, e.g., `celltype_umap_plot.png` or `schizo_celltype_umap_plot.png`.
3. UMAP results saved as `.rds` files, e.g., `celltype_umap_results.rds`.

#### Functions Called

- `load_meth`: For loading cell-type-specific methylation data.
- `format_meth`: For formatting and combining methylation data from different cell types.
- `load_schizo`: For loading schizophrenia methylation data (optional).
- `run_pca`: For running Principal Component Analysis.
- `plot_umap`: For plotting UMAP visualization.

#### Note

This function assumes that the `load_meth`, `format_meth`, `load_schizo`, `run_pca`, and `plot_umap` functions are defined and available in the current R environment. It also assumes that the required packages and global variables are loaded and initialized.
