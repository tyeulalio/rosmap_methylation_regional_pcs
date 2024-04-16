# Differential Methylation Analysis

This pipeline contains a set of functions designed to perform differential methylation analysis and analyze the results.

## Table of Contents

- [01_dmp_analysis.R](#script-01)
    - [load_pheno](#load_pheno1)
    - [load_meth](#load_meth1)
    - [match_cpgs_to_regions](#match_cpgs_to_regions)
    - [get_global_pcs](#get_global_pcs)
    - [run_dmp](#run_dmp)
    - [plot_results](#plot_results)
    - [plot_dmp_counts_barplot](#plot_dmp_counts_barplot)
    - [plot_gene_counts_barplot](#plot_gene_counts_barplot)
    - [main-function](#main-function1)

                                    

## Script 01                                    
## 01_dmp_analysis.R

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

Below are descriptions and usage guidelines for each of the functions in this script.


### <a name="load_pheno1">`load_pheno`</a>

#### Description

The `load_pheno` function is designed to load and preprocess phenotype data specifically formatted for analyses involving different cell types. It reads phenotype data from a predefined RDS file, then selects, transforms, and formats relevant phenotype variables to facilitate downstream analysis. This function handles binary transformations of diagnosis and genotype variables, converts data types appropriately, and filters out records with missing essential data.

#### Parameters

- `cell_type`: A character string specifying the cell type for which phenotype data is to be prepared. Although not directly used to vary the input file in this implementation, it indicates potential flexibility for future adaptations where phenotype data might vary by cell type.

#### Returns

- `sub_pheno`: A data frame with selected and transformed phenotype variables, ready for use in statistical analyses. This includes both numeric and factor (categorical) variables, specifically formatted for further analysis.

#### Details

- The function starts by reading an RDS file containing phenotype data. It assumes a fixed filename and directory, which could be parameterized in future versions.
- Selects relevant phenotype variables and transforms certain key variables for analysis:
  - Converts the diagnosis into a binary format.
  - Transforms APOE genotype into a numeric binary format indicating a specific genotype.
- Adjusts data types across the dataset:
  - Numeric variables are converted from characters for mathematical operations.
  - Categorical variables are factored for statistical analysis compatibility.
- Filters out entries with missing essential demographic data (age at death and sex).
- Adds a derived categorical variable for a specific clinical score to highlight advanced disease stages.
- Outputs the processed data frame for downstream analytical tasks.

#### Example Usage

```R
# Load phenotype data for astrocyte-specific analysis
astro_pheno_data <- load_pheno("astro")
```

---

### <a name="load_meth1">`load_meth`</a>

#### Description

The `load_meth` function loads methylation data from stored RDS files, formatted either as CpG-level methylation values or as summary statistics such as averages or principal components. This function is versatile, handling data files according to specified region types, summary methods, and cell types. It adapts the file-loading process to different formats, ensuring that the output is uniformly structured for subsequent analysis.

#### Parameters

- `region_type`: A string indicating the genomic region of interest (e.g., 'promoter', 'gene_body').
- `summary_type`: A string specifying the type of methylation data summary ('cpgs', 'avgs', 'pcs').
- `cleaned_str`: A string that is part of the file name indicating whether the data is cleaned.
- `cell_type`: A string that specifies the cell type of the data being loaded (e.g., 'astrocytes', 'neurons').

#### Returns

- `meth`: A DataFrame containing methylation data formatted according to the specified summary and region types. It may contain CpG IDs, methylation beta values, or principal component scores, depending on the `summary_type`.

#### Details

- The function constructs the file path based on the input parameters to locate the correct RDS file. This includes handling different file naming conventions for raw CpG-level data and summarized data.
- For CpG summary data, the function directly reads the RDS file.
- For average or principal component summaries, additional data formatting steps are performed to ensure consistency across datasets:
  - For averages, unnecessary columns like duration times are removed.
  - For principal components summaries, non-PC related columns are removed to focus the dataset on relevant variables.
- The data is checked and formatted to handle different structures based on the `summary_type`, ensuring that outputs are always ready for downstream processing without further need for data wrangling.

#### Example Usage

```R
# Load averaged methylation data for astrocytes in the gene body region
astro_gene_body_avg_meth <- load_meth('gene_body', 'avgs', 'cleaned_', 'astrocytes')

# Load principal components of methylation data for neurons in promoter regions
neuron_promoter_pcs_meth <- load_meth('promoter', 'pcs', 'cleaned_', 'neurons')
```

---

### <a name="match_cpgs_to_regions">`match_cpgs_to_regions`</a>

#### Description

The `match_cpgs_to_regions` function is designed to associate CpG sites with specific genomic regions and their corresponding gene mappings. This function facilitates the integration of CpG methylation data with genomic annotations, ensuring that methylation data can be analyzed within the context of genomic regions such as promoters or gene bodies.

#### Parameters

- `meth`: A matrix or DataFrame containing methylation data, where rows are CpG sites and columns are samples.
- `region_type`: A string indicating the specific genomic region (e.g., 'promoter', 'gene_body') for which gene-CpG mappings are required.
- `cell_type`: A string specifying the cell type, used to filter data relevant to specific cellular contexts.

#### Returns

- `meth`: A DataFrame where methylation data has been subset to include only CpGs associated with the specified genomic region and cell type.

#### Details

- The function starts by loading a gene-CpG mapping file specific to the chosen `region_type` and `cell_type`, which contains mappings of CpGs to their respective genes.
- It then loads a CpG position map to fetch the CpG IDs corresponding to their positions.
- Both mapping files are joined to form an extended map which is then used to filter the input methylation data.
- The filtered methylation data contains only those CpGs that are found within the specified genomic regions, making it ready for detailed region-specific analysis.

#### Example Usage

```R
# Load methylation data and match CpGs to promoter regions for astrocytes
astro_promoter_meth <- load_meth(data_matrix, 'promoter', 'astrocytes')
```

---

### <a name="get_global_pcs">`get_global_pcs`</a>

#### Description

The `get_global_pcs` function retrieves global principal components (PCs) that have been pre-computed for methylation data across different cell types. This function is crucial for adjusting methylation data by removing variation attributable to known and unknown confounders. It allows the inclusion or exclusion of global PCs based on their correlation with specific covariates, ensuring that only relevant biological signals are retained for downstream analyses.

#### Parameters

- `cell_type`: A string specifying the cell type for which global PCs are to be retrieved. If 'bulk' is specified, it retrieves global PCs calculated for bulk tissue methylation data.
- `include_global`: A boolean indicating whether to include global PCs in the output. If set to `FALSE`, the function returns `NULL`.

#### Returns

- `global_pcs`: A matrix of principal components that can be used for batch correction in downstream analyses. If `include_global` is `FALSE`, no data is returned.

#### Details

- The function first checks if global PCs should be included based on the `include_global` flag.
- Depending on the `cell_type`, it selects the appropriate file containing the global PCs, which were computed during batch correction stages of data preprocessing.
- It then filters the PCs based on the significance of their correlation with traits of interest and other covariates like batch effects, post-mortem interval (PMI), and study design. This step is crucial to ensure that only PCs affecting the outcome variable are retained, thereby protecting the integrity of the biological signals in the methylation data.

#### Example Usage

```R
# Retrieve global PCs for astrocytes, including them in the analysis
astrocyte_pcs <- get_global_pcs('astrocytes', TRUE)
```

---

### <a name="run_dmp">`run_dmp`</a>

#### Description

The `run_dmp` function runs the analysis to identify differentially methylated probes (DMP) associated with specific traits. It integrates phenotype data, methylation data, and optionally, global principal components (PCs), to perform a comprehensive analysis using the `limma` package. The function ensures data alignment, normalizes the methylation values, constructs a model including necessary covariates, and performs statistical testing.

#### Parameters

- `meth`: A matrix containing methylation beta values, with CpG sites as rows and samples as columns.
- `pheno`: A data frame containing phenotype data, where rows correspond to samples.
- `global_pcs`: A matrix of global principal components to be included as covariates in the model.
- `region_type`: A string indicating the genomic region type being analyzed (e.g., 'promoter', 'gene body').
- `summary_type`: A string specifying the type of summary applied to the methylation data (e.g., 'averages', 'PCs').
- `cleaned_str`: A string to append to filenames indicating the data cleaning level.
- `trait`: A string specifying the trait column in the phenotype data used for analysis.
- `cell_type`: A string specifying the cell type being analyzed.
- `shuffle_pheno`: A boolean indicating whether the phenotype data should be shuffled for permutation testing.
- `shuffle_num`: An integer specifying how many times the phenotype data should be shuffled.

#### Returns

- The function saves the results of the differential methylation analysis to a specified directory and returns a status indicator. The saved results include CpG sites with significant p-values after multiple testing correction.

#### Details

The function follows these steps:

1. Attach global PCs to phenotype data.
2. Ensure the phenotype and methylation data are in the same order.
3. Construct and normalize methylation and phenotype matrices.
4. Formulate the linear model with specified covariates including age, sex, batch effects, and any global PCs.
5. Fit the model using the `limma` package's `lmFit` and `eBayes` functions.
6. Adjust p-values using the Benjamini-Hochberg method.
7. Filter results based on adjusted p-values and save the significant results to a file.

#### Example Usage

```R
# Run DMP analysis for astrocytes using promoter region data
results <- run_dmp(meth_data, pheno_data, global_pcs, "promoter", "PCs", "clean", "Age", "astrocytes", FALSE, 0)
```

---

### <a name="check_results">`check_results`</a>

#### Description

The `check_results` function is designed to assess the results from differential methylation analysis across multiple runs. It facilitates batch processing of different combinations of region types, summary types, traits, cleaning strings, and cell types, thereby enabling a comprehensive review of the analysis results. The function involves several sub-processes: loading result files, formatting results with gene symbols, and checking the overlap of CpG sites across detected regions.

#### Parameters

- None directly; all variables are set within the function.

#### Returns

- The function does not return values but saves formatted result files to a specified directory.

#### Details

The function performs the following steps:

1. **Define Parameter Choices**: It establishes combinations of region types, summary types, traits, and cell types for which results will be reviewed.
2. **Load Gene and CpG Mappings**: It reads gene to CpG mapping files to associate results with gene symbols and positional data.
3. **Format Results**: For each combination of parameters, it loads the results file, formats it by separating principal components from gene IDs, adds gene symbols, and saves the formatted results.
4. **CpG Overlap Analysis**: Checks how CpG sites overlap the regions detected using principal components and other summary types. It maps CpG sites to genes and saves these mappings.
5. **Batch Execution**: Applies the formatting and mapping functions across all parameter combinations using `apply`.

#### Example Usage

```R
# To check and process results after running differential methylation analyses:
check_results()
```

---

Here's the README description for the `plot_results` function, formatted in markdown:

### <a name="plot_results">`plot_results`</a>

#### Description

The `plot_results` function visualizes the results from differential methylation analysis across various data configurations. This function is designed to process and plot results for multiple combinations of region types, summary types, cell types, and traits. It leverages extensive data manipulation to organize, format, and prepare data for visual representation.

#### Parameters

- None directly; variables are set within the function for batch processing.

#### Returns

- The function does not return values but saves various plots and result summaries to specified directories.

#### Details

The function performs the following tasks:

1. **Load Gene and CpG Mappings**: It reads mappings from various sources to associate DMP results with gene symbols and CpG site details.
2. **Data Preparation**: For each combination of parameters, it formats the results, merges them with gene symbols, and checks the overlap of CpG sites with detected regions.
3. **Result Visualization**: Generates various plots to visualize the analysis results, including bar plots for significant DMP counts and Manhattan plots for visualizing the genome-wide distribution of p-values.
4. **Batch Processing**: Applies data preparation and plotting tasks across all combinations of parameters using `apply`.
5. **Save Outputs**: Writes processed results and plots to specified directories for further use or publication.

#### Example Usage

```R
# To visualize DMP analysis results across various data configurations:
plot_results()
```

---

### <a name="plot_dmp_counts_barplot">`plot_dmp_counts_barplot`</a>

#### Description

The `plot_dmp_counts_barplot` function visualizes the counts of significant differential methylation patterns (DMPs) across different configurations such as cell types, summary types, and genomic regions. This visualization helps in assessing the impact of different methylation summarization techniques and the extent of differential methylation across various biological categories.

#### Parameters

- `combined_df`: A data frame containing differential methylation results including features like p-values, genomic regions, cell types, and traits.
- `correction_type`: A string indicating the type of p-value correction applied ('bh' for Benjamini-Hochberg or 'bf' for Bonferroni).
- `summary_type_colors`: A named vector of colors to be used for each summary type in the plots.
- `plotting_dir`: The directory path where the plots should be saved.
- `full_array`: A boolean indicating whether the full array of CpGs was used (`TRUE`) or only those overlapping gene regions (`FALSE`).

#### Returns

- The function does not return a value but generates and saves bar plot images in the specified directory.

#### Details

- **Define Parameters and Format Results**: The function sets up the plotting environment based on the correction type and array usage (full or matched). It prepares sub-directories for saving plots.
- **Count Significant Features**: It calculates the number of significant DMPs for each configuration of summary type, region type, and cell type.
- **Count Total Features**: It computes the total number of features tested for each configuration.
- **Compute Percentage of Significant Results**: It calculates the percentage of significant results for each configuration and formats the data for further analysis.
- **Compare Feature Proportions**: This function compares the proportion of features that were found to be significant across different configurations using Fisher's Test. Adjustments for multiple testing are made to determine significant differences. These comparisons help in understanding how each summary type performs in terms of identifying significant features.
- **Plot Feature Proportions and Counts**: The function plots the percentage of significant features as well as the total counts of features as bar plots. This includes plotting for specific traits across different region types and cell types.
- **Additional Visualizations**: It also plots the full gene counts with traits and generates bar plots comparing significant DM results for the full gene region across all traits. This provides a visual comparison across traits, and these plots are saved in supplementary figures of the regionalpcs manuscript.

#### Example Usage

```R
# To visualize DMP analysis results:
plot_dmp_counts_barplot(combined_df, 'bh', summary_type_colors, 'path/to/plotting_dir', FALSE)
```

---


### <a name="plot_gene_counts_barplot">`plot_gene_counts_barplot`</a>

#### Description

The `plot_gene_counts_barplot` function is designed to visualize the gene-level differential methylation results across different configurations such as cell types, summary types, and genomic regions. This function focuses on the gene count that shows significant differential methylation under various conditions.

#### Parameters

- `combined_df`: A data frame containing differential methylation results including gene ids, p-values, genomic regions, cell types, and traits.
- `correction_type`: A string specifying the type of p-value correction applied ('bh' for Benjamini-Hochberg or 'bf' for Bonferroni).
- `summary_type_colors`: A named vector of colors to be used for each summary type in the plots.
- `plotting_dir`: The directory path where the plots should be saved.
- `full_array`: A boolean indicating whether the full array of CpGs was used (`TRUE`) or only those overlapping gene regions (`FALSE`).

#### Returns

- The function does not return a value but generates and saves bar plot images in the specified directory.

#### Details

1. **Setup Plotting Environment**: Initializes the directory structure based on the correction type and full array usage, and sets the adjusted p-value based on the correction type specified.
2. **Data Filtering**: Filters the results based on the `full_array` parameter to either use all CpGs or just those overlapping gene regions.
3. **Data Preparation for Gene-Level Analysis**: Transforms the data to focus on gene counts rather than CpG or region counts, specifically looking at gene-probe pairs and adjusting for significant CpG counts per gene.
4. **INTACT GWAS Integration**: Revisit the analysis post-INTACT processing to identify overlaps with final GWAS-integration results.
5. **Upset Plot Visualization**: Generate Upset plots to visualize the intersections of DM genes across cell type groups, including results for both 'averages' and 'rPCs' summary types.
6. **Comparison of Region Types**: Employ Upset plots to analyze the differential methylation across various genomic region types.
7. **Summary of Gene-Level Counts**: Summarize counts of significant DM results, total genes tested, and the percentage of significant findings for each analysis run.
8. **Proportional Comparison Across Runs**: Use proportion tests and apply multiple testing corrections to compare the proportion of DM results across different configurations.
9. **Barplot of Significant Gene Counts**: Create a barplot showing the count of significant DM genes for each summary type within the full gene region.
10. **Comparison with Known AD Genes**: Contrast findings with known Alzheimer's disease genes from Open Targets and NIAGDS gene sets.
11. **Correlation Analysis**: Test and plot the correlation between DM p-values and Open Targets overall association scores.
12. **OpenTargets Tile Plot**: Generate a tile plot to display DM genes based on a specific percentile of Open Targets Association Scores, intended for publication figures.


#### Example Usage

```R
# To visualize gene-level DMP analysis results:
plot_gene_counts_barplot(combined_df, 'bh', summary_type_colors, 'path/to/plotting_dir', FALSE)
```

---

### <a name="main-function">Main Function for Differential Methylation Analysis</a>

#### Description

The `main` function orchestrates a comprehensive differential methylation analysis by coordinating data loading, preprocessing, analysis, and visualization across multiple configurations and conditions. It is designed to handle various types of genomic regions, summary metrics, and phenotypic traits.

#### Parameters

The function does not have parameters; it uses predefined settings within its scope for the analysis.

#### Returns

This function does not return a value explicitly but outputs progress logs and generates final plots showing the results of the differential methylation analysis.

#### Details

1. **Parameter Setup**: Initializes the variables for different configurations such as genomic regions (`region_types`), summary types (`summary_types`), traits (`traits`), data cleaning strategies (`cleaned_strs`), and cell types (`cell_types`). Constructs a dataframe `runs` containing all possible combinations of these parameters for exhaustive analysis.

2. **Interactive Testing Setup**: For debugging or demonstration purposes, selects the first configuration from `runs` for immediate processing.

3. **Row-wise Processing**:
   - Each configuration is processed row by row. 
   - Loads phenotypic data specific to the cell type involved.
   - Loads methylation data tailored to the region type, summary type, cleaning string, and cell type.
   - For CpG summary types, matches CpGs to their corresponding genomic regions.
   - Optionally includes global principal components as covariates if specified.
   - Executes differential methylation analysis (`run_dmp`) using the loaded data along with additional parameters like trait and whether to shuffle phenotypes for control tests.

4. **Progress Tracking**: Prints progress updates for each configuration, indicating completion and ongoing status.

5. **Results Checking and Plotting**:
   - After processing all configurations, checks the integrity and completeness of the results.
   - Calls `plot_results` to visualize the outcomes of the differential methylation analysis.

#### Usage Example

```R
# To run the main differential methylation analysis pipeline:
main()
```


