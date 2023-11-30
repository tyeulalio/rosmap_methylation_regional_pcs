# Regional PCs Simulation Analysis


# Table of Contents
- [Overview](#overview)
- [Scripts](#scripts)
- [Simulation Results Raw Data](#simulation-results-raw-data)
    - [Accessing the Data](#accessing-the-data)
    - [Extracting the Data](#extracting-the-data)
- [Description of Simulations](#description-of-simulations)
    - [Simulation Parameters](#simulation-parameters)
    - [Range of Tested Parameters](#range-of-tested-parameters)
    - [Fixed or Dependent Parameters](#fixed-or-dependent-parameters)
    - [Output Files](#output-files)
    - [Description of Results File](#description-of-results-file)
        - [Loading the Results](#loading-the-results)
        - [Examining the Contents](#examining-the-contents)
    - [Contents of the List Object](#contents-of-the-list-object)
        - [`avgs`: Averages](#avgs-averages)
        - [`rpcs`: Regional Principal Components](#rpcs-regional-principal-components)
        - [`raw_data`: Simulated Raw Methylation Data](#raw_data-simulated-raw-methylation-data)
        - [`rpcs_full`: Detailed Regional PC Information](#rpcs_full-detailed-regional-pc-information)
        - [Parameters](#parameters)
- [Differential Methylation Analysis](#differential-methylation-analysis)
    - [Scripts for Differential Methylation Analysis](#scripts-for-differential-methylation-analysis)
    - [Methodology](#methodology)
    - [Normalization Procedure](#normalization-procedure)
    - [Statistical Analysis](#statistical-analysis)
    - [Differential Methylation Output Files](#differential-methylation-output-files)
        - [Accessing the DM Data](#accessing-the-dm-data)
        - [Extracting the DM Data](#extracting-the-dm-data)
        - [Storing Results](#storing-results)
        - [Accessing and Understanding the Results Files](#accessing-and-understanding-the-results-files)
        - [Contents of the Differential Methylation Results List Object](#contents-of-the-differential-methylation-results-list-object)
            - [`avgs_res`: List of Results Using Averages](#avgs_res-list-of-results-using-averages)
            - [`rpcs_res`: List of Results Using rPCs](#rpcs_res-list-of-results-using-rpcs)
            - [`num_sig`: Summarized Results](#num_sig-summarized-results)


# Overview
This README document provides an overview of the Regional PCs (Principal Components) Simulation Analysis. It includes details about the scripts used, accessing and extracting the raw simulation data, and a description of the simulation parameters and results.

# Scripts

- `RRBSmodeling.R`: Modeling script from the RRBS-sim paper.
- `RRBSsimulation.R`: Simulation script from the RRBS-sim paper.
- `00_simulation_test.R`: Script that runs the simulation analysis using the RRBS-sim scripts.
- `01_run_dmp.R`: Runs differential methylation analysis 
- `process_results_example.R`: Example script for processing the simulation results.

# Simulation Results Raw Data

## Accessing the Data
The raw data from the simulation results are available on Google Drive. Access can be requested at the following link: [Google Drive Folder](https://drive.google.com/drive/u/0/folders/19hGbVr4HpASlUPAXpjQJ1N4rXHy8VG-k).

## Extracting the Data
The data is compressed in a tar.gz format. To extract the files, use the following command in the terminal:

```bash
tar -zxf file.tar.gz --directory /path/to/directory
```

For more information on extracting files, visit [CyberCiti: How to extract tar files to a specific directory](https://www.cyberciti.biz/faq/howto-extract-tar-file-to-specific-directory-on-unixlinux/).


# Description of Simulations

## Simulation Parameters
The simulation analysis was conducted using the following parameters:

- `num_sites`: Number of CpG sites to simulate in the region.
- `num_samples`: Number of samples for which data are simulated.
- `percent_meth_difference`: Difference in methylation percentage between control and case at DMR sites.
- `dmr_length`: Length of differentially methylated regions (not equivalent to the number of differentially methylated CpGs).
- `percent_sites_dm`: Percentage of CpGs that should be differentially methylated.

Example settings:
```plaintext
num_sites = 50
num_samples = 500
percent_meth_difference = 0.1
dmr_length = 25
percent_sites_dm = 0.1
```

## Range of Tested Parameters
The following ranges of values were tested for the parameters:

```plaintext
num_sites_range <- c(20, 50)
num_samples_range <- c(50, 500, 5000)
percent_meth_difference_range <- seq(0.1, 0.9, 0.1)
percent_sites_dm_range <- seq(0, 1, 0.25)
```

## Fixed or Dependent Parameters
- `N = 1000`: Number of regions (genes) simulated.
- `dmr_length = num_sites / 2`: Dependent on `num_sites`.

`N` is the number of regions (genes) that were simulated using the defined parameters.


## Output Files
For each parameter combination, we simulated N=1000 gene regions. Each set's results were stored in a separate output file. The file naming convention is as follows:

```R
filename <- paste0("simulated_data",
                    "_numSites", num_sites,
                    "_numSamples", num_samples,
                    "_pct_meth_diff", percent_meth_difference,
                    "_dmr_length", dmr_length,
                    "_pct_sites_dm", percent_sites_dm,
                    "_N", N)
savefile <- paste0(savedir, filename, ".rds")
```

The `filename` contains the parameters for the run, and `savedir` is the directory where the output files are stored.


## Description of Results File

Each results file is stored as an `rds` file, encapsulating a list object with multiple elements related to the simulation results.

### Loading the Results
To load a results file, use the following R command, with `savefile` being the path to your specific `.rds` file:

```R
results <- readRDS(savefile)
```

### Examining the Contents
To understand the structure and components of the loaded `results` list, you can inspect its elements:

```R
names(results)
```

### Contents of the List Object
The `results` list comprises several elements, each offering different insights into the simulated data:

1. **avgs**: A matrix containing average values across all simulated regions. This matrix provides a summary view of the average methylation levels across the gene regions simulated.

2. **rpcs**: This data frame holds the regional principal components (rPCs) for all simulated regions.

3. **raw_data**: A list containing detailed data for each of the simulated regions. 

4. **rpcs_full**: A comprehensive list that includes the full data on regional principal components. 

5. **parameters**: This list element contains the parameters used during the simulation run.

### Detailed Element Description
Each of the above elements plays a vital role in the overall analysis of the simulation data. Below is a detailed description of what each element represents and how it can be utilized:

#### `avgs`: Averages

The `avgs` matrix in the results file provides a summary of the average values calculated for each simulated region. It's structured as follows:

- **Rows**: Each row in this matrix represents a unique region simulated in the study. The total number of rows corresponds to `N`, indicating the total number of regions (genes) that were simulated.
  
- **Columns**: Every column corresponds to a simulated "sample" or individual. The number of columns is equal to `num_samples`, which is the total number of samples included in the simulation. These samples are labeled to indicate their group allocation in the study, with names starting with "case" or "control". 

- **Data Representation**: The values in the matrix represent the average methylation levels. These averages are calculated for each region and each sample. Specifically, the average methylation value is computed across all CpGs that were part of the simulation for that region. 


#### `rpcs`: Regional Principal Components

This data frame encompasses the regional Principal Components (rPCs) for each simulated region. The structure is as follows:

- **Rows**: Each row represents a distinct rPC. It's important to note that a single region may be associated with multiple rPCs. The naming convention for these rows includes the region number followed by the PC number. The region numbers in the row names of the `rpcs` data frame match those in the `avgs` matrix.

- **Columns**: The columns correspond to the simulated "samples" (or individuals). The naming convention for these samples is consistent with that used in the `avgs` matrix. Thus, every column in the `rpcs` data frame has a directly corresponding column in the `avgs` matrix.


#### `raw_data`: Simulated Raw Methylation Data

This section of the results file contains the raw methylation data for each of the simulated regions, as generated by the RRBS-sim tools. The structure is as follows:

- **Structure**: `raw_data` is a list where each element represents a unique region. There are `N` such elements, corresponding to the number of regions simulated.

- **Contents of Each Region**:
  1. **Methylation Data Frame**: The first element within each region is a data frame. Each row corresponds to a specific CpG site, and the columns represent different samples. This data frame includes:
     - `chrpos`: Chromosome position of the CpG site.
     - `dmsite`: Status of differential methylation, indicated as -1 (hypomethylated), 0 (no differential methylation), or 1 (hypermethylated).
     - `mprob` and `mprob.diff`: Methylation probability for cases and controls.
     - Percent Methylation (PM) and Sequencing Coverage (SC) for each CpG site across samples.

  2. **Differentially Methylated Regions (DMRs) Matrix**: The second element is a matrix detailing the differentially methylated regions within the simulated region. Each row represents a DMR with the following columns:
     - `dmr_row_start`: Row position in the Methylation Data Frame where the DMR begins.
     - `dmr_chrpos`: Chromosome position where the DMR starts.
     - `dmr_number_sites`: The number of consecutively differentially methylated sites comprising the DMR.


#### `rpcs_full`: Detailed Regional PC Information

This section comprises comprehensive details on the regional Principal Components (rPCs) for each simulated region. The structure of `rpcs_full` is a list, where each list element corresponds to one region. The index of each element aligns with the respective region number.

Each region's entry in `rpcs_full` is a list containing four distinct components:

1. **Regional PCs Data Frame (`regional_pcs`)**: 
   - This data frame houses the regional PCs for the specific region. 
   - The rows represent individual PCs, while the columns correspond to the samples.
   - The combination of these regional PCs across all regions forms the `rpcs` matrix in the main results file.

2. **Percent Variance Explained (`percent_variance`)**:
   - This component is a data frame that details the percentage of variance each PC explains within the region.
   - Each row represents one PC, providing insights into the relative importance or contribution of each PC to the region's variance.

3. **Loadings (`loadings`)**:
   - A list of loadings for each PC within the region.
   - These loadings offer a quantitative measure of how strongly each CpG site contributes to the respective PC.

4. **Region Information (`info`)**:
   - A matrix that includes essential information about each region:
     - Region name: Identifies the region.
     - Estimated dimension: Determined using the Gavish-Donoho method, which helps discern meaningful PCs from those representing noise.
     - Number of CpGs: The total count of CpG sites in the region.


#### parameters
List containing the parameters used for these simulations. List contains num_sites, num_samples, percent_meth_different, dmr_length, percent_sites_dm, and N.

Your section on Differential Methylation in the README is clear but can be further refined for better clarity and detail. Here's a revised version:


# Differential Methylation Analysis

## Scripts for Differential Methylation Analysis

The differential methylation analysis was conducted using the `01_run_dmp.R` script. This script was applied to two different types of methylation summaries: the averages matrix and the regional principal components (rPCs) matrix. Both analyses were performed separately to assess differential methylation patterns across the summarized regions.

## Methodology

For the analysis, linear models were employed using the `lmFit` function from the `limma` R package. This function is specifically designed to fit linear models to high-dimensional data like methylation data, allowing for efficient and robust analysis of differential methylation across regions.

## Normalization Procedure

Before running the differential methylation analysis, the methylation summaries were normalized. This normalization was carried out using a rank-based inverse normal transformation, implemented via the `RankNorm` function from the `RNOmni` package. This step is crucial to reduce the impact of outliers and ensure that the methylation data conforms more closely to a normal distribution, which is a key assumption of the linear modeling approach used.

## Statistical Analysis

Once the linear models were fitted, an empirical Bayes method was employed to obtain more stable and reliable estimates of the test statistics. This method, integral to the `limma` package, enhances the power of the tests by borrowing information across regions, thereby providing more precise and robust insights into differential methylation patterns.

Your section detailing the output files for differential methylation and how to access them is well-structured. Here’s a slightly refined version to enhance clarity and detail:

## Differential Methylation Output Files

### Accessing the DM Data
The output for the differential methylation analysis can be found on Google Drive. Access can be requested at the following link: [Google Drive Folder](https://drive.google.com/drive/u/0/folders/19hGbVr4HpASlUPAXpjQJ1N4rXHy8VG-k).

### Extracting the DM Data
The data is compressed in a tar.gz format. To extract the files, use the following command in the terminal:

```bash
tar -zxf file.tar.gz --directory /path/to/directory
```

For more information on extracting files, visit [CyberCiti: How to extract tar files to a specific directory](https://www.cyberciti.biz/faq/howto-extract-tar-file-to-specific-directory-on-unixlinux/).

### Storing Results

The differential methylation analysis for each combination of simulation parameters results in the generation of output files for N=1000 simulated gene regions. These output files follow a specific naming convention and are stored in a designated directory. The naming convention for the output files is as follows:

```R
filename <- paste0("DM_results",
                    "_numSites", num_sites,
                    "_numSamples", num_samples,
                    "_pct_meth_diff", percent_meth_difference,
                    "_dmr_length", dmr_length,
                    "_pct_sites_dm", percent_sites_dm,
                    "_N", N)
savedir <- "01_dmp_results/"
savefile <- paste0(savedir, filename, ".rds")
```

In this naming convention, `filename` comprehensively includes all the parameters used in the analysis. The `savedir` variable specifies the directory path where these output files are stored.

## Accessing and Understanding the Results Files

Each differential methylation results file is saved in an `rds` format, which encapsulates a list object. This list object contains multiple elements that provide insights into the differential methylation analysis.

### Loading the Results

To access the data from a results file, you can use the R programming language. The following command demonstrates how to load a results file, with `savefile` being the path to the specific `.rds` file:

```R
results <- readRDS(savefile)
```

### Examining the File Contents

Once the results file is loaded into R, you can explore its structure and contents to understand the various components of the differential methylation analysis. The names of the elements within the `results` list object can be inspected using this command:

```R
names(results)
```

This command will display the different elements stored within the list object, providing an overview of the data structures and results included in the file.

Your README section detailing the contents of the list object in the differential methylation results file is clear and well-structured. Here's a refined version for better clarity and detail:

### Contents of the Differential Methylation Results List Object

The `results` list object, obtained from each differential methylation `.rds` file, comprises several key elements:

1. **avgs_res**: This list contains the results of differential methylation analysis for regions summarized as averages.

2. **rpcs_res**: Similar to `avgs_res`, this list presents the differential methylation results for regions summarized as regional principal components (rPCs).

3. **num_sig**: This element is a tibble (a type of data frame in R) that provides a summarized count of significant differential methylation results. 


### Detailed Element Description for Differential Methylation Results

#### `avgs_res`: List of Results Using Averages
This list comprises two key elements that provide insights into the differential methylation analysis using averages:

1. **full_meth_dims**: A pair of integers detailing the dimensions of the methylation matrix used for differential methylation analysis (formatted as rows, columns). If the dimensions are (N, num_samples), it implies that all simulated regions were included. A discrepancy in these numbers suggests that certain regions were excluded due to NA values, likely resulting from issues in the simulation process. The number of rows indicates the actual number of regions tested.

2. **results**: A data frame with the differential methylation analysis outcome for each region. Derived from the `toptable` function of the `limma` package, it includes several key metrics:
    - `logFC`: Log fold change between cases and controls for the region.
    - `AveExpr`: Average log2 methylation level for the region.
    - `t`: Moderated t-statistic.
    - `P.Value`: Raw p-value.
    - `adj.P.Val`: Raw p-value (no adjustment was specified during run).
    - `B`: Log-odds that the region is differentially methylated.
    - `bh_pval`: Benjamini-Hochberg adjusted p-value.
    - `bf_pval`: Bonferroni adjusted p-value.

#### `rpcs_res`: List of Results Using rPCs
Similarly structured as `avgs_res`, this list contains differential methylation results based on regional PCs:

1. **full_meth_dims**: Dimensions of the methylation matrix, similar to `avgs_res`.

2. **results**: Differential methylation analysis outcome for each region, with same columns as in `avgs_res`.

#### `num_sig`: Summarized Results
This tibble provides a summary of significant differential methylation results, categorized by different adjustment methods and summary types. Each row represents a set of simulation parameters, and the columns include:

1. `num_sites`: Number of sites in the simulated region.
2. `num_samples`: Number of samples in the simulation.
3. `percent_meth_difference`: Methylation difference percentage parameter.
4. `dmr_length`: Length of differentially methylated regions parameter.
5. `percent_sites_dm`: Percentage of differentially methylated sites parameter.
6. `N`: Total number of regions simulated.
7. `num_rpcs`: Total number of rPCs across all simulated regions.
8. `bf_sig_avgs`: Count of significant regions using *Bonferroni adjusted* p-value with *averages* summary.
9. `bf_sig_rpcs`: Count of significant regions using *Bonferroni adjusted* p-value with *rPCs* summary.
10. `bh_sig_avgs`: Count of significant regions using *Benjamini-Hochberg* adjusted p-value with *averages* summary.
11. `bh_sig_rpcs`: Count of significant regions using *Benjamini-Hochberg* adjusted p-value with *rPCs* summary.
12. `uc_sig_avgs`: Count of significant regions using *unadjusted* p-value with *averages* summary.
13. `uc_sig_rpcs`: Count of significant regions using *unadjusted* p-value with *rPCs* summary.

The p-value threshold for significance used in these analyses was set at **0.05**.









