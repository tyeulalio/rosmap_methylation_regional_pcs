# Plotting Scripts

This directory contains scripts to plot figures for the Regionalpcs manuscript. Some script numbers are not included because they weren't used in the final analysis.

## Table of Contents


- [02_plot_qtl_gene_counts.R](#script-01)
- [04_create_upset_plots.R](#script-02)
- [05_plot_pc_qtl.R](#script-03)
- [06_plot_celltype_summary_figure.R](#script-04)
- [07_check_eqtls.R](#script-05)
- [08_plot_dm_figure.R](#script-06)
- [10_plot_supp_figures.R](#script-07)



## Script 01                                    
## 02_plot_qtl_gene_counts.R

This script generates a summary of counts across various steps in the Quantitative Trait Loci (QTL) analysis and Genome-Wide Association Studies (GWAS) integration pipeline. It counts significant QTLs, fine-mapped QTLs, colocalization events, PTWAS results, and INTACT results. Additionally, it creates plots to visualize the progression and results of the analysis.

## Overview

- Defines functions to count significant results and colocalization events at various stages in the QTL and GWAS analysis pipeline.
- Combines and summarizes these counts to provide a comprehensive view of the analysis outcomes.
- Produces plots illustrating the gene counts at different steps in the analysis pipeline.
- Outputs formatted results and supplementary tables for additional analysis or publication.

## Script Execution

### Input Data

- `summary_types`: Types of summaries for analysis.
- `region_types`: Regions of interest in the analysis.
- `cell_types`: Cell types to consider.
- `savedir`: Directory for saving output results.
- `plotdir`: Directory for saving plots.

### Functions

#### `get_total_counts(runname)`

- Counts the total number of genes, features, and variants tested in QTL analysis.
- Outputs a data frame with the counts.

#### `get_sig_mapped_qtls(runname)`

- Counts significant QTLs after QTL mapping.
- Outputs a data frame with counts of genes, features, and variants.

#### `get_finemapped_qtls(runname)`

- Counts significant fine-mapped QTLs.
- Outputs a data frame with counts of genes, features, and variants.

#### `get_colocalization_counts(runname, coloc_thresh)`

- Counts colocalization events with a specified threshold for significance.
- Outputs a data frame with counts of genes, features, and variants.

#### `get_ptwas_counts(runname, ptwas_thresh)`

- Counts significant results from PTWAS.
- Outputs a data frame with counts of genes and features.

#### `get_intact_counts(runname, intact_thresh)`

- Counts significant results from INTACT analysis.
- Outputs a data frame with counts of genes and features.

### Execution Details

1. **Process Runs**:
   - Iterates over a list of runs, counting significant results and colocalization events.
   - Combines and summarizes these counts to create a comprehensive data frame.

2. **Generate Plots**:
   - Creates plots illustrating significant results at various stages of the QTL and GWAS pipeline.
   - Saves the plots to the designated output directory.

### Notes

- Ensure all input files are accessible and correctly formatted.
- Adjust the `summary_types`, `region_types`, `cell_types`, and thresholds as needed.
- Check the output and plot directories to ensure they exist or create them if necessary.

---

## Script 02                                   
## 04_create_upset_plots.R

This script creates Upset plots to visualize overlaps among genes in various analyses, including QTL mapping, fine-mapping, colocalization, PTWAS, and INTACT results. It generates and saves multiple plots that compare cell types and summary types in the context of QTL and GWAS analyses. The script also combines individual plots into comprehensive visualizations.

## Overview

- Defines functions to extract gene lists from different stages of QTL analysis.
- Generates Upset plots to represent overlaps among gene lists for various cell types and summary types.
- Combines multiple Upset plots into single visualizations for broader comparisons.
- Outputs plots in both PDF and PNG formats, allowing easy integration with other analysis tools or reports.

## Script Details

### Input Data

- `savedir`: The directory for saving output plots and data.
- `summary_types`: Types of summaries for analysis.
- `region_types`: Regions of interest.
- `cell_types`: Cell types to consider.
- `coloc_thresh`: Threshold for colocalization significance.
- `intact_thresh`: Threshold for INTACT significance.
- `ptwas_thresh`: Threshold for PTWAS significance.

### Functions

#### `get_total_features(runname)`

- Returns all features tested in the QTL analysis.

#### `get_qtl_genes(cell_type, region_type, summary_type)`

- Returns a list of significant QTL genes for a given cell type, region type, and summary type.

#### `get_finemapped_genes(cell_type, region_type, summary_type)`

- Returns a list of fine-mapped genes for a given cell type, region type, and summary type.

#### `get_coloc_genes(cell_type, region_type, summary_type)`

- Returns a list of colocalized genes based on a given colocalization threshold.

#### `get_ptwas_genes(cell_type, region_type, summary_type)`

- Returns a list of genes from PTWAS analysis with significant z-scores.

#### `get_intact_genes(cell_type, region_type, summary_type)`

- Returns a list of genes from INTACT analysis with significant probabilities.

#### `plot_upset(gene_lists, savefile, colors_map, angle, width, height)`

- Creates and saves an Upset plot for the provided gene lists.
- Utilizes the `ComplexUpset` library for plot generation.

#### `combine_plots(p0, p1, p2, p3, p4, region_type, summary_type, colors_map)`

- Combines multiple plots into a single visualization.
- Saves the combined plot as a PNG file.

### Execution Flow

1. **Define Colors and Maps**:
   - Define colors for various cell types, region types, and summary types.
   - Create maps for easy reference.

2. **Create Upset Plots**:
   - Generate Upset plots for QTL mapping, fine-mapping, colocalization, PTWAS, and INTACT results.
   - Save plots in PDF format for higher resolution.

3. **Combine Plots**:
   - Combine individual Upset plots into a single visualization.
   - Save combined plots in PNG format for easy sharing and integration.

### Usage Notes

- Ensure all necessary libraries are installed and properly loaded.
- Adjust color maps and other settings as needed.
- Provide accurate file paths to ensure correct data loading and plot saving.

---

## Script 03                                    
## 05_plot_pc_qtl.R

This script plots quantitative trait loci (QTL) for variants and methylation levels to visualize significant QTLs identified by principal components (PCs) but not by averages. It creates box plots to compare the distribution of methylation across genotypes for different summary types.

## Overview

- Defines functions to load QTL data, including significant QTLs, variant information, and methylation information.
- Generates box plots to visualize the distribution of methylation levels across genotypes for different summary types.
- Combines variant and methylation information to create comprehensive plots that help identify potential QTL interactions.
- Saves the plots in PNG and RDS formats for further analysis and sharing.

## Script Details

### Functions

#### `load_qtl_list()`

- Loads a pre-filtered list of significant QTLs, focusing on those identified by PCs but not by averages.

#### `select_qtl(qtl_list)`

- Selects a random QTL from the provided list for visualization.

#### `load_variant_info(qtl)`

- Loads genotype information for the selected QTL's variant.

#### `load_methylation_info(qtl)`

- Loads methylation data for the selected QTL's gene and summary type.

#### `plot_info(variant_info, meth_info, qtl)`

- Joins variant and methylation information.
- Creates a box plot to visualize the distribution of methylation levels across genotypes.
- Saves the plot in PNG format and as an RDS object for further analysis.

### Execution Flow

1. **Load QTL List**:
   - Loads a list of significant QTLs to select from.

2. **Select a QTL**:
   - Randomly selects a QTL to plot.

3. **Load Variant Information**:
   - Retrieves genotype data for the selected variant.

4. **Load Methylation Information**:
   - Retrieves methylation data for the selected gene.

5. **Plot QTL Information**:
   - Combines variant and methylation information.
   - Generates box plots to visualize the distribution of methylation levels across genotypes for different summary types.
   - Saves plots in PNG and RDS formats.

### Usage Notes

- Ensure necessary libraries are installed and loaded.
- Provide accurate file paths for data loading and plot saving.
- Adjust the number of iterations and seed for consistent results.
- Customize plot aesthetics and output formats as needed.

---

## Script 04                                    
## 06_plot_celltype_summary_figure.R

This script creates a manuscript figure showcasing cell-type and summary-type sections, visualizing data related to gene expression and cell-type proportions. It combines multiple plots into a cohesive figure for analysis and presentation in research papers.

## Overview

- Defines functions to load various figures and data related to cell-type proportions, UMAPs, gene summaries, and CpG-PC counts.
- Generates box plots, scatter plots, and violin plots to visualize different aspects of gene expression and methylation.
- Combines these plots into a single figure, which is saved in PNG format.

## Script Details

### Functions

#### `load_celltype_props()`

- Loads a pre-generated plot of cell-type proportions.

#### `load_celltype_umap()`

- Loads a pre-generated plot of a cell-type UMAP.

#### `combine_plots(celltype_props_fig, celltype_umap_fig, summarytype_fig, cpg_pc_fig)`

- Combines various plots into a single figure.
- Formats and standardizes plot themes for consistency.
- Combines cell-type proportions, UMAP, summary-type, and CpG-PC plots.
- Saves the combined plot in PNG format.

#### `load_summarized_genes_data(cell_type, region_type)`

- Loads data on summarized genes for specified cell types and region types.
- Extracts necessary information for plotting.

#### `get_st_genes(matched_dat)`

- Retrieves gene symbols and their annotations, including Alzheimer's disease (AD) gene association scores.

#### `make_summarytype_fig()`

- Creates a figure illustrating differences in methylation levels for specified genes across summary types.
- Generates violin plots and adds relevant annotations.
- Saves plots in PNG and RDS formats.

#### `make_cpg_pc_fig()`

- Creates a scatter plot of CpG counts against rPC counts.
- Fits a linear model and visualizes the relationship.
- Returns the scatter plot.

### Execution Flow

1. **Load Data**:
   - Loads cell-type proportions and UMAP plots.
   - Loads summarized gene data for specified cell types and region types.
   - Retrieves gene symbols and annotations.

2. **Generate Plots**:
   - Creates summary-type plots for specified genes.
   - Generates a scatter plot illustrating the relationship between CpG counts and rPC counts.

3. **Combine Plots**:
   - Combines cell-type proportions, UMAP, summary-type, and CpG-PC plots into a single figure.
   - Standardizes plot themes for consistency.
   - Saves the combined plot in PNG format.

### Usage Notes

- Ensure all necessary data files and plot objects are available.
- Adjust file paths, cell types, and region types as needed.
- Modify plot aesthetics and themes to fit manuscript requirements.
- Check for consistency in color schemes and legends across plots.

---

## Script 05                                    
## 07_check_eqtls.R

This script is designed to check for Alzheimer's disease (AD) genes across various stages of QTL analysis and assess their proximity to the APOE gene, a well-known marker for Alzheimer's risk.

## Overview

The script includes functions to:

- **Retrieve Total Features**: Get total genes tested from specific datasets.
- **Check for QTL Genes**: Identify significant genes at various stages of QTL analysis, including colocalization, PTWAS, and INTACT.
- **Calculate Hypergeometric Test**: Compute p-values to assess the significance of gene overlaps.
- **Load APOE CentiMorgan Map**: Get centiMorgan data for APOE and other regions of interest.
- **Load Region Annotations**: Retrieve annotations for specified gene regions.
- **Find Nearby Genes**: Find genes that are close to APOE and other loci.
- **Create and Save Plots**: Plot various data and save results for further analysis.

## Script Details

### Functions

#### `get_total_features(runname)`

- Returns a list of all unique features (genes) tested in a given run.

#### `check_qtls(old_runname, total_genes, ad_genes)`

- Checks significant QTLs after QTL mapping.
- Computes hypergeometric p-values to assess overlap between significant genes and AD-related genes.

#### `check_finemapped(runname, total_genes, ad_genes)`

- Checks significant QTLs after fine-mapping and computes overlap with AD genes.

#### `check_colocalization(runname, coloc_thresh, gene_symbols)`

- Checks for colocalized genes above a specified threshold.

#### `check_ptwas(runname, ptwas_thresh, gene_symbols)`

- Checks significant PTWAS genes using a specific threshold.

#### `check_intact(runname, intact_thresh, gene_symbols)`

- Checks significant INTACT genes based on a specific threshold.

#### `plot_counts(counts_df)`

- Plots counts of feature variants, QTL significance, and fine-mapped QTLs.

#### `get_pval(total_genes, ad_genes, test_genes)`

- Computes hypergeometric p-values for gene overlaps.

#### `load_centimorgan_map()`

- Loads centiMorgan map data and identifies regions of interest.

#### `load_region_annots(region_type)`

- Loads gene annotations for a specified region type.

#### `find_nearby_genes(centi_map, region_annots)`

- Finds genes located near specified centiMorgan regions.

### Execution Flow

1. **Load Data**:
   - Loads gene symbols and other necessary data for the analysis.

2. **Check Genes at Different Stages**:
   - Checks for significant genes across different stages of QTL analysis (colocalization, PTWAS, INTACT).
   - Computes hypergeometric p-values to assess overlap with AD-related genes.

3. **Analyze Proximity to APOE**:
   - Loads centiMorgan data and finds genes near APOE and other regions of interest.

4. **Save and Plot Results**:
   - Saves gene data and plots results in various formats for further analysis.

### Usage Notes

- Ensure the necessary data files are available, including gene symbols, centiMorgan maps, and QTL results.
- Adjust thresholds and other parameters as needed.
- The script includes sections for detailed statistical analysis and proximity checks for specific genes.
- Consider additional validation methods to confirm significant results.

---

## Script 06                                    
## 08_plot_dm_figure.R

This script is designed to combine multiple plots into a single comprehensive figure for a manuscript. It involves the following key elements:

## Overview

- **Data Loading**: Load different plots saved as `.rds` files.
- **Plot Formatting**: Standardize the formatting of the plots to ensure consistency across the combined figure.
- **Plot Combination**: Combine the plots into a structured arrangement using `ggpubr` and `gridExtra`.

## Script Details

### Functions

#### `combine_plots(plot1, plot2, plot3)`

- Combines three plots into a single figure.
- Applies common formatting to ensure consistency across plots.
- Uses `ggarrange` to organize plots in a specified layout.

#### `main()`

- Loads the individual plots from their respective directories.
- Invokes `combine_plots` to create the final figure.
- Saves the combined plot as a `.png` file in the specified output directory.

### Execution Flow

1. **Load Plots**:
   - Loads three different plots from specified locations. These plots represent various aspects of a larger analysis.
   - The plots are in `.rds` format, allowing easy retrieval of `ggplot` objects.

2. **Format Plots**:
   - Standardizes the appearance of the plots to ensure a consistent look in the combined figure.
   - Customizes aspects like text size, legend position, and color scale.

3. **Combine Plots**:
   - Combines the plots into a multi-panel figure.
   - Uses `ggarrange` to create a top row with two plots and a second row with one plot.
   - This arrangement allows for a clear and structured presentation of multiple plots.

4. **Save Combined Plot**:
   - Saves the final combined plot as a `.png` file in the specified output directory.
   - Adjusts plot dimensions to ensure the entire figure is visible without distortion.

### Usage Notes

- Ensure the necessary `.rds` files are available in the specified directories.
- Adjust plot formatting and layout as needed.
- Consider additional elements for inclusion in the combined plot, depending on the manuscript's requirements.
- The script uses `ggpubr` and `gridExtra` for plot combination; ensure these packages are installed.
- Modify the `savedir` variable to specify the desired output location for the final combined plot.

---

## Script 07                                    
## 10_plot_supp_figures.R

This script generates and formats supplementary figures for a manuscript, focusing on loading, formatting, and combining plots. It follows these steps:

## Script Overview

- **Data Loading**: Load individual plots from saved `.rds` files.
- **Plot Formatting**: Customize the appearance of the plots.
- **Plot Combination**: Use `ggpubr` to combine plots into a single supplementary figure.

## Script Breakdown

### Functions

#### `plot_supp4()`

- Loads individual plots from specific directories.
- Applies additional formatting to ensure a consistent appearance.
- Combines multiple plots into one figure.
- Saves the combined figure to a specified output directory.

#### `main()`

- Invokes `plot_supp4` to create and save the supplementary figure.

### Detailed Steps

1. **Load Plots**:
   - Load two plots from the given paths. These plots represent gene counts and feature proportions comparisons.
   - They are loaded as `ggplot` objects from `.rds` files.

2. **Format Plots**:
   - Adjust text angle and legend positioning to improve readability.
   - Customize the legend spacing and axis text to ensure a clear presentation.

3. **Combine Plots**:
   - Use `ggarrange` to combine the formatted plots into a single figure.
   - Arrange the plots horizontally with shared legends for a cohesive look.
   - Label the plots for easy reference.

4. **Save Combined Plot**:
   - Save the combined plot to the specified output directory as a `.png` file.
   - Adjust the plot size to ensure proper visibility and clarity.

### Usage Notes

- Ensure that the source `.rds` files are in the specified directories.
- Adjust plot formatting as needed to match the manuscript's style and requirements.
- Modify the `savedir` variable to specify the desired output location for the final supplementary figure.
- The script uses `ggpubr` for plot combination; ensure this package is installed before running the script.