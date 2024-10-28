# rosmap_methylation_regional_pcs

# regionalpcs improve discovery of DNA methylation associations with complex traits

[![DOI](https://zenodo.org/badge/638786676.svg)](https://doi.org/10.5281/zenodo.14004153)

## Overview
This release provides the code and resources used in the Nature Communications paper, "regionalpcs improve discovery of DNA methylation associations with complex traits." The repository includes functions for simulations, methylation data preprocessing, cell type deconvolution, regional principal components (rPCs) generation, differential methylation analysis, and methylation QTL mapping, as well as code to integrate methylation data with GWAS for causal inference.

## Key Features:
- RegionalPCS Method: Summarize methylation at the gene-level using PCA, including functions to select optimal principal components.
- Differential Methylation Analysis: Identify differentially methylated regions using rPCs, averages, or individual CpGs.
- meQTL Mapping: Map methylation quantitative trait loci (meQTL) using regional summaries and fine-mapping for causal variant identification.
- Integration with GWAS: Link methylation changes to disease risk using colocalization and instrumental variable analysis.
  
## Getting Started:
Follow the instructions in README.md wihtin each folder for main steps of each script.

## Citing regionalpcs

[to do]

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
