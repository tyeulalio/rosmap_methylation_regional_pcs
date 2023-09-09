# Methylation Preprocessing Pipeline

This pipeline contains a set of functions designed to process and perform quality control on methylation data.

## Table of Contents

- [Usage](#usage)
  - [load_samples](#load_samples)
  - [qualitycontrol_samples](#qualitycontrol_samples)
  - [preprocess_pipeline](#preprocess_pipeline)
  - [bmiq_norm](#bmiq_norm)
  - [filter_snps](#filter_snps)
  - [qc_probes_norm_remove_reps](#qc_probes_norm_remove_reps)
  - [get_filtered_rgset](#get_filtered_rgset)
- [Running the analysis](#running-the-analysis)


## Usage

Below are descriptions and usage guidelines for each of the functions.

### `load_samples`

#### Description

1. Reads sample metadata from a specified file.
2. Filters samples based on data subtype and file format.
3. Formats annotations for further processing.

#### Parameters

None.

#### Returns

`formatted_annotations`: a dataframe containing filtered and formatted filenames


---

### `qualitycontrol_samples`

#### Description

1. Filters samples based on bisulfite conversion efficiency.
2. Uses control metrics from `ewastools` to remove failed samples.
3. Saves the filtered datasets as RDS files for future use.

#### Parameters

- `targets`: A data frame containing the target sample IDs.
- `rgset`: An rgset object containing the methylation data.

#### Returns

`res`: a list containing the quality control-filtered methylation and target sample data. List elements include:

- `rgset_qcfilt``: QC'ed rgset methylation data
- `targets_qcfile`: QC'ed targets sample data

---

### `preprocess_pipeline`

#### Description

1. Calls the `qualitycontrol_samples` function to perform quality control on sample data.
2. Sets the newly filtered results for further analysis.
3. Calls the `qc_probes_norm_remove_reps` function to perform probe quality control normalization, and removal of replicates.

#### Parameters

- `targets`: A data frame containing metadata and annotations for the samples.
- `rgset`: A `GenomicRatioSet` object that contains the raw methylation data.

#### Returns

None.

---


### `bmiq_norm`

#### Description

1. Performs BMIQ normalization on the beta-values of the methylation data.
2. Saves the normalized beta-values and target data as RDS files.

#### Parameters

- `targets_champ`: A data frame containing the target sample IDs.
- `betas_champ`: A matrix containing the beta-values for methylation data.

#### Returns

None.

---


### `filter_snps`

#### Description

1. Maps methylation probes to their genomic positions.
2. Reads ROSMAP annotated VCF files to gather Allele Frequencies (AF) for SNPs across multiple chromosomes.
3. Filters out SNPs with an AF greater than a certain threshold.
4. Filters the beta-values matrix to remove probes overlapping with such SNPs.
5. Saves the filtered datasets and probe position mappings as RDS files for future use.

#### Parameters

- `betas`: A matrix containing the beta-values for methylation data.
- `targets`: A data frame containing the target sample IDs.

#### Returns

The function does not return any value but saves the following datasets:

- `probe_position_map_hg37.rds`: An RDS file containing the probe to position mappings, saved in the specified directory.
- `remove_snps_af_<thresh>.rds`: An RDS file containing the SNPs to be removed, saved in the specified directory.
- `snp<thresh>_window0_filtered_betas.rds`: An RDS file containing the filtered beta-values, saved in the specified directory.
- `snp<thresh>_window0_filtered_targets.rds`: An RDS file containing the filtered target sample IDs, saved in the specified directory.

#### Usage Messages

- Prints messages indicating the processing status for each chromosome.
- Prints a summary of allele frequencies.
- Prints the number of SNPs to be removed.

---

### `qc_probes_norm_remove_reps`

#### Description

1. Performs quality control, normalization, and replicate removal on Illumina 450K DNA methylation data.
2. Incorporates noob dye bias correction, sex chromosome filtering, SNP filtering, cross-reactive probe exclusion, and BMIQ normalization.
3. Saves multiple stages of the processed data as RDS files.

#### Parameters

- `targets`: A data frame containing sample metadata such as sample IDs.
- `rgset`: An object holding the red and green channel methylation data, commonly of class `RGChannelSet`.

#### Returns

None. The function saves multiple `.rds` files at various stages of the processing:

1. `noob_preprocessed_meth.rds`: After noob dye bias correction.
2. `snp_filtered_thresh{thresh}.rds`: After SNP filtering.
3. `detP_ewas_raw.rds`: Dataset with raw detection p-values.
4. `betas_filtered.rds`: Final filtered beta values.
5. `targets_filtered.rds`: Final filtered targets.
6. `betas_bmiq.rds`: After BMIQ normalization.

---

### `get_filtered_rgset`

#### Description

1. Loads a pre-existing RGChannelSet object from an RDS file which contains red and green channel methylation data.
2. Reads filtered beta-values from another RDS file.
3. Maps the CpG IDs to CpG positions to align the RGChannelSet with the filtered beta-values.
4. Filters the RGChannelSet to only include the final set of CpGs that were included in the filtered beta-values.
5. Estimates cell counts using the `estimateCellCounts` method.
6. Saves the estimated cell counts as an RDS file.

#### Parameters

None. The function relies on global variables like `savedir` for specifying directory paths, and `thresh` for specifying the SNP threshold, among other things. These should be properly set before running the function.

#### Returns

None. The function saves multiple `.rds` files:

1. `neun_estimated_cell_counts.rds`: After estimating the cell counts.

Note: This function uses several saved RDS files like `rgset_qcsfilt_hg19.rds` and `snp{thresh}_window0_filtered_betas.rds`, which should exist in the specified directory before running the function. The names and directories of these files are specified by the `savedir` and `thresh` variables.

#### Side Effects

The function modifies the global state by reading and saving various `.rds` files and also by altering the state of variables like `rgset`, `formatted_betas`, `filtered_rgset`, etc. Make sure to understand the dependencies and side effects before running this function.


## Running the analysis

### Description

This part of the code performs the following tasks:

1. Retrieves annotation data for the Illumina Human Methylation 450k array, which is essential for mapping probes to the genome.
2. Loads the sample sheet that contains phenotype and array information.
3. Reads in raw IDAT files containing the methylation data. This step can take a long time.
4. Assigns column names to the raw RGChannelSet (`rgset_raw`) based on the sample sheet.
5. Calls a preprocessing pipeline that performs a comprehensive analysis including sample quality control, probe quality control, data normalization, and removal of replicate samples.
