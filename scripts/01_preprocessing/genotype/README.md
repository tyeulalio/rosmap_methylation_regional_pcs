# Processing Genotypes

## Table of Contents

- [VCF LiftOver Pipeline](#vcf-liftover-pipeline)
- [Usage](#usage)
- [Workflow](#workflow)
- [Parameters](#parameters)
- [Contributing](#contributing)
- [License](#license)

## VCF LiftOver Pipeline

This pipeline provides a step-by-step process to lift over VCF files from `hg19` to `hg38` using GATK's LiftoverVcf tool.

### Data Source

The data was sourced from the [ROMSAP harmonization whole genome sequencing as VCF files](https://www.synapse.org/#!Synapse:syn11707420). To download these files:

1. Navigate to the provided link.

2. Click on "Download Options" to find detailed instructions.

Note: Accessing this data requires a signed Data Use Agreement.

### Prerequisites

1. **Modules**: Before using the script, ensure the following modules are installed and available on your system:
    - `gatk/4.0.10.0`
    - `bcftools`

2. **Files & Directories**:
    - Source VCF files situated in the designated paths.
    - Reference fasta files for `hg19` and `hg38` must be correctly pointed in the script.
    - Chain file for liftOver, moving from `hg19` to `hg38`.

### Usage

1. Clone this repository to your local machine:

```bash
git clone [repository-url]
cd [repository-dir]
```

### Run the script:

```bash
./liftover_script.sh
```

### Workflow

- **Concatenation**: Initially, chromosomes from the ROMSAP VCFs are concatenated. The code prepends "chr" in front of each chromosome to match the annotations in the LiftOver chain file.
- **Normalization**: The VCF undergoes normalization using the `bcftools norm` function. This step uses the `hg19 GATK resource bundle` for accuracy, retains multiallelic sites together (avoiding issues later in the pipeline), and swaps ref/alt alleles if incorrectly placed.
- **LiftOver**: This step employs GATK/4.0.10.0's `LiftoverVcf` tool to transform the VCF from `hg19` to `hg38`. The tool uses the `hg19tohg38` chain and references the `hg38` fasta genome file located in the SCG reference folder. It's critical to adjust the VCF header to include "chr" before chromosome definitions to ensure compatibility. The liftOver operation occasionally faces "out of memory" issues. To mitigate this, set `MAX_RECORDS_IN_RAM` to 500. If the problem persists, rerunning the script typically resolves it.
- **Final Normalization**: A final round of normalization is performed post LiftOver.

### Parameters

Inside the script, specific parameters are set at the beginning. These include paths to the required VCF files, reference genomes, and the chain file for LiftOver. To ensure the script works correctly in your environment, adjust these parameters as necessary.

### Additional Details

- The main code for this pipeline is found in `01_liftover_vcf.sh`.
- For reference sequences, the pipeline utilizes `/reference/RefGenomes/GATK_Resource_Bundle/hg38/hg38.fa`.
- If you encounter memory issues, rerunning the script is usually a quick solution, though it's not guaranteed to always work.

## Contributing

Should you wish to contribute, whether it's through bug reports, feature requests, or direct code contributions, you're welcome to do so! Simply open an issue or send in a pull request.

## License
[Add licensing details here]
