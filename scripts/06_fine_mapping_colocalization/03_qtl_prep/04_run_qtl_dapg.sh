#!/usr/bin/env bash

# run the deterministic approximation of posteriors (DAP)
# to get pobabilistic QTL annotations for fastENLOC
# following GTEx instructions here: https://github.com/xqwen/dap/tree/master/gtex_v8_analysis

# -d sbam file
# -p prior file from torus
# --all forces output of all SNPs, not just noteworthy ones
# -t number of parallel threads
# -ld_control lowest LD threshold (R^2) to admit a SNP into a signal cluster


# variables to set
SUMMARY_TYPE=$2
CELL_TYPE=$3
REGION_TYPE=$4

PROPROTION_TYPE="with_proportions"

# create the string for this run
RUNSTR="${SUMMARY_TYPE}_${CELL_TYPE}_${REGION_TYPE}"

# get file with genes to process
GENESDIR="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/03_dapg_input/${PROPROTION_TYPE}/${RUNSTR}"
GENESFILE="${GENESDIR}/filtered_genes_list.txt"

# get gene number from command line
i=$1

# used this command to filter out genes with missing priors
if [ $i == 1 ]; then
    echo "creating filtered genes list file"
    comm -23 <(sort ${GENESDIR}/genes_list.txt) <(sort ${GENESDIR}/missing_priors_file.txt) > $GENESFILE
fi


# -- input and output directory
# create the output directory
OUTDIR="${HOME}/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/03_qtl_prep/04_dapg" 
mkdir -p $OUTDIR

# create directory specific to te run
OUTDIR="${OUTDIR}/${PROPROTION_TYPE}"
mkdir -p $OUTDIR
OUTDIR="${OUTDIR}/${RUNSTR}"
mkdir -p $OUTDIR


# get the gene from the file
GENE=$(sed "${i}q;d" $GENESFILE)

echo "$GENE $i"

# format gene properly
FORMATTED_GENE=$(echo "$GENE" | tr \\. - )

OUTFILE="${OUTDIR}/${FORMATTED_GENE}.dapg" 

# skip ifoutput file exists
if [ -f "$OUTFILE" ]; then
    echo "output file exists"
    exit
fi

# run dap-g
dap-g -d "${GENESDIR}/sbams/${GENE}_fastqtl_singl_snp_output.sbam" \
    -p "${GENESDIR}/priors/${GENE}_priors.txt" \
    -ld_control 0.25 \
    --all \
    -t 4 \
    -n ${GENE} \
    > $OUTFILE
