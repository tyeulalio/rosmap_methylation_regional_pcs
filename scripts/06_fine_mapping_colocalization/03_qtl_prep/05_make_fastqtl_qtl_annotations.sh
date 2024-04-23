
# use fastqtl to create annoations from the fine-mapped QTL data
# following the tutorial on the github https://github.com/xqwen/fastenloc/tree/master/tutorial/
# formats the QTL DAP-G output for FastQTL
# the summarize_dap2enlco.pl script come from the FastEnloc github
# this script requires a SNP annotation file in VCF form - created by script 05/02 using plink

# arguments
# dir = directory with DAP fine-mapped qtl results
# vcf = file to annotate all SNP's positions

SUMMARYTYPE=$1
CELLTYPE=$2
REGIONTYPE=$3

PROPORTION_TYPE="with_proportions"

RUNNAME="${SUMMARYTYPE}_${CELLTYPE}_${REGIONTYPE}"

# create output directory
OUTDIR="/path/to/output/06_fine_mapping_colocalization/03_qtl_prep/05_fastqtl_qtl_annotations"
mkdir -p ${OUTDIR}

OUTDIR="${OUTDIR}/${PROPORTION_TYPE}" 
mkdir -p ${OUTDIR}
OUTDIR="${OUTDIR}/${RUNNAME}" 
mkdir -p ${OUTDIR}

DAPDIR="/path/to/output/06_fine_mapping_colocalization/03_qtl_prep/04_dapg/${PROPORTION_TYPE}/${RUNNAME}"

time /path/to/programs/fastenloc/src/summarize_dap2enloc.pl \
    -dir $DAPDIR \
    -vcf "/path/to/data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms.vcf.gz" \
    | gzip - > "${OUTDIR}/fastenloc.qtl.annotation.vcf.gz"
