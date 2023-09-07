#!/bin/bash

# using plink on scg
module load plink
module load plink2

# run plink to generate files 
# need to run all three parts here, some are commented out during testing


# try to merge files, this will throw an error because of multiallelic sites
# use the output of this to remove those sites
echo "-- Finding multiallelic sites --"
plink \
    --make-bed \
    --biallelic-only \
    --merge-list "../../output/05_qtl_analysis/05_filter_mishap/${CHROMDIR}bim_files_hg38.txt" \
    --exclude "../../output/05_qtl_analysis/05_filter_mishap/${CHROMDIR}exclude_variants_hg38.txt" \
    --keep "../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38"

# excluding multiallelic sites found
for CHROM in {1..22}
do
    echo "removing multiallelic sites for chromosome ${CHROM}"
        plink2 \
            --pfile "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_chr${CHROM}" \
            --make-pgen \
            --exclude "../../output/05_qtl_analysis/06_merge_bims/all_chroms_hg38-merge.missnp" \
            --out "../../output/05_qtl_analysis/06_merge_bims/rosmap_wgs_hg38_chr${CHROM}_refiltered"
done

# now merge the files
echo "-- Merging files --"
plink2 \
    --make-pgen \
    --pfile "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms" \
    --exclude "../../output/05_qtl_analysis/05_filter_mishap/${CHROMDIR}exclude_variants_hg38.txt" \
    --keep "../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38"

plink2 \
    --make-bed \
    --pfile "../../data/rosmap_wgs_harmonization_plink/rosmap_wgs_hg38_all_chroms" \
    --exclude "../../output/05_qtl_analysis/05_filter_mishap/${CHROMDIR}exclude_variants_hg38.txt" \
    --keep "../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38"
    #--pmerge-list "../../output/05_qtl_analysis/05_filter_mishap/${CHROMDIR}bim_files_hg38_refiltered.txt" \

# create counts of reference allele
echo "-- Creating counts file"
plink2 \
    --pfile "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_ref_counts_hg38" \
    --export A-transpose
    

#create ped file
echo "-- Creating ped file"
plink2 \
    --pfile "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_hg38" \
    --out "../../output/05_qtl_analysis/06_merge_bims/${CHROMDIR}all_chroms_ref_counts_hg38" \
    --export ped
