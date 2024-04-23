library(data.table)
library(GenomicRanges)
library(liftOver)
library(ggplot2)
library(tidyverse)

# check that the annotations for gwas and qtl SNPs matched
# look for the main() function at the bottom to see the flow of this script

# save this mapping from chromosome position to the var name used in the QTL analysis
savedir <- paste0("/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/")
dir.create(savedir)


gwas_types <- c('wightman', 'kunkle', 'trubetskoy')
gwas_type <- gwas_types[1]
savedir <- paste0("/path/to/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/", gwas_type, "/")
dir.create(savedir)

# check on the variants that failed to map
check_failed <- function(){
    # read in the rejected records
    datafile <- paste0("../../../output/06_fine_mapping_colocalization/02_gwas_prep/01_liftover_gwas/rejected_records.vcf")
    rejected <- read_tsv(datafile, comment='##')
    head(rejected)

    chainfile <- paste0("/home/eulalio/deconvolution/data/hg19ToHg38.over.chain")
    chain <- import.chain(chainfile)

    # format the data
    formatted_rejected <- rejected %>%
        rename(chrom="#CHROM",
               pos=POS,
               name=ID,
               ref=REF,
               alt=ALT,
               qual=QUAL,
               filter=FILTER,
               info=INFO
        ) %>%
        mutate(start=as.numeric(pos), end=start+1)
    head(formatted_rejected)

    # create genomic ranges 
    rejected_gr <- makeGRangesFromDataFrame(formatted_rejected, keep.extra.columns=TRUE)
    rejected_gr

    lifted <- liftOver(rejected_gr, chain) %>%
        as.data.frame()
    head(lifted)

}



load_qtl <- function(){
    # load qtl snps
    datafile <- paste0("../../../output/05_qtl_analysis/09_remove_outliers/all_chroms_hg38.bim")
    qtl <- read_tsv(datafile, col_names=c('chrom', 'var', 'group', 'pos', 'ref', 'alt'))
    head(qtl)

    # select columns that we care about
    qtl_snps <- qtl %>%
        select(chrom, var, pos, ref, alt)
    head(qtl_snps)
    dim(qtl_snps)

    qtl_snps
}


join_gwas_qtls <- function(gwas, qtl_snps, gwas_type){
    # join the two together on chromosome and position
    head(qtl_snps)
    head(gwas)
    
    qtl_gwas_snps <- qtl_snps %>%
        inner_join(gwas, by=c('chrom', 'pos'), suffix=c('.qtl', '.gwas')) 
    head(qtl_gwas_snps)


    # how many did we match?
    head(qtl_gwas_snps)
    print(paste("qtl:", nrow(qtl_snps),
                "gwas:", nrow(gwas),
                "matched:", nrow(qtl_gwas_snps)))
    dim(qtl_snps)
    dim(gwas)
    dim(qtl_gwas_snps)


    # save this mapping from chromosome position to the var name used in the QTL analysis
    savefile <- paste0(savedir, "chrom_pos_mapped_to_variant_name.tsv")
    savefile
    saveRDS(qtl_gwas_snps, savefile)

    # find positions with swapped ref/alt alleles - need to adjust gwas beta for these
    head(qtl_gwas_snps)
    swapped <- qtl_gwas_snps %>%
        separate(var, c('var_chrom', 'var_pos', 'var_ref', 'var_alt'), sep='_', remove=FALSE) %>%
        filter(var_ref == alt.gwas,
               var_alt == ref.gwas
        )
    head(swapped)
    dim(swapped)
    print(paste("swapped:", nrow(swapped)))


    # find variants with either ref or alt not matching - need to remove these
    head(qtl_gwas_snps)
    mismatched <- qtl_gwas_snps %>%
        filter(!gwas_name %in% swapped$gwas_name) %>%
        separate(var, c('var_chrom', 'var_pos', 'var_ref', 'var_alt'), sep='_', remove=FALSE) %>%
        filter((var_ref != ref.gwas) | (var_alt != alt.gwas)) 
    head(mismatched)
    dim(mismatched)
    print(paste("mismatched:", nrow(mismatched)))

    # check what the mismatched look like
    mismatched %>%
        select(chrom, pos, var_ref, ref.gwas, var_alt, alt.gwas)


    # remove the mismatched from final set
    dim(qtl_gwas_snps)
    head(qtl_gwas_snps)
    final_matched <- qtl_gwas_snps %>%
        filter(!var %in% mismatched$var)
    head(final_matched)
    dim(final_matched)

    # save this mapping from chromosome position to the var name used in the QTL analysis
    savefile <- paste0(savedir, "matched_chrom_pos_mapped_to_variant_name.tsv")
    savefile
    saveRDS(final_matched, savefile)
    
    print(paste("final matched", nrow(final_matched)))

    #final_matched <- readRDS(savefile)
    final_matched
}


## ----- CREATE Z-SCORES -----

load_gwas_scores <- function(final_matched){
    # load the hg37 gwas data with scores
    if (gwas_type == 'wightman'){
        datafile <- paste0("/path/to/data/ad_gwas/wightman_etal_2021_results.txt.gz")
        gwas_scores_hg37 <- read_tsv(datafile) 
        head(gwas_scores_hg37)
        head(final_matched)

        # filter the gwas scores using the final matched snps
        filtered_scores_hg37 <- gwas_scores_hg37 %>%
            mutate(Chromosome = paste0("chr", chr)) %>%
            unite(gwas_name, Chromosome, PosGRCh37, testedAllele, otherAllele, sep='_') %>%
            inner_join(final_matched)
        head(filtered_scores_hg37)

        dim(final_matched)
        dim(filtered_scores_hg37)

        # check if these are in kunkle
        no_inf <- filtered_scores_hg37 %>%
            filter(!z %in% c(Inf, -Inf)) 
        head(no_inf)

        (max_z <- max(no_inf$z))

        summary(no_inf$z)
        filtered_scores_hg37 %>%
            filter(var == 'chr7_100192194_T_C')

        head(filtered_scores_hg37)
        # just need var id and zscore
        final_scores <- filtered_scores_hg37 %>%
            select(var, z) %>%
            mutate(fixed_z = ifelse(z == Inf, max_z, z)) %>%
            mutate(fixed_z = ifelse(z == -Inf, -max_z, fixed_z)) %>%
            select(var, zscore=fixed_z)
        head(final_scores)

        savefile <- paste0(savedir, "gwas_updated_scores.csv")
        write_csv(final_scores, savefile)

        summary(final_scores$zscore)
        1

    }
    dim(final_scores)
    dim(final_matched)
    head(final_scores)


    # now remove multi-SNP locations
    # select the one with best z-score
    #ld_map <- read_tsv(paste0(savedir, "gwas_snp_ld_block_annots.tsv"))
    #head(ld_map)

    #ld_scores <- ld_map %>%
        #select(ld_block=phenotype_id, var=variant_id) %>%
        #left_join(final_scores)
    #head(ld_scores)

    #sub_ld_scores <- ld_scores %>%
        #filter(ld_block == 'ld_8')
    #head(sub_ld_scores)
    #dim(sub_ld_scores)
    #length(unique(sub_ld_scores$var))
    head(final_scores)
    multi <- final_scores %>%
        count(var) %>%
        arrange(-n)
    head(multi)

    # select highest score
    head(final_scores)
    filt_scores <- final_scores %>%
        group_by(var) %>%
        arrange(-(abs(zscore))) 
    filt_scores %>%
        filter(var == multi[1,]$var)

    sub_scores <- filt_scores %>%
        group_by(var) %>%
        slice(1)
    head(sub_scores)
    filt_multi <- filt_scores %>%
        count(var) %>%
        arrange(-n)
    head(filt_multi)

    dim(sub_scores)

    (savefile <- paste0(savedir, "formatted_signed_dgap_gwas_scores.txt"))
    write_delim(sub_scores, savefile, delim=" ", col_names=FALSE)
    1


    sub_scores
}



check_zscores <- function(){
    # format dgap input
    # need 2 columns:
    # - snp name
    # - unsigned zscore
    head(updated_scores)

    plotdir <- paste0(savedir, "plots/")
    dir.create(plotdir)

    # plot pval density, then plot zscore density to make sure we're not messing anything up
    head(updated_scores)
    dim(updated_scores)
    formatted_scores <- updated_scores %>%
        mutate(bins = cut(gwas_pval, breaks=10))
    head(formatted_scores)

    table(formatted_scores$bins)
    table(formatted_scores$gwas_pval < 0.05)


    # compute zscore using 2 tails
    dgap_scores <- updated_scores %>%
        unite('var_id', chrom:alt, sep='_') %>%
        select(var_id, gwas_pval, updated_beta, gwas_SE) %>%
        mutate(gwas_zscore = qnorm(gwas_pval/2, lower.tail=FALSE)) %>% # compute zscore
        mutate(pval = 2 * pnorm(q=gwas_zscore, lower.tail=FALSE)) %>% # make sure zscore is correct
        mutate(zscore2 = updated_beta / gwas_SE) %>%
        mutate(pval2 = 2 * pnorm(q=zscore2, lower.tail=FALSE)) %>% # make sure zscore is correct
        mutate(gwas_zscore = gwas_zscore * sign(updated_beta)) 
    head(dgap_scores, 10)
    dim(dgap_scores)


    formatted_scores <- dgap_scores %>%
        mutate(bins = cut(abs(gwas_zscore), breaks=10))
    head(formatted_scores)

    table(formatted_scores$bins)
    table(abs(formatted_scores$zscore2) > 1.96)
    table(abs(formatted_scores$gwas_zscore) > 1.96)


    # comparing two ways of computing zscore
    dgap_scores %>%
        mutate(diff = abs(gwas_zscore - zscore2)) %>%
        arrange(-diff)

    summary(dgap_scores$gwas_zscore - dgap_scores$zscore2)

    # no Inf values
    summary(dgap_scores$gwas_zscore)

    summary(dgap_scores$zscore2)

    final_scores = dgap_scores %>%
        select(var_id, gwas_zscore)
        #select(var_id, zscore2)

    savefile <- paste0(savedir, "formatted_signed_dgap_gwas_scores.txt")
    write_delim(final_scores, savefile, delim=" ", col_names=FALSE)
    1
}


create_ld_annotations <- function(){
    ### ---- create LD matrix annotations for snps --- ##

    # get the snps from the zscores files
    savefile <- paste0(savedir, "formatted_signed_dgap_gwas_scores.txt")
    zscores <- read_table(savefile, col_names=c('var', 'zscore'))
    head(zscores)


    # searate variant into postion
    formatted_snps <- zscores %>%
        select(var) %>%
        unique() %>%
        separate(var, c('chrom', 'start', 'ref', 'alt'), sep='_') %>%
        mutate(end = as.numeric(start)+1)
    head(formatted_snps)

    # load in the LD bloack
    datafile <- paste0("/path/to/data/ld_blocks/deCODE_EUR_LD_blocks.bed")
    ld <- read_table(datafile)
    head(ld)

    formatted_ld <- ld %>%
        mutate(block = paste0("ld_", row_number()))
    head(formatted_ld)

    # use genomic ranges to overlap snps with blocks
    snps_gr <- makeGRangesFromDataFrame(formatted_snps, keep.extra=TRUE)
    ld_gr <- makeGRangesFromDataFrame(formatted_ld, keep.extra=TRUE)

    overlaps <- findOverlaps(subject=snps_gr, query=ld_gr) %>%
        as.data.frame()
    head(overlaps)

    ld_snps <- cbind(snps=as.data.frame(snps_gr)[overlaps$subjectHits,], ld=as.data.frame(ld_gr)[overlaps$queryHits,])
    head(ld_snps)

    length(unique(ld_snps$ld.block))
    counts <- ld_snps %>%
        count(ld.block)
    summary(counts$n)

    # format the output
    formatted_ld_snps <- ld_snps %>%
        mutate(variant_id = paste(snps.seqnames, snps.start, snps.ref, snps.alt, sep='_')) %>%
        rename(phenotype_id = ld.block)
    head(formatted_ld_snps)

    savefile <- paste0(savedir, "gwas_snp_ld_block_annots.tsv")
    savefile
    write_tsv(formatted_ld_snps, savefile)
    1
}


load_gwas_snps <- function(gwas_type){
    # load gwas data
    print(paste("Loading gwas", gwas_type))
    datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/01_liftover_gwas/wightman_etal_2021_results_hg38.vcf")

    gwas <- read_tsv(datafile, comment='##') %>%
            select(chrom="#CHROM", pos=POS, gwas_name=ID, ref=REF, alt=ALT)
    head(gwas)
    length(unique(gwas$gwas_name))

    gwas
}

main <- function(){
    # match gwas, qtl snps
    # create ld gwas annotations

    # load gwas data
    gwas_type
    gwas <- load_gwas_snps(gwas_type)

    # load qtl snps
    qtl_snps <- load_qtl()

    # join qtls and gwas 
    final_matched <- join_gwas_qtls(gwas, qtl_snps, gwas_type)

    final_scores <- load_gwas_scores(final_matched)

    # 
    create_ld_annotations()
}

1
