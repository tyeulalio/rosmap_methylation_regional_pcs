library(coloc)
library(tidyverse)

# format the gwas data

savedir <- paste0("/path/to/new_rosmap/output/07_ptwas/04_formatted_gwas/")
dir.create(savedir)

# select the gwas to use
gwas_type <- "wightman"

if (gwas_type == 'wightman'){
    # get other data fields from gwas
    savefile <- paste0(savedir, gwas_type, "_joined_full_gwas.rds")
    if (!file.exists(savefile)){
        gwas_file <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/",
                            gwas_type,
                            "/formatted_signed_dgap_gwas_scores.txt")
        gwas <- read_table(gwas_file, col_names=c('var', 'zscore'))
        head(gwas)

        # get map to original variant names
        gwas_file <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/",
                            gwas_type, "/matched_chrom_pos_mapped_to_variant_name.tsv")
        gwas_map <- readRDS(gwas_file)
        head(gwas_map)

        # read in original data
        datafile <- paste0("/path/to/data/ad_gwas/wightman_etal_2021_results.txt.gz")
        full_gwas <- read_table(datafile)
        head(full_gwas)

        # join the additional columns to zscore
        joined_gwas <- gwas %>%
            left_join(gwas_map)
        head(joined_gwas)
        dim(gwas)
        dim(joined_gwas)

        head(joined_gwas)
        head(full_gwas)
        
        # now join full gwas
        formatted_joined <- joined_gwas %>%
            separate(gwas_name, c('chr', 'PosGRCh37', 'testedAllele', 'otherAllele'), sep='_') 
        head(formatted_joined)

        joined_full_gwas <- formatted_joined %>%
            mutate(chr = str_remove(chr, 'chr')) %>%
            mutate(chr = as.numeric(chr)) %>%
            mutate(PosGRCh37 = as.numeric(PosGRCh37)) %>%
            left_join(full_gwas)
        head(joined_full_gwas)

        saveRDS(joined_full_gwas, savefile)
    } else{
        joined_full_gwas <- readRDS(savefile)
    }

    # need to get MAFs for these variants
    datafile <- paste0("/path/to/new_rosmap/data/rosmap_wgs_liftover_hg38/MAFs/rosmap_wgs_hg38_all_chroms_MAFs.txt")
    mafs <- read_table(datafile, col_names=c('chrom', 'pos', 'ref', 'alt', 'maf'))
    head(mafs)

    head(joined_full_gwas)

    joined_gwas_mafs <- joined_full_gwas %>%
        mutate(ref=ref.gwas,
               alt=alt.gwas
        ) %>%
        inner_join(mafs) 
    head(joined_gwas_mafs)

    dim(joined_full_gwas)
    dim(joined_gwas_mafs)

    # try this method to recover effect size and SE from the zscore and p-value
    get_stats <- function(){
        # get SE from zscore and pvalue
        # get beta from zscore and se
        # https://www.biostars.org/p/319584/#414335
        head(joined_gwas_mafs)
        stats_df <- joined_gwas_mafs %>%
            mutate(denom = sqrt(2*maf*(1-maf)*(N + zscore^2))) %>%
            mutate(SE = 1 / denom) %>%
            mutate(beta = zscore / denom)
        head(stats_df)

        sub_stats <- stats_df %>%
            select(var, zscore, p, SE, beta, N) %>%
            mutate(new_zscore=beta/SE)
        head(sub_stats)

        savefile <- paste0(savedir, gwas_type, "_estimated_effect_sizes.tsv")

        sub_stats <- read_tsv(savefile)
        head(sub_stats)

        # make some plots to make sure the estimations make sense
        p <- ggplot(sub_stats) +
            geom_point(aes(x=zscore, y=beta))
        savefile <- paste0(savedir, "zscore_beta_scatterplot.png")
        savefile
        ggsave(p, file=savefile)

        write_tsv(sub_stats, savefile)
    }


    head(sub_stats)
    # get columns that we need
    formatted_gwas <- sub_stats %>%
        separate(var, c('chrom', 'pos', 'ref', 'alt'), sep='_', remove=FALSE) %>%
        # remove chr from chrom
        mutate(chrom = str_remove(chrom, "chr")) %>%
        # make numeric and arrange by chrom and pos
        mutate(chrom=as.numeric(chrom),
               pos=as.numeric(pos)
        ) %>%
        arrange(chrom, pos) %>%
        mutate(ANNO='.') %>%
        # match header
        #select(chrom, POS=pos, REF=ref, ALT=alt, SNP_ID=var, N, ZSCORE=zscore, ANNO) 
        select(chrom, pos, ref, alt, snp_id=var, N, zscore, ANNO)
    head(formatted_gwas)

    # make sure it's in sorted order
    sorted_gwas <- formatted_gwas %>%
        arrange(chrom, pos) %>%
        mutate(chrom = paste0("chr", chrom))
    head(sorted_gwas)

    # check that we are overlapping variants with qtl
    check_qtl <- FALSE
    if (check_qtl){
        qtl_file <- paste0("/path/to/new_rosmap/output/07_ptwas/03_gambit_db/with_proportions/avgs_bulk_full_gene/all_gene.ptwas_weights.gambit.vcf")
        qtl_weights <- read_table(qtl_file, comment="##")
        head(qtl_weights)

        sub_qtls <- qtl_weights %>%
            select(chrom="#chr", pos, ref, alt) %>%
            mutate(qtl=TRUE) %>%
            unique()
        head(sub_qtls)
        dim(sub_qtls)

        joined <- sorted_gwas %>%
            inner_join(sub_qtls)
        head(joined)
        dim(joined)
    }


    # write to an output file
    (savefile <- paste0(savedir, gwas_type, "_gwas_file.vcf"))
    write_tsv(sorted_gwas, savefile, col_names=FALSE)

    # gzip the file
    cmd <- paste0("bgzip ", savefile)
    system(cmd)

    # try running tabix
    cmd <- paste0("tabix -p vcf ", savefile, ".gz")
    system(cmd)

    1
}


# read in the gwas
gwas_file <- paste0("/path/to/new_rosmap/output/06_fine_mapping_colocalization/02_gwas_prep/02_check_liftover/",
                    gwas_type, "/formatted_signed_dgap_gwas_scores.txt")
gwas <- read_table(gwas_file, col_names=c('var', 'zscore'))
head(gwas)



# get columns that we need
formatted_gwas <- gwas %>%
    separate(var, c('chrom', 'pos', 'ref', 'alt'), sep='_', remove=FALSE) %>%
    # remove chr from chrom
    mutate(chrom = str_remove(chrom, "chr")) %>%
    # make numeric and arrange by chrom and pos
    mutate(chrom=as.numeric(chrom),
           pos=as.numeric(pos)
    ) %>%
    arrange(chrom, pos) %>%
    mutate(N=94437, ANNO='.') %>%
    # match header
    #select(chrom, POS=pos, REF=ref, ALT=alt, SNP_ID=var, N, ZSCORE=zscore, ANNO) 
    select(chrom, pos, ref, alt, snp_id=var, N, zscore, ANNO)
head(formatted_gwas)

# make sure it's in sorted order
sorted_gwas <- formatted_gwas %>%
    arrange(chrom, pos) %>%
    mutate(chrom = paste0("chr", chrom))
head(sorted_gwas)

# check that we are overlapping variants with qtl
check_qtl <- FALSE
if (check_qtl){
    qtl_file <- paste0("/path/to/new_rosmap/output/07_ptwas/03_gambit_db/all_gene.ptwas_weights.gambit.vcf.gz")
    qtl_weights <- read_table(qtl_file, comment="##")
    head(qtl_weights)

    sub_qtls <- qtl_weights %>%
        select(chrom=`#chr`, pos, ref, alt) %>%
        mutate(qtl=TRUE) %>%
        unique()
    head(sub_qtls)
    dim(sub_qtls)

    joined <- sorted_gwas %>%
        inner_join(sub_qtls)
    head(joined)
    dim(joined)
}


# write to an output file

(savefile <- paste0(savedir, gwas_type, "_gwas_file.vcf"))
write_tsv(sorted_gwas, savefile, col_names=FALSE)


# gzip the file
cmd <- paste0("bgzip ", savefile)
system(cmd)

# try running tabix
cmd <- paste0("tabix -p vcf ", savefile, ".gz")
system(cmd)


