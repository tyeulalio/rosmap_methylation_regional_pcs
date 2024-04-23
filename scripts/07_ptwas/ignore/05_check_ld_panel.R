library(liftOver)
library(GenomicRanges)
library(tidyverse)

# check if we have any matches between LD panel and GWAS 
# LD panel probably has to be lifted

gwas_type="wightman"
savedir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/05_lifted_ld_panel/", gwas_type, "/")
dir.create(savedir)

ld_savedir <- paste0(savedir, "G1K_EUR_3V5_hg38/")
dir.create(ld_savedir)

gwas_savedir <- paste0(savedir, "updated_gwas/")
dir.create(gwas_savedir)

# read in the gwas file
datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/04_formatted_gwas/", gwas_type, "_gwas_file.vcf.gz")
full_gwas <- read_table(datafile, col_names=c('chrom', 'pos', 'ref', 'alt', 'var', 'n', 'zscore', 'anno'))
head(full_gwas)

chromosome <- 2 # for testing
args <- commandArgs(trailingOnly=TRUE)
chromosome <- args[1] 

process_chromosome <- function(chromosome){ 
    print(paste("Processing chromosome", chromosome, Sys.time()))

    # filter gwas for this chromosome
    gwas <- full_gwas %>%
        filter(chrom == paste0('chr', chromosome))
    head(gwas)

    # read in chromosome 22 LD panel
    print(paste("Reading LD panel", Sys.time()))
    ld_file <- paste0("/home/eulalio/labs/jachang4/ptwas_quickstart/G1K_EUR_3V5/chr", chromosome, ".vcf.gz")
    ld_panel <- read_table(ld_file, comment="##")
    head(ld_panel)

    # try lifting ld panel to hg38
    print(paste("Lifting LD panel", Sys.time()))
    head(ld_panel)
    colnames(ld_panel)[1] <- "chrom"
    formatted_ld <- ld_panel %>%
        mutate(start=as.numeric(POS),
               end = start+1
        ) %>%
        mutate(chrom = paste0("chr", chrom)) 
    head(formatted_ld)

    ld_gr <- makeGRangesFromDataFrame(formatted_ld, keep.extra=TRUE)

    chain <- import.chain("/home/eulalio/deconvolution/data/hg19ToHg38.over.chain")
    ld_hg38 <- liftOver(ld_gr, chain=chain) %>%
        as.data.frame() %>%
        select(-group, -group_name, -width, -strand, -POS, -end) %>%
        arrange(start) %>%
        #update the ID with the new position
        mutate(ID = paste(seqnames, start, REF, ALT, sep='_'))
    head(ld_hg38)[1:10]
    

    # check for duplicate variants
    dups <- ld_hg38 %>%
        unique() %>%
        group_by(ID) %>%
        summarise(count=n()) %>%
        arrange(-count) %>%
        filter(count > 1)
    head(dups)
    dim(dups)

    id <- dups[1,]$ID
    ld_hg38[ld_hg38$ID == id,1:10]

    # these variants mapped to the same position after lifting over
    # just delete them, there's only 4 for chrom 22
    print(paste("FILTERING OUT SNPS THAT MAPPED TO ANOTHER CHROMOSOME"))
    ld_nodups <- ld_hg38 %>%
        filter(!ID %in% dups$ID) %>%
        unique() 
    head(ld_nodups)[1:10]
    dim(ld_nodups)

    unique(ld_nodups$ID) %>% length()

    head(ld_nodups)[1:10]
    ld_chrom <- ld_nodups %>%
        filter(seqnames  == paste0("chr", chromosome)) %>%
        arrange(start)
    table(ld_chrom$seqnames)

    colnames(ld_chrom)[1] <- "#CHROM"
    colnames(ld_chrom)[2] <- "POS"
    head(ld_chrom)[1:10]



    # write this lifted version to output
    (savefile <- paste0(ld_savedir, "chr", chromosome, ".vcf"))
    # write header to file first
    (cmd <- paste0("zcat ", ld_file, " | grep '##' > ", savefile))
    system(cmd)

    write_tsv(ld_chrom, savefile, append=TRUE, col_names=TRUE)

    # gzip the file
    cmd <- paste0("bgzip -f ", savefile)
    system(cmd)

    # try running tabix
    cmd <- paste0("tabix -p vcf ", savefile, ".gz")
    system(cmd)
    savefile




    print(paste("Matching GWAS and LD", Sys.time()))
    # now try connecting them 
    head(gwas)
    head(ld_chrom)[1:10]

    # get only the columns we need now
    colnames(ld_chrom)[1] <- "chrom"
    sub_ld <- ld_chrom %>%
        select(chrom, pos=POS, ref=REF, alt=ALT, ID)
    head(sub_ld)
    
    joined_hg38 <- gwas %>%
        inner_join(sub_ld, by=c('chrom', 'pos'), suffix=c('.gwas', '.ld'))
    head(joined_hg38)
    dim(joined_hg38) # now we match a lot of them
    dim(gwas)
    dim(ld_nodups)

    # check how many have swapped ref/alt
    dim(joined_hg38)
    head(joined_hg38)

    matched <- joined_hg38 %>%
        filter(ref.gwas == ref.ld, 
               alt.gwas == alt.ld
        )
    head(matched)
    dim(matched)
    dim(joined_hg38) # a lot of them match exactly

    swapped <- joined_hg38 %>%
        filter(ref.gwas == alt.ld, alt.gwas == ref.ld)
    head(swapped)
    dim(swapped)

    # write these to a file 
    savefile <- paste0(savedir, "chr", chromosome, "_swapped_variants.rds")
    saveRDS(swapped, savefile)
    #swapped <- readRDS(savefile)

    total <- nrow(matched) + nrow(swapped)
    nrow(joined_hg38) - total # number of mismatched ref/alt - not too bad, about 100 - just drop these

    # create a file for swapping using plink
    # columns should be:
    # 1. variant id
    # 2. old allele1
    # 3. old allele2
    # 4. new allele1
    # 5. new allele2
    head(swapped)
    plink_swapped <- swapped %>%
        select(ID, old1=alt.gwas, old2=ref.gwas, new1=ref.gwas, new2=alt.gwas) 
    head(plink_swapped)

    savefile <- paste0(savedir, "plink_formatted_swapped_chr", chromosome, ".txt")
    write_delim(plink_swapped, savefile, col_names=FALSE, delim=' ')

    1
}

res <- process_chromosome(chromosome)
#res <- lapply(1:22, process_chromosome)
