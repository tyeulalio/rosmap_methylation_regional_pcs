library(RNOmni)
library(tidyverse)

# format the methylation regional PCs in bed format
# for tensorqtl

datadir <- "/path/to/data"
savedir <- "/path/to/save/output"
include_study=FALSE

dir.create(savedir)

cleaned_str = "cleaned_global_"

load_methylation <- function(cell_type, region_type, summary_type, match_cpgs){
    # load methylation data
    if (summary_type == 'cpgs'){
        (datafile <- paste0(datadir, cleaned_str, cell_type, "_lifted_hg38_mvals.rds"))
        meth <- readRDS(datafile)
        head(meth)
    } else{
        basefile <- paste(region_type, cell_type, summary_type, sep='_')
        datafile <- paste0(datadir, cleaned_str, basefile, ".rds")
        datafile
        meth <- readRDS(datafile)

        head(meth)[1:10]
        tail(meth)[1:10]

        # remove columns that we don't need
        if (summary_type == 'avgs'){
            meth <- meth %>%
                select(-duration_secs)
        }

        if (summary_type == 'pcs'){
            # remove columns that are not pcs
            head(meth)[,1:10]
            if (region_type %in% c('full_gene', 'promoters', 'preTSS', 'gene_body')){
                meth <- meth %>%
                    select(!percent_variance_explained)
            } else{
                # some genes failed -- remove them
                meth = meth %>%
                    filter(!is.na(gv_dim)) %>%
                    select(-num_cpgs, -gv_dim, -mp_dim, -error, -duration_secs, -percent_variance)
                dim(meth)
            }
        }
    }
    head(meth)
    head(meth)[1:10]


    if (match_cpgs & summary_type == 'cpgs'){
        # match the cpgs to the regions
        # load cpg-gene map for region
        datafile <- paste0(datadir, region_type, '_', cell_type, '_gene_cpg_map.rds')
        cpg_map <- readRDS(datafile)
        head(cpg_map)

        head(meth)

        # only keep probes within gene region boundaries
        keep_probes <- intersect(rownames(meth), cpg_map$cpg_id)
        length(keep_probes)

        sub_meth <- meth[keep_probes,]
        dim(meth)
        dim(sub_meth)
        head(sub_meth)

        meth <- sub_meth
    }
    head(meth)

    meth
}

load_annotations <- function(region_type){
    
    datafile <- paste0("/path/to/annotate_regions/01_annotate_sites/", region_type, "_hg38_annotations.rds")
    if (region_type == 'full_gene'){
        region_annots <- readRDS(datafile)
        head(region_annots)

        unique(region_annots$type)

        # remove X and Y chromosomes
        # select all regions that fall within the promoter to the 5' UTR
        formatted_annots <- region_annots %>%
            filter(width != 0,
                   !seqnames %in% c('chrX', 'chrY')) %>%
            filter(type %in% c('hg38_genes_5UTRs', 'hg38_genes_promoters')) %>%
            mutate(type = gsub('hg38_genes_', '', type)) %>%
            gather('range', 'pos', start:end) %>%
            mutate(region_pos = paste(type, range, sep='_')) %>%
            select(gencode_gene_id, region_pos, seqnames, pos, strand) %>%
            group_by(region_pos, gencode_gene_id) %>%
            filter(row_number() == 1) %>%
            ungroup() %>%
            spread(region_pos, pos) %>%
            mutate(tss=ifelse(strand == '-', promoters_end, promoters_start))

            savefile <- paste0(savedir, region_type, "_tss_positions.rds")
            savefile
            saveRDS(formatted_annots, savefile)
    }
    if (region_type == 'promoters'){
        region_annots <- readRDS(datafile)
        head(region_annots)
        # add center of region position
        head(region_annots)
        formatted_annots <- region_annots %>%
            rowwise() %>%
            mutate(tss = max(start, end)) 
    }
    if (!region_type %in% c('full_gene', 'astro_enh', 'promoters')){
        # load TSS
        savefile <- paste0(savedir, 'full_gene', "_tss_positions.rds")
        formatted_annots <- readRDS(savefile)
        head(formatted_annots)
    }
    head(formatted_annots)


    # get the qtl window start and end
    tss_annots <- formatted_annots %>%
        select(gene_id=gencode_gene_id, chrom=seqnames, start=tss) %>%
        mutate(end=start+1) 
    head(tss_annots)

    tss_annots
}

convert_ids <- function(meth){
    # convert sample IDs to match genotype bed file
    samples <- names(meth)
    head(samples)

    # load matched clinical
    datafile <- paste0("../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
    matched_clinical <- readRDS(datafile)

    head(meth)

    dim(matched_clinical)
    dim(meth)


    # put the samples in the same order
    head(samples)

    head(matched_clinical)
    ordered_clinical <- matched_clinical[match(samples, matched_clinical$individualID),]
    head(ordered_clinical)

    identical(ordered_clinical$individualID, names(meth))

    # replace meth sample names with wgs.specimenID
    names(meth) <- ordered_clinical$wgs.specimenID
    head(meth)

    list(meth=meth, clinical=ordered_clinical)
}

create_bedfile <- function(converted_dat, annots, summary_type){
    # conver the data to a bedfile format
    # need columns: chrom, start, end, ... phenotypes for each sample
    names(converted_dat)

    meth <- converted_dat$meth %>%
        na.omit()

    head(meth)
    head(annots)


    # rank normalize using inverse normal transform (INT)
    #int_meth <- apply(meth, 1, RankNorm) %>% t()
    int_meth <- meth 
    head(int_meth)
    dim(int_meth)
    dim(meth)


    # create bed with normalized methylation
    if (summary_type == 'avgs'){
        bed_formatted <- int_meth %>%
            as.data.frame() %>%
            rownames_to_column('gene_id') %>%
            left_join(annots) %>%
            select(chrom, start, end, gene_id, everything()) %>%
            rename('#chr'=chrom)
    }
    if (summary_type == 'pcs'){
        bed_formatted <- int_meth %>%
            as.data.frame() %>%
            rownames_to_column('gene_pc') %>%
            separate(gene_pc, c('gene_id', 'pc'), sep='-') %>%
            left_join(annots) %>%
            unite(gene_pc, gene_id:pc, sep='-') %>%
            select(chrom, start, end, gene_id=gene_pc, everything()) %>%
            rename('#chr'=chrom)
    }
    if (summary_type == 'cpgs'){
        bed_formatted <- int_meth %>%
            as.data.frame() %>%
            rownames_to_column('cpg_id') %>%
            separate(cpg_id, c('chrom', 'start', 'end'), sep='_', remove=FALSE) %>%
            select(chrom, start, end, gene_id=cpg_id, everything()) %>%
            mutate(end = as.numeric(start)+1) %>%
            rename('#chr'=chrom)
    }
    head(bed_formatted)[,1:10]


    bed_formatted
}

save_bed <- function(formatted_bed, cell_type, region_type, summary_type){
    # save the bedfile

    savefile <- paste(cell_type, region_type, "formatted", summary_type, sep='_')
    savepath <- paste0(savedir, cleaned_str, savefile, ".bed")
    savepath

    write_tsv(formatted_bed, savepath)

    # gzip the output file
    if (file.exists(paste0(savepath, ".gz"))) system(paste0("rm ", savepath, ".gz"))
    system(paste("gzip", savepath))

    1
}

main <- function(){
    # load the methylation PCs
    region_types <- c('preTSS', 'gene_body', 'full_gene', 'promoters')
    cell_types <- c('astro', 'bulk', 'endo', 'neuron', 'oligo_opc')
    summary_types <- c('avgs', 'pcs')

    region_type <- region_types[1]
    cell_type <- cell_types[1]
    summary_type <- summary_types[1]

    # get all combinations of parameters to run
    runs <- expand.grid(region_type = region_types,
                        cell_type = cell_types,
                        summary_type = summary_types
    ) %>%
        mutate(row_number= row_number())
    head(runs)
    runs

    run = runs[1,]
    process_run <- function(run){
        cell_type <- run[['cell_type']]
        region_type <- run[['region_type']]
        summary_type <- run[['summary_type']]
        rn <- run[['row_number']]

        print(paste("Processing", cell_type, region_type, summary_type, rn, "out of", nrow(runs), Sys.time()))

        # load methylation data
        print(paste("Loading methylation", Sys.time()))
        meth <- load_methylation(cell_type, region_type, summary_type, match_cpgs=FALSE)
        head(meth)

        # load region annotations
        print(paste("Loading annotations", Sys.time()))
        annots <- load_annotations(region_type)

        # convert IDs
        print(paste("Converting IDs", Sys.time()))
        converted_dat <- convert_ids(meth)

        # format bed files
        print(paste("Formatting bed", Sys.time()))
        formatted_bed <- create_bedfile(converted_dat, annots, summary_type) 
        head(formatted_bed)[1:10]
        
        # save the bed file
        print(paste("Saving bed", Sys.time()))
        save_bed(formatted_bed, cell_type, region_type, summary_type)

        1
    }

    apply(runs, 1, process_run)


    # do one run for genome-wide CpGs
    cell_type <- 'astro'
    process_run <- function(cell_type){
        summary_type <- "cpgs"
        region_type <- "genome"

        print(paste("Processing genome Cpgs", cell_type, region_type, summary_type, Sys.time()))

        # load methylation data
        print(paste("Loading methylation", Sys.time()))
        meth <- load_methylation(cell_type, region_type, summary_type, match_cpgs=FALSE)
        head(meth)
        dim(meth)

        # convert IDs
        print(paste("Converting IDs", Sys.time()))
        converted_dat <- convert_ids(meth)

        # format bed files
        print(paste("Formatting bed", Sys.time()))
        formatted_bed <- create_bedfile(converted_dat, annots=NA, summary_type) 
        head(formatted_bed)
        dim(formatted_bed)
        
        # save the bed file
        print(paste("Saving bed", Sys.time()))
        save_bed(formatted_bed, cell_type, region_type, summary_type)

        1
    }

    #process_run(cell_type)
    lapply(cell_types, process_run)
    1
}

main()
