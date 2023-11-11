library(RnBeads.hg19)
library(RnBeads.hg38)
library(RnBeads)
library(GenomicRanges)
library(liftOver)
library(biomaRt)
library(PCAtools)
library(RNOmni)
library(tidyverse)



# summarise gene regions using cpgs, averages, and the pc method


sub_chroms = FALSE
c1_covariates = FALSE
c2_covariates = TRUE
reestimate_props = TRUE


jobnum=10

#args = commandArgs(trailingOnly=TRUE)

#jobnum <- args[1] %>%
    #as.numeric()


main_datadir <- paste0("/home/eulalio/deconvolution/new_rosmap/output/02_deconvolution/02_svd_batch_correction/")
main_savedir <- paste0("../../output/03_summarised_genes/01_summarised_genes/")
dir.create(main_savedir, showWarnings=FALSE)
dir.create(main_savedir, showWarnings=FALSE)

paths_file <- paste0("/home/eulalio/deconvolution/new_rosmap/scripts/shared_functions/01_filepaths.R")
include_study=FALSE # set this because we always include it anyways
source(paths_file)

datadir
savedir

# remove randomize_cpgs from datadir
datadir <- str_remove(datadir, "randomize_cpgs/")
datadir

cleaned_mvals = TRUE
if (cleaned_mvals){
    cleaned_str = "cleaned_"
} else{
    cleaned_str = ""
}

clean_global <- TRUE
if (clean_global){#  # Tue May 23 18:16:32 2023 ------------------------------
    cleaned_str <- "cleaned_global_"
}

precleaned = FALSE
if (precleaned){
    precleaned_str = "precleaned_"
    savedir <- paste0(savedir, "precleaned/")
    dir.create(savedir)
} else{
    precleaned_str = ""
}

datadir
savedir

summarize_results = FALSE



load_methylation <- function(cell_type){
    # read in the mvals before and after batch correction

    if (cell_type == 'bulk' | cleaned_mvals){
        # use the cleaned bulk data - batch removed with removeBatchEffect
        list.files(datadir)
        datafile <- paste0(datadir, precleaned_str, cleaned_str, cell_type, "_cleaned_meth.rds")
    } else{
        # load appropriate methylation
        datafile <- paste0(datadir, precleaned_str, cell_type, "_formatted_meth.rds")
    }
    if (clean_global){
        datafile <- paste0(datadir, precleaned_str, cleaned_str, cell_type, "_cleaned_meth.rds")
    }

    datafile
    mvals <- readRDS(datafile)
    head(mvals)
    dim(mvals)

    mvals
}

lift_methylation <- function(mvals, cell_type){
    # liftOver the methylation data from 
    # hg19 to hg38


    # import the chain file
    datafile <- paste0("/home/eulalio/deconvolution/data/hg19ToHg38.over.chain")
    ch = import.chain(datafile)


    # use the rnbeads hg19 annotations    
    cpg_annots <- rnb.annotation2data.frame(rnb.get.annotation(type='probes450', assembly='hg19'))
    head(cpg_annots)

    head(mvals)

    cpgs <- rownames(mvals)
    head(cpgs)

    #if (cell_type == 'bulk'){
        tmp <- cpg_annots %>%
            mutate(cpg = paste(Chromosome, Start, sep='_')) %>%
            rownames_to_column('rn') %>%
            column_to_rownames('cpg')
        head(tmp)

        sub_annots <- tmp[cpgs,] 
        head(sub_annots)
    #} else{
        #sub_annots <- cpg_annots[cpgs,] 
        #head(sub_annots)
    #}


    dim(sub_annots)
    dim(mvals)


    # lift over cpg annotations to hg38
    head(sub_annots)
    mvals_gr <- makeGRangesFromDataFrame(sub_annots)
    mvals_gr

    mvals38 <- liftOver(mvals_gr, ch) %>% as.data.frame()
    head(mvals38)

    # we lose some cpgs after the liftOver
    overlapped_cpgs <- intersect(mvals38$group_name, rownames(mvals))
    length(overlapped_cpgs)

    # create a column for cpg position - this will be the new cpg name
    cpgs38 <- mvals38 %>%
        mutate(cpg=paste(seqnames, start, end, sep='_')) %>%
        rename(cpg_gr37 = group_name)
    head(cpgs38)

    savefile <- paste0(savedir, cell_type, "_cpg_position_map.rds")
    savefile
    saveRDS(cpgs38, savefile)

    # order both data the same
    ordered_cpgs <- cpgs38[match(overlapped_cpgs, cpgs38$cpg_gr37),]
    ordered_mvals <- mvals[overlapped_cpgs,]
    
    # these need to match - both should be TRUE
    identical(ordered_cpgs$cpg_gr37, overlapped_cpgs)
    identical(rownames(ordered_mvals), overlapped_cpgs)

    # attach the positional info as the rownames for cpgs
    rownames(ordered_mvals) <- ordered_cpgs$cpg
    head(ordered_mvals)
    
    ordered_mvals
}

get_annotations <- function(region_type){
    # use the rnbeads hg38 annotations    
    if (region_type == 'enhancers'){
        datadir <- paste0("/home/eulalio/deconvolution/rosmap/output/annotate_regions/01_annotate_sites/")
        datafiles <- list.files(datadir) %>% str_subset('enhancer_annots')
        datafiles

        # load all of the enhancer annots
        datapaths <- paste0(datadir, datafiles)
        annots <- lapply(datapaths, readRDS)
        annots_df <- do.call(rbind, annots)
        head(annots_df)

        length(annots_df$enh)
        length(unique(annots_df$enh))


        # separate enhancers by type by adding type to gene name
        formatted_df <- annots_df %>%
            mutate(gene_id = paste(gene_id, type, sep='-')) 
        head(formatted_df)
        dim(formatted_df)

        return(formatted_df)

    } else{
        if (region_type == 'gene_body'){
            # load the full gene results
            datafile <- paste0("/home/eulalio/deconvolution/rosmap/output/annotate_regions/01_annotate_sites/", "full_gene", "_hg38_annotations.rds")
            annots <- readRDS(datafile)
            head(annots)

            # filter to regions after the TSS
            unique(annots$type)

            filtered_annots <- annots %>%
                filter(!type %in% c('hg38_genes_1to5kb', 'hg38_genes_promoters'))
            unique(filtered_annots$type)

            annots <- filtered_annots
        } else{
            if (region_type == 'preTSS'){
                # load the full gene results
                datafile <- paste0("/home/eulalio/deconvolution/rosmap/output/annotate_regions/01_annotate_sites/", "full_gene", "_hg38_annotations.rds")
                annots <- readRDS(datafile)
                head(annots)

                # filter to regions after the TSS
                unique(annots$type)

                filtered_annots <- annots %>%
                    filter(type %in% c('hg38_genes_1to5kb', 'hg38_genes_promoters'))
                unique(filtered_annots$type)

                annots <- filtered_annots
            } else{
                # load the compiled annotations
                datafile <- paste0("/home/eulalio/deconvolution/rosmap/output/annotate_regions/01_annotate_sites/", region_type, "_hg38_annotations.rds")
                annots <- readRDS(datafile)
                head(annots)
            }

        }
    }




    #if (region_type == 'expanded_promoter'){
        #annots <- rnb.annotation2data.frame(rnb.get.annotation(type='promoters', assembly='hg38'))
        #head(annots)
        #enhancer_range = 1e6

         ## get granged for annotations
        #annots_gr <- makeGRangesFromDataFrame(annots, keep.extra.columns=TRUE)
        #head(annots_gr)

        ## expand the range

        ## use this to extend the range of the genomic regions
        #extend <- function(x, upstream=0, downstream=0)     
        #{
            #if (any(strand(x) == "*"))
                #warning("'*' ranges were treated as '+'")
            #on_plus <- strand(x) == "+" | strand(x) == "*"
            #new_start <- start(x) - ifelse(on_plus, upstream, downstream)
            #new_end <- end(x) + ifelse(on_plus, downstream, upstream)
            #ranges(x) <- IRanges(new_start, new_end)
            #trim(x)
        #}


        ## expand the range for each promoter region
        #extended_gr <- extend(annots_gr, enhancer_range, enhancer_range)
        #extended_gr

         #annots_gr

        ## get dataframe from gr
        #annots <- data.frame(extended_gr)

    #} else{
        #annots <- rnb.annotation2data.frame(rnb.get.annotation(type=region_type, assembly='hg38'))
    #}

    head(annots)
    unique(annots$type)

    # attach gene biotype info 
    #savefile <- paste0(savedir, "gene_biomart_annots.csv")
    savefile <- paste0("/home/eulalio/deconvolution/data/gene_biomart_table.csv")
    if (!file.exists(savefile)){
        ensembl <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
        latt <- listAttributes(ensembl)
        latt[1:100,]
        atts <- c('ensembl_gene_id_version', 'ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype')

        full_bm <- getBM(
                    attributes = atts,
                    mart = ensembl
                    )
        head(full_bm)

        # rename bm columns to match
        gene_annots <- full_bm %>%
            rename(gene_id=ensembl_gene_id, symbol=hgnc_symbol, biotype=gene_biotype) 
        head(gene_annots)

        #savefile <- paste0(savedir, "gene_biomart_annots.csv")
        #savefile
        write_csv(gene_annots, savefile)
    }

    gene_annots <- read_csv(savefile)

    head(gene_annots)
    head(annots)

    if (region_type == 'astro_enh'){
        full_annots <- annots %>%
            #rename(ensembl_gene_id_version = gencode_gene_id) %>%
            left_join(gene_annots) %>%
            select(-gene_id, gene_id = ensembl_gene_id_version) %>%
            #separate(gencode_region, c('chrom', 'start', 'end'), sep='_') %>%
            #select(chrom, start, end, everything()) %>%
            select(-strand)
    } else{
        full_annots <- annots %>%
            rename(ensembl_gene_id_version = gencode_gene_id) %>%
            left_join(gene_annots) %>%
            select(-gene_id, gene_id = ensembl_gene_id_version) %>%
            #separate(gencode_region, c('chrom', 'start', 'end'), sep='_') %>%
            #select(chrom, start, end, everything()) %>%
            select(-strand)
    }
    head(full_annots)

    table(full_annots$biotype, useNA='ifany')
    region_type

    # filter down to protein coding for now
    #keep_annots <- full_annots %>%
        #filter(biotype == 'protein_coding')
    #head(keep_annots)
    #dim(keep_annots)

    annots_savedir <- paste0(savedir, "annotated_regions/")
    dir.create(annots_savedir, showWarnings=FALSE)
    savefile <- paste0(annots_savedir, region_type, "_annotations.rds")
    savefile
    saveRDS(full_annots, savefile)

    full_annots
}



#get_annotations <- function(region_type){
    ## use the rnbeads hg38 annotations    
    #if (region_type == 'enhancers'){
        #savefile <- paste0(savedir, "abc_enhancers_hg38.csv")
        #savefile
        #if (file.exists(savefile)){
            #print(paste("Enhancer annotations exist. Loading file."))
            #keep_annots <- read_csv(savefile)
            #return(keep_annots)
        #}
        ## read in abc enhancer file
        #datafile <- paste0("/home/eulalio/deconvolution/data/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt")
        #abc_enhancers <- read_tsv(datafile)
        #head(abc_enhancers)

        ## celltypes of interest
        #abc_celltypes <- c('astrocyte-ENCODE', 'endothelial_cell_of_umbilical_vein-Roadmap', 
                           #'bipolar_neuron_from_iPSC-ENCODE', 'H1_Derived_Neuronal_Progenitor_Cultured_Cells-Roadmap')
        #filtered_abc <- abc_enhancers %>%
            #filter(CellType %in% abc_celltypes)
        #head(filtered_abc)
        
        #summary(filtered_abc$ABC.Score)

        ## use cut-offs from paper
        ## ABC.Score >= 0.02
        #scored_abc <- filtered_abc %>%
            #filter(ABC.Score >= 0.02)
        #dim(filtered_abc)
        #dim(scored_abc)

        ## need to liftover to hg38 from hg19
        #chain <- import.chain("/home/eulalio/deconvolution/data/hg19ToHg38.over.chain")
        #gr_abc <- makeGRangesFromDataFrame(scored_abc, keep.extra.columns=TRUE)
        #head(gr_abc)

        #gr_abc_hg38 <- liftOver(gr_abc, chain)
        #head(gr_abc_hg38)

        ## convert lifted data to data frame
        #abc_hg38 <- as.data.frame(gr_abc_hg38)
        #head(abc_hg38)

        ## convert target gnee symbols to enseml ids
        #bm_file <- paste0("/home/eulalio/deconvolution/data/full_biomart_table.csv")
        #bm <- read_csv(bm_file)
        #head(bm)
        #head(abc_hg38)

        #abc_hg38_genes <- abc_hg38 %>%
            #mutate(Symbol = TargetGene) %>%
            #left_join(bm)
        #head(abc_hg38_genes)

        #unique(abc_hg38_genes$Gene_biotype)

        ## filter for protein coding genes
        #abc_hg38_formatted <- abc_hg38_genes %>%
            #filter(Gene_biotype == 'protein_coding') %>%
            #select(-group, -group_name)
        #head(abc_hg38_formatted)
        #dim(abc_hg38_formatted)
        #dim(abc_hg38_genes)



        ## get the important columns to match the other annotations
        #head(abc_hg38_formatted)
        ##head(promoter_annots)

        #keep_annots <- abc_hg38_formatted %>%
            #select(Chromosome=seqnames, Start=start, End=end, Strand=strand,
                    #symbol=Symbol, gene_id=Ensembl_ID, biotype=Gene_biotype, celltype=CellType) %>%
            #mutate(gene_id = paste0(gene_id, "__", celltype))
        #head(keep_annots)

        ## save the lifted file
        #savefile <- paste0(savedir, "abc_enhancers_hg38.csv")
        #savefile
        #write_csv(keep_annots, savefile)

        #return(keep_annots)
    #}

    #if (region_type == 'expanded_promoter'){
        #annots <- rnb.annotation2data.frame(rnb.get.annotation(type='promoters', assembly='hg38'))
        #head(annots)
        #enhancer_range = 1e6

         ## get granged for annotations
        #annots_gr <- makeGRangesFromDataFrame(annots, keep.extra.columns=TRUE)
        #head(annots_gr)

        ## expand the range

        ## use this to extend the range of the genomic regions
        #extend <- function(x, upstream=0, downstream=0)     
        #{
            #if (any(strand(x) == "*"))
                #warning("'*' ranges were treated as '+'")
            #on_plus <- strand(x) == "+" | strand(x) == "*"
            #new_start <- start(x) - ifelse(on_plus, upstream, downstream)
            #new_end <- end(x) + ifelse(on_plus, downstream, upstream)
            #ranges(x) <- IRanges(new_start, new_end)
            #trim(x)
        #}


        ## expand the range for each promoter region
        #extended_gr <- extend(annots_gr, enhancer_range, enhancer_range)
        #extended_gr

         #annots_gr

        ## get dataframe from gr
        #annots <- data.frame(extended_gr)

    #} else{
        #annots <- rnb.annotation2data.frame(rnb.get.annotation(type=region_type, assembly='hg38'))
    #}

    #head(annots)

    ## attach gene biotype info 
    #savefile <- paste0(savedir, "gene_biomart_annots.csv")
    #if (!file.exists(savefile)){
        #ensembl <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
        #latt <- listAttributes(ensembl)
        #latt[1:100,]
        #atts <- c('ensembl_gene_id_version', 'ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype')

        #full_bm <- getBM(
                    #attributes = atts,
                    #mart = ensembl
                    #)
        #head(full_bm)

        ## rename bm columns to match
        #gene_annots <- full_bm %>%
            #rename(gene_id=ensembl_gene_id, symbol=hgnc_symbol, biotype=gene_biotype) 
        #head(gene_annots)

        #savefile <- paste0(savedir, "gene_biomart_annots.csv")
        #savefile
        #write_csv(gene_annots, savefile)
    #}

    #gene_annots <- read_csv(savefile)

    #head(gene_annots)
    #head(annots)

    #full_annots <- annots %>%
        #rename(gene_id = ID, rnbeads_symbol=symbol) %>%
        #left_join(gene_annots)
    #head(full_annots)

    #table(full_annots$biotype, useNA='ifany')
    #region_type

    ## filter down to protein coding for now
    #keep_annots <- full_annots %>%
        #filter(biotype == 'protein_coding')
    #head(keep_annots)
    #dim(keep_annots)


    #keep_annots
#}


get_pcs <- function(gene_mvals){
    # estimate dims and get sig pcs for this gene
    start_time = Sys.time()
    
    head(gene_mvals)
    dim(gene_mvals)

    # get PCs
    # rows = variables
    # cols = samples
    pca_res <- PCAtools::pca(gene_mvals)
    names(pca_res)

    # get dim using gavish donoho
    # rows = variables
    # cols = observations
    gv_res <- chooseGavishDonoho(gene_mvals, var.explained=pca_res$sdev^2, noise=1)
    gv_res

    mp_res <- chooseMarchenkoPastur(gene_mvals, var.explained=pca_res$sdev^2, noise=1)
    mp_res

    #print(paste("gavish donoho dim:", gv_res[[1]], "marchenko pastur dim:", mp_res[[1]]))

    names(pca_res)
    pcs <- pca_res$rotated
    head(pcs)

    # get the loadings also
    loads <- pca_res$loadings
    head(loads)

    est_dim <- gv_res[[1]]
    est_dim
    

    # filter down to estimated dimensions
    sig_pcs <- pcs[,1:est_dim,drop=FALSE] %>%
        as.data.frame()
    head(sig_pcs)
    dim(sig_pcs)
    sig_loads <- loads[,1:est_dim,drop=FALSE] %>%
        as.data.frame()
    head(sig_loads)


    # get pct var explained
    eig_sq <- pca_res$sdev ^ 2
    pct_var <- eig_sq / sum(eig_sq)
    subset_pct_var <- data.frame(percent_variance_explained=pct_var[1:est_dim, drop=FALSE]) %>%
        t() %>%
        as.data.frame()
    names(subset_pct_var) <- names(sig_pcs)

    head(subset_pct_var)
    head(sig_pcs)

    formatted_pcs <- rbind(subset_pct_var, sig_pcs)
    head(formatted_pcs)

    mp_var <- sum(pct_var[1:mp_res[[1]]])
    #print(paste("Variance explained, GV:", round(sum(subset_pct_var), 2), "MP:", round(sum(mp_var), 2)))

    end_time = Sys.time()
    duration = as.numeric(end_time - start_time, "secs")

    # return sig pcs and estimated dim
    pc_res <- list(sig_pcs=formatted_pcs, loadings=sig_loads, gv_dim=gv_res[[1]], mp_dim=mp_res[[1]], error=NA, duration_secs = duration)
    pc_res
}

#mvals <- cleaned_mvals
compute_averages <- function(mvals, gene_map, region_type, cell_type, cleaned_str){
    # summarise the regions using averages

    #savefile <- paste0(savedir, cleaned_str, region_type, cell_type, "_avg_mvals.rds")
    #if (file.exists(savefile)){
       #return(NA) 
    #}

    savefile <- paste0(savedir, cleaned_str, region_type, cell_type, "_avg_mvals.rds")
    savefile
    if (file.exists(savefile)){
        return(1)        
    }

    head(mvals)
    head(gene_map)


    # grab genes and samples
    genes <- unique(gene_map$gene_id)
    samples <- colnames(mvals)

    # mvals are NOT centered and scaled prior
    # to computing averages
    gene <- genes[6]
    get_gene_avg <- function(gene){
        gene_idx = which(gene == genes)
        if (gene_idx %% 100 == 0) 
            print(paste("Computing average for gene", gene_idx, "out of", length(genes), Sys.time()))
        # get cpgs in this gene
        gene_cpgs <- gene_map %>%
            filter(cpg_id %in% rownames(mvals)) %>%
            filter(gene_id == gene)
        gene_cpgs

        if (nrow(gene_cpgs) < 1) return(NA)

        # get the mvals for this gene
        gene_mvals <- mvals[gene_cpgs$cpg_id,,drop=FALSE]
        head(gene_mvals)

        if (nrow(gene_mvals) == 1){
            return_mvals <- (gene_mvals[,samples,drop=FALSE])
            rownames(return_mvals) <- gene
        }

        # compute averages
        gene_avgs <- data.frame(apply(gene_mvals, 2, mean)) %>% t()
        gene_avgs

        # make sure the samples are ordered
        ordered_avgs <- gene_avgs[,samples,drop=FALSE]
        rownames(ordered_avgs) <- gene
        ordered_avgs
    }

    avgs <- mclapply(genes, get_gene_avg, mc.cores=16)
    avgs_df <- do.call(rbind, avgs)

    # compute region averages
    #region_mvals <- mvals %>%
        #as.data.frame() %>%
        #rownames_to_column('cpg_id') %>%
        #inner_join(gene_map) %>%
        #gather('specimenID', 'mval', -gene_id, -cpg_id)
    #head(region_mvals)

    #avgs <- region_mvals %>%
        #group_by(gene_id, specimenID) %>%
        #summarise(avg=mean(mval))
    #head(avgs)

    avgs_df <- na.omit(avgs_df)
    head(avgs_df)
    dim(avgs_df)

    # save these averages
    savefile <- paste0(savedir, cleaned_str, region_type, cell_type, "_avg_mvals.rds")
    savefile
    saveRDS(avgs_df, savefile)

    1
}

select_genes <- function(array_idx, mvals, gene_map, region_type, cell_type, cleaned_str){
    # summarise the regions using averages and PCs

    head(mvals)
    head(gene_map)


    #pc_dir <- paste0(savedir, cleaned_str, region_type, "_", cell_type, "_summaries/")
    #savefile <- paste0(pc_dir, cleaned_str, region_type, "_", cell_type, "_gene_avgs_", array_idx, ".rds")
    #savefile
    #if (file.exists(savefile)){
        #print("output exists. returning.")
        #return(NA)
    #}


    # mvals are NOT centered and scaled prior
    # to computing averages

    # compute region averages
    head(mvals)
    region_mvals <- mvals %>%
        as.data.frame() %>%
        rownames_to_column('cpg_id') %>%
        select(cpg_id) 
    head(region_mvals)

    head(gene_map)

    region_map <- gene_map %>% 
        filter(cpg_id %in% region_mvals$cpg_id)
    dim(region_map)
    dim(gene_map)
    head(region_map)


    # figure out which genes to process
    genes <- unique(region_map$gene_id)
    num_genes <- length(genes)
    num_genes

    genes_per_job <- ceiling(num_genes/total_jobs)
    genes_per_job

    start_idx <- array_idx * genes_per_job + 1
    end_idx <- min(start_idx + genes_per_job - 1, num_genes) # don't go past the end of the gene array

    if (start_idx > num_genes){
        print("All jobs complete.")
        return(NA)
    }

    job_genes <- genes[start_idx:end_idx]

    # filter the map down to the genes we want to summarize
    job_map <- region_map %>%
        filter(gene_id %in% job_genes)
    head(job_map)
    length(unique(job_map$gene_id))

    print(paste("total jobs", total_jobs, "genes per job", genes_per_job))
    print(paste("Job", array_idx, "processing genes", start_idx, "to", end_idx, "out of", num_genes))

    job_map
}

#summary_type = "pcs"
summarise_genes <- function(mvals, job_map, summary_type){ 
    print(paste("Summarizing genes with", summary_type))

    genes <- unique(job_map$gene_id)

    # for testing
    gene <- genes[3]
    gene_cpgs <- job_map %>%
        filter(gene_id == gene)
    gene_cpgs

    # compute region PCs
    ordered_samples <- colnames(mvals)
    get_gene_summaries <- function(gene){

        print(paste("Computing", summary_type, "for gene", which(gene == genes), "out of", length(genes), Sys.time()))

        gene_cpgs <- job_map %>%
            filter(gene_id == gene) %>%
            unique()
        gene_cpgs
        dim(gene_cpgs)

        # grab the mvals for this gene
        gene_mvals <- as.data.frame(mvals)[gene_cpgs$cpg_id,,drop=FALSE]
        head(gene_mvals)
        dim(gene_mvals)

        num_cpgs <- nrow(gene_mvals)

        if (summary_type == 'avgs'){
            # get this gene's avgs
            start_time = Sys.time()
            avgs <- data.frame(avg=colMeans(gene_mvals)) %>% t() %>% as.data.frame()
            #rownames(avgs) <- gene
            avgs$gene_id = gene
            head(avgs)

            end_time = Sys.time()
            duration = as.numeric(end_time - start_time, "secs")
            avgs$duration_secs = duration
            rownames(avgs) = NULL
            head(avgs)
            return(avgs)
        }

        if (summary_type == 'pcs'){
            # get the pcs for this gene
            pc_res <- get_pcs(gene_mvals)
            names(pc_res)

            # no error occurred
            if (all(is.na(pc_res$error))){
                pc_res$gene_id = gene

                # format the pcs to combine later
                gene_pcs <- pc_res$sig_pcs
                head(gene_pcs)

                formatted_pcs <- gene_pcs[c('percent_variance_explained', ordered_samples),,drop=FALSE] %>% 
                    t() %>%
                    as.data.frame()
                rownames(formatted_pcs) <- paste(gene, rownames(formatted_pcs), sep='-')
                head(formatted_pcs)

                # do the same for loadings
                names(pc_res)
                gene_loads <- pc_res$loadings
                colnames(gene_loads) <- paste(gene, colnames(gene_loads), sep='-')
                head(gene_loads)

                # gather the other info into a df
                pc_info = data.frame(gene_id=gene,
                                     gv_dim=pc_res$gv_dim,
                                     mp_dim=pc_res$mp_dim,
                                     duration_secs=pc_res$duration_secs
                )
                head(pc_info)


                return(list(formatted_pcs=formatted_pcs,
                            loadings=gene_loads,
                            pc_info=pc_info,
                            gene_id=gene))
            } 
            #else{ # error occurred
                #formatted_pcs = data.frame(row.names=c(ordered_samples, 'num_cpgs', 'est_dim', 'error', 'duration_secs'), 
                                           #gene=rep(NA, length(ordered_samples)+4)) %>% t() %>% as.data.frame()
                #formatted_pcs$error = as.character(pc_res$error)
                #formatted_pcs$num_cpgs = num_cpgs
                #rownames(formatted_pcs) = gene
                #head(formatted_pcs)
            #}
        }
    }

    # run for all of the job's genes
    all_summaries <- lapply(genes, get_gene_summaries)

    print(paste("Combining gene summaries"))
    if (summary_type == 'avgs'){
        # combine all of the averages together
        summary_df <- do.call(rbind, all_summaries) %>%
            column_to_rownames('gene_id')
        head(summary_df)

        return(summary_df)
    }

    if (summary_type == 'pcs'){
        head(all_summaries[[1]])
        names(all_summaries[[1]])

        # combine the pcs into one df
        head(all_summaries[[1]]$formatted_pcs)
        all_pcs <- do.call(rbind, lapply(all_summaries, function(x) x$formatted_pcs))
        head(all_pcs)

        # combine loadings into one df
        head(all_summaries[[1]]$loadings)
        genes <- lapply(all_summaries, function(x) x$gene_id) %>%
            unlist()
        head(genes)
        all_loadings <- lapply(all_summaries, function(x) x$loadings) 
        names(all_loadings) <- genes
        head(all_loadings)

        # combine other info into one df
        all_info <- do.call(rbind, lapply(all_summaries, function(x) x$pc_info))
        head(all_info)

        return(list(pcs=all_pcs,
                    loadings=all_loadings,
                    pc_info=all_info))

    }



    1
}


combine_results <- function(){
    region_types <- c(
                      'promoters', 'exons', 'introns',
                      '1to5kb', '5UTRs', '3UTRs',
                      'cpg_islands', 'cpg_shores',
                      'cpg_shelves', 'full_gene',
                      'astro_enh', 'gene_body', 'preTSS'
                        )
    region_idx <- 12
    region_types <- c(region_types[region_idx])
    # for testing
    #region_types <- c('introns')
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')

    #cell_types <- c('astro')
    region_types

    #region_idx <- 1
    #celltype_idx <- 1

    #region_type <- region_types[region_idx]
    #cell_type <- cell_types[celltype_idx]

    
    total_jobs <- 100
    num_jobs = total_jobs

    #process_celltype(cell_type)

    runs <- expand.grid(cell_type=cell_types,
                        region_type=region_types
    )
    head(runs)

    array_idx = 1
    run <- runs[1,]
    process_celltype <- function(run){
        region_type <- run[['region_type']]
        cell_type <- run[['cell_type']]
            
        print(paste("Processing", region_type, cell_type, Sys.time()))
        #pc_dir <- paste0(savedir, cleaned_str, region_type, "_", cell_type, "_summaries/")
        #pc_dir

        #list.files(pc_dir)

        ## -- load pcs
        load_pcs <- function(array_idx){
            #savefile <- paste0(pc_dir, cleaned_str, region_type, "_", cell_type, "_gene_pcs_", array_idx, ".rds")
            #savefile
            summary_type <- "pcs"
            savepath <- paste(summary_type, region_type, cell_type, sep='_')
            pc_dir <- paste0(savedir, savepath, "/")
            savefile <- paste0(pc_dir, savepath, array_idx, "_summaries.rds")
            savefile

            if (!file.exists(savefile)){
                print(paste("Job file", array_idx, "does not exists"))
                return(NA)
            }

            pcs <- readRDS(savefile)
            head(pcs)

            pcs
        }

        all_pcs <- lapply(0:(num_jobs-1), load_pcs)

        names(all_pcs[[1]])

        pcs_df <- do.call(rbind, lapply(all_pcs, function (x) x$pcs))
        head(pcs_df)
        dim(pcs_df)

        savefile <- paste0(savedir, cleaned_str, region_type, "_", cell_type, "_pcs.rds")
        savefile
        saveRDS(pcs_df, savefile)


        # do the same for loadings
        type(all_pcs[[1]]$loadings)
        all_loadings <- do.call(c, lapply(all_pcs, function(x) x$loadings))
        head(all_loadings)
        length(all_loadings)

        savefile <- paste0(savedir, cleaned_str, region_type, "_", cell_type, "_pc_loadings.rds")
        savefile
        saveRDS(all_loadings, savefile)

        # do the same for pc info
        head(all_pcs[[1]]$pc_info)
        all_pc_info <- do.call(rbind, lapply(all_pcs, function(x) x$pc_info))
        head(all_pc_info)
        dim(all_pc_info)

        savefile <- paste0(savedir, cleaned_str, region_type, "_", cell_type, "_pc_info.rds")
        savefile
        saveRDS(all_pc_info, savefile)


        # -- load avgs
        load_avgs <- function(array_idx){
            summary_type <- "avgs"
            savepath <- paste(summary_type, region_type, cell_type, sep='_')
            avg_dir <- paste0(savedir, savepath, "/")
            savefile <- paste0(avg_dir, savepath, array_idx, "_summaries.rds")
            savefile
            #savefile <- paste0(pc_dir, cleaned_str, region_type, "_", cell_type, "_gene_avgs_", array_idx, ".rds")
            #savefile

            if (!file.exists(savefile)){
                print(paste("Avgs file", array_idx, "does not exists"))
                return(NA)
            }

            avgs <- readRDS(savefile)
            head(avgs)

            avgs
        }

        all_avgs <- lapply(0:(num_jobs-1), load_avgs)
        avg_df <- do.call(rbind, all_avgs) 
        head(avg_df)
        dim(avg_df)

        savefile <- paste0(savedir, cleaned_str, region_type, "_", cell_type, "_avgs.rds")
        savefile
        saveRDS(avg_df, savefile)
        1
    }

    runs
    process_celltype(runs[5,])
    apply(runs, 1, process_celltype)

    1
}


#mvals <- cleaned_mvals
compute_pcs <- function(mvals, gene_map, region_type, array_idx, cell_type, cleaned_str){
    # summarise the regions using averages and PCs

    head(mvals)
    head(gene_map)


    pc_dir <- paste0(savedir, cleaned_str, region_type, cell_type, "_pcs/")
    savefile <- paste0(pc_dir, cleaned_str, region_type, cell_type, "_gene_pcs_", array_idx, ".rds")
    if (file.exists(savefile)) return(NA)


    # mvals are NOT centered and scaled prior
    # to computing averages

    # compute region averages
    head(mvals)
    region_mvals <- mvals %>%
        as.data.frame() %>%
        rownames_to_column('cpg_id') %>%
        select(cpg_id) 
    head(region_mvals)

    head(gene_map)

    region_map <- gene_map %>% 
        filter(cpg_id %in% region_mvals$cpg_id)
    dim(region_map)
    dim(gene_map)


    # figure out which genes to process
    total_jobs <- 800 # for testing
    #total_jobs <- 100 # seems reasonable - should take about 
    genes <- unique(region_map$gene_id)
    num_genes <- length(genes)
    num_genes

    genes_per_job <- ceiling(num_genes/total_jobs)
    genes_per_job

    start_idx <- array_idx * genes_per_job + 1
    end_idx <- min((array_idx+1) * genes_per_job, num_genes) # don't go past the end of the gene array

    print(paste("Job", array_idx, "processing genes", start_idx, "to", end_idx, "out of", num_genes))

    job_genes <- genes[start_idx:end_idx]
    length(job_genes)


    # for testing
    gene <- genes[3]
    gene_cpgs <- region_map %>%
        filter(gene_id == gene)
    gene_cpgs

    # compute region PCs
    ordered_samples <- colnames(mvals)
    get_gene_pcs <- function(gene){

        print(paste("Processing gene", which(gene == genes), "out of", length(job_genes), Sys.time()))

        gene_cpgs <- region_map %>%
            filter(gene_id == gene) 
        gene_cpgs
        dim(gene_cpgs)

        # grab the mvals for this gene
        gene_mvals <- as.data.frame(mvals)[gene_cpgs$cpg_id,,drop=FALSE]
        head(gene_mvals)
        dim(gene_mvals)

        num_cpgs <- nrow(gene_mvals)

        # get the pcs for this gene
        pc_res <- get_pcs(gene_mvals)
        names(pc_res)

        # no error occurred
        if (all(is.na(pc_res$error))){
            gene_pcs <- pc_res$sig_pcs
            head(gene_pcs)

            # format the pcs to combine later
            formatted_pcs <- gene_pcs[ordered_samples,,drop=FALSE] %>% 
                t() %>%
                as.data.frame()
            rownames(formatted_pcs) <- paste(gene, rownames(formatted_pcs), sep='-')
            head(formatted_pcs)

            # add number of cpgs and estiamted dims
            formatted_pcs <- formatted_pcs %>%
                mutate(num_cpgs=num_cpgs,
                       est_dim=pc_res$est_dim,
                       error=NA
                )
            
        } else{ # error occurred
            formatted_pcs = data.frame(row.names=c(ordered_samples, 'num_cpgs', 'est_dim', 'error'), 
                                       gene=rep(NA, length(ordered_samples)+3)) %>% t() %>% as.data.frame()
            formatted_pcs$error = as.character(pc_res$error)
            formatted_pcs$num_cpgs = num_cpgs
            rownames(formatted_pcs) = gene
            head(formatted_pcs)
        }

        formatted_pcs
    }

    # run for all of the job's genes
    all_pcs <- lapply(job_genes, get_gene_pcs)

    # combine data together
    pcs_df <- do.call(rbind, all_pcs) 
    head(pcs_df)
    dim(pcs_df)

    # to save job pcs
    pc_dir <- paste0(savedir, cleaned_str, region_type, cell_type, "_pcs/")
    pc_dir
    dir.create(pc_dir, showWarnings=FALSE)

    savefile <- paste0(pc_dir, cleaned_str, region_type, cell_type, "_gene_pcs_", array_idx, ".rds")
    savefile
    saveRDS(pcs_df, savefile)

    1
}


map_cpgs <- function(lifted_mvals, annots, region_type, cell_type){
    # map cpgs to the annotation regions

    # save the gene map
    savefile <- paste0(savedir, region_type, "_", cell_type, "_gene_cpg_map.rds")
    savefile
    if (file.exists(savefile)){
        gene_map <- readRDS(savefile)
        head(gene_map)
        return(gene_map)
    }


    #names(lifted_mvals)

    mvals <- lifted_mvals
    head(mvals)



    # create gr annots and for cpgs
    cpgs <- data.frame(cpgs=(rownames(mvals))) %>%
        separate(cpgs, c('Chromosome', 'Start', 'End'), sep='_', remove=FALSE)
    head(cpgs)

    cpgs_gr <- makeGRangesFromDataFrame(cpgs)
    head(cpgs_gr)


    head(annots)
    annots_gr <- makeGRangesFromDataFrame(annots, keep.extra.columns=TRUE)
    head(annots_gr)


    # get overlap of cpgs with regions
    overlaps <- findOverlaps(annots_gr, cpgs_gr, ignore.strand=TRUE) %>% as.data.frame()
    head(overlaps)

    head(annots)
    head(cpgs)

    # connect overlaps together
    gene_map <- data.frame(gene_id=annots[overlaps$queryHits,]$gene_id,
                            cpg_id=cpgs[overlaps$subjectHits,]$cpgs)
    head(gene_map)

    # save the gene map
    savefile <- paste0(savedir, region_type, "_", cell_type, "_gene_cpg_map.rds")
    savefile
    saveRDS(gene_map, savefile)

    gene_map_ext <- data.frame(enh=annots[overlaps$queryHits,],
                            cpg=cpgs[overlaps$subjectHits,])
    head(gene_map_ext)

    savefile <- paste0(savedir, region_type, "_", cell_type, "_gene_cpg_map_extended.rds")
    savefile
    saveRDS(gene_map_ext, savefile)

    gene_map
}

# -- NOT IN USE -- USE NEXT DMP SCRIPT!!! --- 
consolidate_pcs <- function(){
    # do not call this function 
    # just run on interactive session to consolidate
    # results from all of the jobs
    region_idx = 2
    region_types <- c('genes', 'promoters', 'expanded_promoter', 'enhancers')
    region_type <- region_types[region_idx]

    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')
    #cleaned_strs <- c("", "cleaned_")

    cell_type <- cell_types[celltype_idx]


    num_jobs = 100


    # load the data and combine
    load_data <- function(array_idx){
        datafile <- paste0(pc_dir, cleaned_str, region_type, cell_type, "_pcs.rds")
        datafile

        if (!file.exists(datafile)) return(NA)

        pcs <- readRDS(datafile)
        head(pcs)

        pcs
    }


    # -- mvals before cleaning first
    cleaned_str = ""
    "../../output/episcore/03_summarised_genes/includeCovariates/genesastro_pcs"
    pc_dir <- paste0(savedir, cleaned_str, region_type, "_pcs/")
    all_pcs <- lapply(1:num_jobs, load_data)
    pcs_df <- do.call(rbind, all_pcs)
    head(pcs_df)
    dim(pcs_df)

    savefile <- paste0(savedir, cleaned_str, region_type, "_combined_pcs.rds")
    savefile
    saveRDS(pcs_df, savefile)


    # -- mvals after cleaning 
    cleaned_str = "cleaned_"
    pc_dir <- paste0(savedir, cleaned_str, region_type, "_pcs/")
    all_pcs <- lapply(1:num_jobs, load_data)
    pcs_df <- do.call(rbind, all_pcs)
    head(pcs_df)
    dim(pcs_df)

    savefile <- paste0(savedir, cleaned_str, region_type, "_combined_pcs.rds")
    savefile
    saveRDS(pcs_df, savefile)
}

match_geno <- function(lifted_mvals){
    # keep samples from geno only
    datafile <- paste0("/home/eulalio/deconvolution/rosmap/output/qtl_analysis/09_remove_outliers/all_chroms_hg38.fam")
    geno_fam <- read_tsv(datafile, col_names=FALSE)
    head(geno_fam)

    keep_samples <- unique(geno_fam$X2)
    head(keep_samples)
    length(keep_samples)

    head(lifted_mvals)
    dim(lifted_mvals)


    # get the mapping of sample names
    datafile <- paste0("../../output/02_deconvolution/02_svd_batch_correction/restimate_proportions/c1_c2_covariates/astro_formatted_clinical.rds")
    clinical <- readRDS(datafile)
    head(clinical)

    filtered_clinical <- clinical %>%
        filter(individualID %in% colnames(lifted_mvals)) %>%
        filter(wgs.specimenID %in% keep_samples)
    head(filtered_clinical)
    dim(filtered_clinical)

    length(intersect(filtered_clinical$individualID, colnames(lifted_mvals)))

    filtered_mvals <- lifted_mvals[,filtered_clinical$individualID]
    head(filtered_mvals)
    dim(filtered_mvals)

    filtered_mvals
}

save_results <- function(gene_summaries, summary_type, region_type, cell_type, array_idx){ 
    print(paste("Saving results for", summary_type, region_type, cell_type, array_idx))
    head(gene_summaries)
    savepath <- paste(summary_type, region_type, cell_type, sep='_')
    summarydir <- paste0(savedir, savepath, "/")
    summarydir
    dir.create(summarydir)
    savefile <- paste0(summarydir, savepath, array_idx, "_summaries.rds")
    savefile

    saveRDS(gene_summaries, savefile)
    1
}


main <- function(array_idx, region_idx, celltype_idx){
    # load methylation data 
    cell_types <- c('astro', 'endo', 'neuron', 'oligo_opc', 'bulk')
    #cell_types <- c('endo', 'neuron', 'oligo_opc', 'bulk')
    region_types <- c(
                      'promoters', 'exons', 'introns',
                      '1to5kb', '5UTRs', '3UTRs',
                      'cpg_islands', 'cpg_shores',
                      'cpg_shelves', 'full_gene',
                      'astro_enh', 'gene_body', 'preTSS'
                        )

    # for testing
    testing=FALSE
    if (testing){
        cell_type <- cell_types[1]
        region_idx <- 13
    }

    process_celltype <- function(cell_type){ 

        print(paste("Loading data for", cell_type, Sys.time()))

        # lift cpgs to hg38 if we didn't do this yet
        savefile <- paste0(savedir, precleaned_str, cleaned_str, cell_type, "_lifted_hg38_mvals.rds")
        savefile
        if (file.exists(savefile)){
            print(paste("LOADING EXISTING LIFTED MVALS"))
            lifted_mvals <- readRDS(savefile)
        } else{
            mvals <- load_methylation(cell_type)

            # convert cpgs to hg38 from hg19
            lifted_mvals <- lift_methylation(mvals, cell_type)
            head(lifted_mvals)

            savefile
            saveRDS(lifted_mvals, savefile)
        }

        dim(lifted_mvals)

        if (randomize_cpgs){
            print(paste("Randomizing CpG order"))
            # permute the CpG labels
            head(lifted_mvals)[,1:10]

            set.seed(1174117171)
            permuted_cpgs <- rownames(lifted_mvals)[sample(1:nrow(lifted_mvals), nrow(lifted_mvals), replace=FALSE)]
            head(permuted_cpgs)
            head(rownames(lifted_mvals))

            rownames(lifted_mvals) <- permuted_cpgs
            head(rownames(lifted_mvals))
        }



        # assign the region type here
        region_type <- region_types[region_idx]

        # create rnbeads object and get annotations
        print(paste("Loading annotations for", region_type, Sys.time()))
        annots <- get_annotations(region_type)
        head(annots)
        print(paste("Number of genes in", region_type, length(unique(annots$gene_id))))


        # map cpgs to regions
        print(paste("Mapping cpgs", Sys.time()))
        gene_map <- map_cpgs(lifted_mvals, annots, region_type, cell_type)
        head(gene_map)


        # select the mvals to use
        mvals <- lifted_mvals
        head(mvals)
        dim(mvals)


        # filter down map to the genes that we will process
        job_map <- select_genes(array_idx, mvals, gene_map, region_type, cell_type, cleaned_str)

        print(paste("Summarising genes", Sys.time()))
        summary_type = 'avgs'
        gene_summaries <- summarise_genes(mvals, job_map, summary_type)
        # save results
        save_results(gene_summaries, summary_type, region_type, cell_type, array_idx)


        summary_type = 'pcs'
        gene_summaries <- summarise_genes(mvals, job_map, summary_type)
        # save results
        save_results(gene_summaries, summary_type, region_type, cell_type, array_idx)
        1
    }

    res <- lapply(cell_types, process_celltype)
    #process_celltype(cell_types[celltype_idx])

}

# for testing
#region_types
#cell_types

region_idx = 13
array_idx = 0
celltype_idx = 5
total_jobs = 100


# use this to run from sbatch
args = commandArgs(trailingOnly=TRUE)
print(args)
region_idx <- as.numeric(args[1])
array_idx <- as.numeric(args[2])
celltype_idx <- as.numeric(args[3])
total_jobs <- as.numeric(args[4])

if (summarize_results){
    print(paste("Summarizing results"))
    combine_results()
    quit()
}

#for (array_idx in 1:99){
    #print(paste("Processing", region_types[region_idx]))
    main(array_idx, region_idx, celltype_idx)
#}


quit()

# run jobs
combine_results 


