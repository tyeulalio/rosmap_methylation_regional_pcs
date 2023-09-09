#library(isva) # required for smartSVA
library(wateRmelon)
library(EpiSCORE)
library(TCA)
library(RNOmni)
library(RColorBrewer)
library(ggplot2)
#library(limma)
library(tidyverse)

source("/home/eulalio/deconvolution/new_rosmap/scripts/shared_functions/02_plot_colors.R")

#library(RhpcBLASctl)
#blas_set_num_threads(16) #limit number of threads for openblas to 16

#library(parallel)

# Run this on betas from the normalized and QC'ed beta values

sub_chroms = FALSE
keep_chroms = paste0("chr", 22)

c1_covariates = FALSE
c2_covariates = TRUE
reestimate_props = TRUE

main_savedir <- paste0("../../output/02_deconvolution/01_deconvolved_data/")
main_savedir
dir.create(main_savedir, showWarnings=FALSE)


if (reestimate_props){
    main_savedir <- paste0(main_savedir, "restimate_proportions/")
    dir.create(main_savedir, showWarnings=FALSE)
}


if (sub_chroms){
    chr_savedir <- paste0(main_savedir, "sub_chroms/")
    chr_savedir
    dir.create(chr_savedir)
} else{
    chr_savedir <- main_savedir
}

if (c2_covariates & !c1_covariates){
    savedir <- paste0(chr_savedir, "c2_covariates/")
    savedir
    dir.create(savedir)
}

if (c2_covariates & c1_covariates){
    savedir <- paste0(chr_savedir, "c1_c2_covariates/")
    savedir
    dir.create(savedir)
}

if (!c1_covariates & !c2_covariates){
    savedir <- paste0(chr_savedir, "no_covariates/")
    savedir
    dir.create(savedir)
}

if (c1_covariates & !c2_covariates){
    savedir <- paste0(chr_savedir, "c1_covariates/")
    savedir
    dir.create(savedir)
}

savedir

#precleaned_str = "precleaned_"
precleaned_str = ""

include_study = TRUE
if (include_study){
    savedir = paste0(savedir, "include_study/")
    dir.create(savedir)
}



load_methylation <- function(thresh){

    if (precleaned_str == "precleaned_"){
        print("Loading precleaned data")
        # use the cleaned methylation data
        datafile <- paste0("../../output/01_preprocessing/02_svd_batch_correction/cleaned_betas.rds")
    } else{
        datafile <- paste0("../../output/01_preprocessing/methylation/01_preprocessing_methylation/snp", thresh, "_window0_filtered_betas.rds")
    }

    betas <- readRDS(datafile)
    head(betas)
    dim(betas)

    betas
}

load_clinical <- function(){
    # load phenotype data
    #datafile <- paste0("../../output/01_preprocessing/methylation/01_preprocessing_methylation/snp0.01_window10_filtered_targets.rds")
    datafile <- paste0("../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
    targets <- readRDS(datafile)
    head(targets)


    #datafile <- paste0("../../output/06_svd_batch_correction/formatted_clinical.rds")
    #datafile <- paste0("../../output/01_preprocessing/")
    #datafile <- paste0("../../output/05_check_covariates/clinical_metadata.rds") 
    #clinical <- readRDS(datafile)
    #head(clinical)

    #return(targets)

    # grab map between specimenID and individualID
    #datafile <- paste0("../../output/02_connect_samples/exp_meth_overlapped_meta.rds")
    #meta <- readRDS(datafile)
    #names(meta)
    #meth_meta <- meta$meth_meta
    #head(meth_meta)

    # grab mapping from meth colname and specimenID
    datafile <- paste0("../../data/ROSMAP_data/Epigenetics/Epigenetics (DNA methylation array)/SYNAPSE_METADATA_MANIFEST.tsv")
    manifest <- read_tsv(datafile)
    head(manifest)

    # map colname to specimenID
    name_map <- manifest %>%
        filter(fileFormat == 'idat') %>%
        select(specimenID, name) %>%
        mutate(name = gsub("_[A-Za-z]{3}\\.idat", "", name)) %>%
        unique() %>%
        filter(!is.na(specimenID))
    head(name_map)

    # map specimen to individualID
    #head(meth_meta)
    specimen_map <- targets %>%
        select(individualID, specimenID=meth.specimenID) %>%
        inner_join(name_map) %>%
        rename(meth.specimenID=specimenID)
    head(specimen_map)

    # connect to clincal
    mapped_clinical <- targets %>%
        left_join(specimen_map)
    head(mapped_clinical)


    # connect batch info to metadata
    #datafile <- paste0("../../data/ROSMAP_data/Metadata/ROSMAP_assay_methylationArray_metadata.csv")
    #meta <- read_csv(datafile)
    #head(meta)

    #table(mapped_clinical$Sample_Plate, useNA='ifany')

    #sub_meta <- meta %>%
        #select(specimenID, platform, Sentrix_ID, Sentrix_Row_Column, Sample_Plate, Sample_Well, batch)

    #mapped_clinical <- mapped_clinical %>%
        #left_join(sub_meta)
    #head(mapped_clinical)

    savefile <- paste0(savedir, "mapped_clinical.rds")
    savefile
    saveRDS(mapped_clinical, savefile)

    dim(mapped_clinical)

    mapped_clinical
}

format_data <- function(betas, mapped_clinical){
    # format the betas data and clinical data
    head(betas)
    dim(betas)
    head(mapped_clinical)

    if (precleaned_str == "precleaned_"){
        datafile <- paste0("../../output/01_preprocessing/02_svd_batch_correction/formatted_clinical.rds")
        formatted_clinical <- readRDS(datafile)
        formatted_betas <- betas
        formatted_dat <- list(formatted_betas=formatted_betas, formatted_clinical=formatted_clinical)
        return(formatted_dat)
    }

    if (sub_chroms){
        # subset to chrom 22
        datafile <- paste0("../../output/episcore/03_summarised_genes/includeCovariates2/bulk_cpg_position_map.rds")
        cpg_map <- readRDS(datafile)
        head(cpg_map)

        chr22_map <- cpg_map %>%
            filter(seqnames %in% keep_chroms)
        dim(chr22_map)
        head(chr22_map)

        chr22_betas <- betas[chr22_map$group_name,]
        head(chr22_betas)
        sum(is.na(chr22_betas))

        betas <- chr22_betas
    }

    dim(betas)
    #rmeans <- rowMeans(betas)
    #summary(rmeans)

    #lower_thresh <- 0.05
    #upper_thresh <- 0.95

    #sum(rmeans<lower_thresh)
    #sum(rmeans>upper_thresh)

    ## remove sites with rowmeans > 0.9 and <0.1 as per TCA paper
    #filtered_betas <- betas[(rmeans > lower_thresh) & (rmeans < upper_thresh),]
    #dim(filtered_betas)
    #dim(betas)

    #summary(rowMeans(filtered_betas))


    #betas <- filtered_betas

    # put both data in same order
    head(betas)
    head(mapped_clinical)

    overlapped_names <- intersect(colnames(betas), mapped_clinical$name) %>% na.omit()
    length(overlapped_names)

    ordered_clinical <- mapped_clinical[match(overlapped_names, mapped_clinical$name),]
    ordered_betas <- betas[,overlapped_names]


    if (!identical(colnames(ordered_betas), ordered_clinical$name)) stop("ERROR: ORDERING DOES NOT MATCH")
    

    formatted_betas <- ordered_betas
    colnames(formatted_betas) <- ordered_clinical$individualID
    head(formatted_betas)

    head(ordered_clinical)
    names(ordered_clinical)



    # determine relationship between covariates
    factor_covariates <- c('sex', 'race', 'spanish', 'cogdx', 'diagnosis', 'apoe4', 'meth.batch', 
                           'meth.platform', 'meth.Sentrix_ID', 'meth.Sentrix_Row_Column', 'meth.Sample_Plate', 'meth.Sample_Well')
    cont_covariates <- c('age_death', 'pmi', 'educ')

    ordered_clinical <- data.frame(lapply(ordered_clinical,function(x){x <- sapply(x,function(y){str_replace_all(as.character(y),'\\+','')})}))
    rownames(ordered_clinical) <- ordered_clinical$individualID
    head(ordered_clinical)

    # Convert factor covariates to factors
    ordered_clinical[,factor_covariates] = lapply(ordered_clinical[,factor_covariates], factor)
    ordered_clinical[,cont_covariates] = lapply(ordered_clinical[,cont_covariates], as.character)
    ordered_clinical[,cont_covariates] = lapply(ordered_clinical[,cont_covariates], as.numeric)

    head(ordered_clinical)

    formatted_clinical <- ordered_clinical

    savefile <- paste0(savedir, "formatted_betas.rds")
    savefile
    saveRDS(formatted_betas, savefile)

    formatted_dat <- list(formatted_betas=formatted_betas, formatted_clinical=formatted_clinical)
    formatted_dat
}


get_cell_proportions <- function(formatted_betas, formatted_clinical, thresh){
    # get the proportions of cell types in the data
    head(formatted_betas)

    # use the noob preprocessed data
    load_noob <- function(){
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/01_preprocessing/methylation/01_preprocessing_methylation/noob_preprocessed_meth.rds")
        if (precleaned_str != ""){
            datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/01_preprocessing/02_svd_batch_correction/cleaned_betas.rds")
        }
        datafile
        noob_meth <- readRDS(datafile)
        head(noob_meth)

        # adjust the sample names
        head(formatted_clinical)
        ordered_samples <- intersect(colnames(noob_meth), formatted_clinical$name)
        #ordered_samples <- intersect(colnames(noob_meth), formatted_clinical$individualID)
        length(ordered_samples)

        # match the order of betas and clinical
        ordered_meth <- noob_meth[,ordered_samples]
        ordered_clinical <- formatted_clinical[match(ordered_samples, formatted_clinical$name),]
        #ordered_clinical <- formatted_clinical[match(ordered_samples, formatted_clinical$individualID),]

        identical(colnames(ordered_meth), ordered_clinical$name)
        #identical(colnames(ordered_meth), ordered_clinical$individualID)

        colnames(ordered_meth) <- ordered_clinical$individualID
        head(ordered_meth)

        ordered_meth
    }
    ordered_meth <- load_noob()
    head(ordered_meth)


    # load reference matrix for brain
    data(BrainRef) # should load exprefBrain.m and mrefBrain.m
    head(mrefBrain.m)
    head(exprefBrain.m)

    brain_ref <- mrefBrain.m
    head(brain_ref)

    # summarise betas 
    # computes average beta across cpgs 200bp upstream of TSS
    # cols = samples
    # rows = cpgs
    use_formatted <- function(){
        #head(formatted_betas)

        # need to get probe names again
        datafile <- paste0("../../output/01_preprocessing/methylation/01_preprocessing_methylation/probe_position_map_hg37.rds")
        pos_map <- readRDS(datafile)
        head(pos_map)

        head(formatted_betas)
        ordered_meth <- formatted_betas %>%
            as.data.frame() %>%
            rownames_to_column('cpg') %>%
            left_join(pos_map[c('cpg', 'Name')]) %>%
            column_to_rownames('Name') %>%
            select(-cpg)
        head(ordered_meth)
        ordered_meth
    }
    #ordered_meth <- use_formatted()

   #inf <- probeInfo450k.lv
    #str(inf)
    #inf850 <- probeInfo850k.lv
    #str(inf850)

    #cpg <- "cg00000029"
    #cpg %in% inf$probeID
    #cpg %in% inf850$probeID

    #idx = which(cpg == inf$probeID)
    #idx850 = which(cpg == inf850$probeID)

    #inf$EID[idx]
    #inf850$EID[idx850]

    #length(intersect(rownames(pos_betas), inf$probeID))
    #length(intersect(rownames(pos_betas), inf850$probeID))

    #pos_betas <- formatted_betas

    head(ordered_meth )
    # try remove low variance probes here
    sds <- rowSds(as.matrix(ordered_meth))
    summary(sds)

    sd_thresh = 0
    table(sds >= sd_thresh)

    var_betas <- ordered_meth[sds >= sd_thresh,]

    betas_mat <- as.matrix(var_betas)
    dim(betas_mat)
    head(betas_mat)

    arraytype = "850k"
    avg_betas = constAvBetaTSS(betas_mat, type=arraytype)
    #dim(pos_betas)
    head(avg_betas)


    # decompose data
    head(brain_ref)
    celltype_proportions <- wRPC(avg_betas, brain_ref, maxit=500)

    
    names(celltype_proportions)
    head(celltype_proportions$estF)
    head(celltype_proportions$ref)

    savefile <- paste0(savedir, precleaned_str, "snpThresh", thresh, "_episcore_estimated_celltypes.rds")
    savefile 
    saveRDS(celltype_proportions$estF, savefile)

    savefile <- paste0(savedir, precleaned_str, "snpThresh", thresh, "_episcore_reference_panel.rds")
    savefile 
    saveRDS(celltype_proportions$ref, savefile)


    plot_props <- function(){
        plotdir <- paste0(savedir, "plots/")
        dir.create(plotdir, showWarnings=FALSE)

        # plot the estiamted props
        episcore_props <- celltype_proportions$estF
        head(episcore_props)
        colnames(episcore_props) <- c('neuron', 'oligo', 'astro', 'OPC', 'endo', 'microglia')
        p <- episcore_props %>%
            as.data.frame() %>%
            rownames_to_column('sample_id') %>%
            mutate(oligo_opc = oligo + OPC) %>%
            gather('cell_type', 'episcore', -sample_id) %>%
            ggplot() +
            geom_boxplot(aes(x=cell_type, y=episcore, fill=cell_type),
                         show.legend=FALSE) +
            ggtitle(paste("Cell proportion comparison, thresh =", thresh, ", sd thresh =", sd_thresh, arraytype)) +
            theme_bw()
        savefile <- paste0(plotdir, precleaned_str, "snpThresh", thresh, "_episcore_boxplot.png")
        savefile
        ggsave(plot=p, file=savefile)

        # compare with neun proportions
        datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/01_preprocessing/methylation/01_preprocessing_methylation/neun_estimated_cell_counts.rds")
        neun_props <- readRDS(datafile)
        head(neun_props)

        # change sample names
        head(formatted_clinical)
        formatted_neun <- neun_props %>%
            as.data.frame() %>%
            rownames_to_column('name') %>%
            left_join(formatted_clinical[c('name', 'individualID')]) %>%
            select(-name)
        head(formatted_neun)

        head(episcore_props)

        joined_props <- episcore_props %>%
            as.data.frame() %>%
            rownames_to_column('individualID') %>%
            full_join(formatted_neun)
        head(joined_props)

        joined_props %>%
            filter(is.na(neuron)) %>%
            head()

        p <- ggplot(joined_props) +
            geom_point(aes(x=neuron, y=NeuN_pos)) +
            geom_abline(slope=1) +
            coord_cartesian(xlim=c(0,0.6), ylim=c(0,0.6)) +
            ylab(paste("NeuN estimated cell proportion")) +
            xlab(paste("Episcore estimated cell proportion")) +
            ggtitle(paste("Cell proportion comparison, thresh =", thresh, ", sd thresh =", sd_thresh, arraytype))
        savefile <- paste0(plotdir, "cell_proportions_comparions_scatter.png")
        savefile
        ggsave(p, file=savefile)


        # plot the proportion neurons with braak score
        head(joined_props)
        head(formatted_clinical)

        #joined_clin <- joined_props %>%
            #left_join(formatted_clinical[c('individualID', 'braaksc', 'diagnosis', 'cogdx', 'ceradsc', 'age_death')])
            ##gather('method', 'neuron_proportion', NeuN_pos, neuron)
        #head(joined_clin)

        #trait <- "braaksc"
        #res <- lm(joined_clin$neuron~ joined_clin$age_death)
        #summary(res)

        #trait <- 'diagnosis'
        #joined_clin['trait'] <- joined_clin[trait]
        #p <- ggplot(joined_clin) +
            #geom_boxplot(aes(x=as.factor(trait), y=neuron_proportion)) +
            #facet_wrap(method ~ .) +
            #ggtitle(paste("trait~neuron prop", trait)) +
            #xlab(trait)
        #savefile <- paste0(plotdir, trait, "_trait_neuron_boxplot.png")
        #ggsave(p, file=savefile)

        
    }
    plot_props()

    celltype_proportions$estF
}

deconvolve_data <- function(formatted_betas, proportions_df, formatted_clinical){
    head(proportions_df)
    summary(proportions_df$Oligo) # very small, combine with oligo
    summary(proportions_df$OPC) # very small, combine with oligo
    summary(proportions_df$Microglia) # zero microliga -- remove

    proportions_df <- proportions_df %>%
        mutate(oligo_opc = Oligo + OPC)
    head(proportions_df)

    proportions_df %>%
        filter(oligo_opc != OPC) %>%
        head() # only 4 samples with Oligo

    sub_props <- proportions_df %>%
        select(neuron, oligo_opc, astro=Astro, endo=Endo)
    head(sub_props)

    savefile <- paste0(savedir, precleaned_str, "formatted_celltype_proportions.rds")
    saveRDS(sub_props, savefile)

    head(formatted_betas)
    head(sub_props)



    # decide whether to add batch in the tca model
    # c1 = covariates taht affect methylation at the cell type level (age, gender)
    # c2 = covariates that are not expected to affect meth at the celll type level,
    # but rather may affect tissue-level meth (batch, technical variation)
    if (c1_covariates | c2_covariates){
        head(formatted_clinical)
        identical(colnames(formatted_betas), formatted_clinical$individualID)

        sub_clinical <- formatted_clinical %>%
            filter(!is.na(age_death)) %>%
            mutate(study_num = as.numeric(as.factor(Study)))

        table(sub_clinical$study_num)

        c1_vars = c('age_death', 'msex')
        c2_vars = c('meth.batch', 'pmi')

        if (include_study){
            c2_vars = c('meth.batch', 'pmi', 'study_num')
        }
        
        # set up c1 design matrix
        f = paste0("~ ", paste(c1_vars, collapse = ' + ')) %>% as.formula()
        f
        c1_design <- model.matrix(f, data=sub_clinical)[,-1]
        head(c1_design)

        rownames(c1_design) = sub_clinical$individualID
        colnames(c1_design)


        # set up c2 design matrix
        # replace missing pmi with the mean
        sub_clinical <- sub_clinical %>%
            mutate(pmi = ifelse(is.na(pmi), mean(pmi, na.rm=TRUE), pmi))

        f = paste0("~ ", paste(c2_vars, collapse = ' + ')) %>% as.formula()
        f
        c2_design <- model.matrix(f, data=sub_clinical)[,-1]
        head(c2_design)

        rownames(c2_design) = sub_clinical$individualID
        colnames(c2_design)

        head(c2_design)
        dim(c2_design)


        sub_betas = formatted_betas[,sub_clinical$individualID]
        dim(sub_betas)
        
        sub_props <- sub_props[sub_clinical$individualID,]
        head(sub_props)

    } else{
        sub_betas = formatted_betas
    }


    # run the tca deconvolution
    # betas - cols=samples, rows=cpgs
    # proportions = cols=celltypes, rows=samples
    if (!c1_covariates & c2_covariates){
        print("!c1 and c2 covariates")
        tca_out <- tca(X=sub_betas,
                       W=sub_props,
                       #C1=c1_design,
                       C2=c2_design,
                       #parallel = TRUE, num_cores = 16
                       constrain_mu=TRUE,
                       refit_W=TRUE
        )
    }
    if (c1_covariates & c2_covariates){
        print("c1 and c2 covariates")
        tca_out <- tca(X=sub_betas,
                       W=sub_props,
                       C1=c1_design,
                       C2=c2_design,
                       #parallel = TRUE, num_cores = 16
                       constrain_mu=TRUE,
                       refit_W=TRUE
        )
    }
    if (!c1_covariates & !c2_covariates){
        print("!c1 and !c2 covariates")
        tca_out <- tca(X=sub_betas,
                       W=sub_props,
                       #C1=c1_design,
                       #C2=c2_design
                       #parallel = TRUE, num_cores = 16
                       constrain_mu=TRUE,
                       refit_W=TRUE
        )
    }
    if (c1_covariates & !c2_covariates){
        print("c1 and !c2 covariates")
        tca_out <- tca(X=sub_betas,
                       W=sub_props,
                       C1=c1_design,
                       #C2=c2_design
                       #parallel = TRUE, num_cores = 16
                       constrain_mu=TRUE,
                       refit_W=TRUE
        )
    }

    savefile <- paste0(savedir, precleaned_str, "tca_out.rds")
    savefile
    #if (include_batch == TRUE){
        #savefile <- paste0(savedir, "batch_included_tca_out.rds")
    #}

    saveRDS(tca_out, savefile)


    # project from original mixed betas matrix to celltype-specific tensor
    # this step uses multilinear algebra
    head(sub_betas)
    betas_mat <- as.matrix(sub_betas)
    tca_tensor <- tensor(betas_mat, tca_out, parallel=FALSE)

    test <- tca_tensor[[1]]

    savefile <- paste0(savedir, precleaned_str, "tca_tensor.rds")
    savefile
    saveRDS(tca_tensor, savefile)

    1
}

plot_proportions <- function(thresh, formatted_clinical){
    # compare the original and the reestimated proportions

    # load episcore proportions
    #datafile <- paste0("../../output/episcore/01_deconvolved_data/includeCovariates2/", "episcore_estimated_celltypes.rds")
    #datafile <- paste0(savedir, "episcore_estimated_celltypes.rds")
    thresh <- 0
    datafile <- paste0(savedir, precleaned_str, "snpThresh", thresh, "_episcore_estimated_celltypes.rds")
    episcore_props <- readRDS(datafile)

    # load reestimated props
    datafile <- paste0(savedir, "tca_out.rds")
    tca_out <- readRDS(datafile)
    names(tca_out)

    tca_props <- tca_out$W

    head(tca_props)
    head(episcore_props)

    episcore_df <- episcore_props %>%
        as.data.frame()
    head(episcore_df)
    summary(episcore_df$Oligo)
    table(episcore_df$Oligo > 0)

    identical(rownames(tca_props), rownames(episcore_df))
    head(tca_props)
    head(episcore_df)

    episcore_matched <- episcore_df %>%
        mutate(oligo_opc = Oligo + OPC) %>%
        select(neuron, oligo_opc, astro=Astro, endo=Endo)
    head(episcore_matched)
    head(tca_props)

    corrs <- cor.test(as.matrix(episcore_matched),as.matrix(tca_props))
    corrs$p.value

    # make long versions
    tca_long <- tca_props %>%
        as.data.frame() %>%
        rownames_to_column('sample_id') %>%
        gather('cell_type', 'tca', -sample_id)
    head(tca_long)


    colnames(episcore_props) <- c('neuron', 'oligo', 'astro', 'OPC', 'endo', 'microglia')
    epi_long <- episcore_props %>%
        as.data.frame() %>%
        rownames_to_column('sample_id') %>%
        mutate(oligo_opc = oligo + OPC) %>%
        gather('cell_type', 'episcore', -sample_id)
    head(epi_long)

    props <- inner_join(tca_long, epi_long)
    head(props)

    plotdir <- paste0(savedir, "plots/")
    dir.create(plotdir, showWarnings=FALSE)

    # plot boxplots of proportions
    formatted_props <- props %>%
        left_join(cell_type_map)
    head(formatted_props)

    # plot tca estimates
    p <- ggplot(formatted_props) +
        geom_boxplot(aes(x=formatted_ct, y=tca, 
                         fill=formatted_ct),
                     show.legend=FALSE,
                     linewidth=1.5) +
        scale_fill_manual(values=cell_type_colors) +
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              axis.title.x = element_blank(),
              text=element_text(size=20)
        ) +
        ylab(paste("Estimated cell type proportions"))
    (savefile <- paste0(plotdir, "tca_celltype_proportions_boxplot.png"))
    ggsave(p, file=savefile, width=7, height=7)

    (savefile <- paste0(plotdir, "tca_celltype_proportions_boxplot.rds"))
    saveRDS(p, savefile)

    # plot episcore estimates
    p <- ggplot(formatted_props) +
        geom_boxplot(aes(x=formatted_ct, y=episcore, 
                         fill=formatted_ct),
                     show.legend=FALSE,
                     linewidth=1.5) +
        scale_fill_manual(values=cell_type_colors) +
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border = element_rect(color='black', fill=NA, linewidth=2),
              axis.title.x = element_blank(),
              text=element_text(size=30),
              axis.text.x = element_text(angle=25, hjust=1, vjust=1)
        ) +
        ylab(paste("Estimated proportions"))
    (savefile <- paste0(plotdir, "episcore_celltype_proportions_boxplot.png"))
    ggsave(p, file=savefile, width=7, height=6)


    head(formatted_clinical)
    clin_props <- props %>%
        rename(individualID=sample_id) %>%
        left_join(formatted_clinical)
    head(clin_props)
    savefile <- paste0(savedir, "clinical_props.rds")
    saveRDS(clin_props, savefile)

    str(clin_props)

    # plot the cell type proportions by age
    head(clin_props)
    sub_clin_props <- clin_props %>%
        select(cell_type, tca, episcore, age_death, diagnosis, sex, braaksc:dcfdx_lv) %>%
        mutate(age_brackets = cut(age_death, breaks=seq(60,100, by=5), right=FALSE, include.lowest =TRUE)) %>%
        filter(age_death > 70)
    head(sub_clin_props)
    summary(sub_clin_props$age_death)

    p <- ggplot(sub_clin_props) +
        geom_boxplot(aes(x=age_brackets, y=episcore, fill=cell_type)) +
        facet_wrap(. ~ diagnosis, ncol=1) +
        theme_bw()
    (savefile <- paste0(plotdir, "episcore_celltype_by_age_boxplot.png"))
    ggsave(plot=p, file=savefile, height=7, width=7)

    # check a line plot
    head(sub_clin_props)
    p <- ggplot(sub_clin_props, aes(x=age_death, y=episcore, color=cell_type)) +
        geom_point(alpha=0.5) +
        geom_smooth(se=FALSE, method='lm') +
        facet_grid(sex ~ ceradsc) +
        theme_bw()
    (savefile <- paste0(plotdir, "episcore_celltype_by_age_line.png"))
    ggsave(plot=p, file=savefile, height=7, width=7)

    # plot the estiamted props
    head(episcore_props)
    p <- episcore_props %>%
        as.data.frame() %>%
        rownames_to_column('sample_id') %>%
        mutate(oligo_opc = oligo + OPC) %>%
        gather('cell_type', 'episcore', -sample_id) %>%
        ggplot() +
        geom_boxplot(aes(x=cell_type, y=episcore, fill=cell_type))
    savefile <- paste0(plotdir, "episcore_boxplot.png")
    ggsave(plot=p, file=savefile)

    # plot the reestimated props
    p <- props %>%
        ggplot() +
        geom_boxplot(aes(x=cell_type, y=tca, fill=cell_type),
                     show.legend=FALSE) +
        theme_bw() +
        theme(text=element_text(size=20)) +
        ylab(paste("Estimated proportion")) +
        xlab(paste("Cell type"))
    savefile <- paste0(plotdir, "reestimated_boxplot.png")
    ggsave(plot=p, file=savefile)

    # plot the props together
    p <- ggplot(props) +
        geom_point(aes(x=tca, y=episcore, color=cell_type)) +
        geom_abline(slope=1, color='black') +
        facet_wrap(. ~ cell_type)

    savefile <- paste0(plotdir, "proportions_scatter.png")
    savefile
    ggsave(plot=p, filename=savefile)


    ## -- check the proportion of neurons with age
    head(props)
    head(formatted_clinical)

    sub_clinical <- formatted_clinical %>%
        select(sample_id=individualID, age_death, braaksc, ceradsc, batch=meth.batch, study=Study, msex,
               pmi, diagnosis, cogdx, dcfdx_lv, apoe4
               ) %>%
        mutate(braaksc = as.numeric(braaksc)) %>%
        mutate(ceradsc= as.numeric(ceradsc)) %>%
        mutate(cat_braaksc = as.numeric(braaksc >= 5))
    head(sub_clinical)




    ct <- 'neuron'
    tr <- 'braaksc'
    row <- runs[2,]
    get_corr <- function(row){
        ct <- row[['cell_type']] %>% as.character()
        tr <- row[['trait']] %>% as.character()
        print(paste("Processing", ct, tr))

        head(props)
        str(props)
        props_clinical <- props %>%
            left_join(sub_clinical) %>%
            filter(cell_type == ct) %>%
            mutate(scaled_props = (tca * (nrow(props) - 1) + 0.5) / nrow(props)) %>%
            mutate(logit_tca = log2(scaled_props/(1-scaled_props))) 
        head(props_clinical)
        summary(props$logit_tca)

        #if (tr == 'braaksc'){
            #props_clinical <- props_clinical %>%
                #filter(braaksc > 0)
        #}

        res <- lm(props_clinical$logit_tca ~ 
                  props_clinical[[tr]] +
                  props_clinical$age_death +
                  props_clinical$msex +
                  props_clinical$study +
                  props_clinical$pmi +
                  props_clinical$batch 
                  #props_clinical$apoe4
        )
        summary(res)
        names(res)

        sres <- (summary(res))
        names(sres)
        sres$coefficients

        i <- 2
        data.frame(trait=tr,
                   cell_type=ct,
                   effect_size=res$coefficients[[i]],
                   standard_error=sqrt(diag(vcov(res)))[[i]],
                   pval=sres$coefficients[[i,4]]
        ) %>%
        mutate(sig = ifelse(pval < 0.05, "*", ""))
    }

    
    cell_types <- unique(props$cell_type)
    traits <- c('ceradsc', 'braaksc', 'cat_braaksc')
    runs <- expand.grid(cell_type=cell_types, trait=traits)
    head(runs)
    corrs <- apply(runs, 1, get_corr)
    corrs_df <- do.call(rbind, corrs)
    corrs_df


    # is age and proportions correlated?
    check_age <- function(ct){
        head(props)
        head(sub_clinical)
        num_samples <- length(unique(props$sample_id))
        props_clinical <- props %>%
            left_join(sub_clinical) %>%
            filter(cell_type == ct) %>%
            #mutate(scaled_props = (tca-min(tca)) / (max(tca) - min(tca))) %>%
            mutate(scaled_props = (tca * (num_samples - 1) + 0.5) / num_samples) %>%
            mutate(logit_tca = log2(scaled_props/(1-scaled_props))) %>%
            #mutate(binned_age = cut(age_death, breaks=c(60,70,80,90)))
            mutate(binned_age = ntile(age_death, 20))
            #filter(age_death != 90)
        head(props_clinical)


        summary(props_clinical$tca)
        table(props_clinical$binned_age)
        summary(props_clinical$scaled_props)

        res <- lm(logit_tca ~ 
                  binned_age,
              data=props_clinical
        )
        sres <- summary(res)
        sres
        names(sres)
        sres$coefficients
        sres
    }
    check_age('oligo_opc')
        
    ct <- "neuron"
    props_clinical <- props %>%
        left_join(sub_clinical) %>%
        filter(cell_type == ct) %>%
        mutate(int_tca = RankNorm(tca)) %>%
        mutate(int_age = RankNorm(age_death))
    head(props_clinical)
    summary(props_clinical$int_age)
    summary(props_clinical$int_tca)

    lmfit <- lm(episcore ~ 0 + diagnosis, data=props_clinical)
    summary(lmfit)

        1
}

main <- function(){

    print("Loading data")

    # load normalized methylation
    thresh <- 0
    betas <- load_methylation(thresh)
    dim(betas)

    # load the phenotype data
    mapped_clinical <- load_clinical()

    # format the data 
    formatted_dat <- format_data(betas, mapped_clinical)
    names(formatted_dat)
    formatted_betas <- formatted_dat$formatted_betas
    formatted_clinical <- formatted_dat$formatted_clinical

    head(formatted_clinical)
    dim(formatted_betas)

    savefile <- paste0(savedir, "formatted_clinical.rds")
    saveRDS(formatted_clinical, savefile)

    formatted_clinical <- readRDS(savefile)

    #formatted_betas <- betas

    # get the cell type proportions 
    print(paste("Getting cell type proportions", Sys.time()))
    celltype_proportions <- get_cell_proportions(formatted_betas, formatted_clinical, thresh)


    # deconvolve the data
    print(paste("Deconvolving data", Sys.time()))
    proportions_df <- as.data.frame(celltype_proportions)
    head(proportions_df)

    deconvolve_data(formatted_betas, proportions_df, formatted_clinical)

    plot_proportions(thresh, formatted_clinical)
}

main()
