library(isva) # required for smartSVq
library(uwot)
library(RNOmni)
library(PCAtools)
library(leidenAlg)
library(cluster)
library(bluster)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(GenomicRanges)
library(liftOver)
library(wateRmelon)
library(tidyverse)

# plot the methylation values for cell types in a UMAP together
# use m-values for UMAP plotting

datadir = "/path/to/data"

load_meth <- function(cell_type, downsample){
    # downsample argument allows for testing on smaller datasets

    # load methylation data
    print(paste("Loading methylation for", cell_type))

    # data file for this cell type
    datafile <- paste0(datadir, cell_type, "_formatted_meth.rds")
    datafile

    # read in the methylation data
    # append celltype to column names
    mvals <- readRDS(datafile)
    colnames(mvals) <- paste0(colnames(mvals), "_", cell_type)
    head(mvals)

    if (downsample){
        print(paste("DOWNSAMPLING METH DATA"))
        # just keep 50 people for testing
        sub_mvals <- mvals[,1:50]
        head(sub_mvals)
        dim(sub_mvals)
        mvals <- sub_mvals
    }

    mvals
}

format_meth <- function(ct_meth){
    # cominbe the cell type methylation together
    head(ct_meth[[1]])

    # get the intersection of CpGs
    cpgs <- Reduce(intersect, (lapply(ct_meth, function(x) rownames(x))))
    head(cpgs)
    length(cpgs)

    # combine brain region data
    ordered_meth <- lapply(ct_meth, function(x) x[cpgs,])
    formatted_meth <- do.call(cbind, ordered_meth)
    head(formatted_meth)
    dim(formatted_meth)
    
    formatted_meth
}

run_pca <- function(formatted_meth, cleaned_mvals){ 
    # estimate dims and get sig pcs 
    head(formatted_meth)
    dim(formatted_meth)

    # remove zero var cpgs - not needed
    var_meth <- formatted_meth[apply(formatted_meth, 1, var) != 0,]
    dim(var_meth)

    # get PCs
    # rows = variables
    # cols = samples
    pca_res <- PCAtools::pca(var_meth)

    # get dim using gavish donoho
    # rows = variables
    # cols = observations
    gv_res <- chooseGavishDonoho(var_meth, var.explained=pca_res$sdev^2, noise=1)
    gv_res

    mp_res <- chooseMarchenkoPastur(var_meth, var.explained=pca_res$sdev^2, noise=1)
    mp_res

    print(paste("gavish donoho dim:", gv_res[[1]], "marchenko pastur dim:", mp_res[[1]]))

    # get the pcs
    names(pca_res)
    pcs <- pca_res$rotated
    head(pcs)

    # get estimated dimension
    est_dim <- gv_res[[1]]
    est_dim

    # subset pcs to estimated dimension 
    sig_pcs <- pcs[,1:est_dim,drop=FALSE] %>%
        as.data.frame()
    head(sig_pcs)
    dim(sig_pcs)

    # get pct var explained
    eig_sq <- pca_res$sdev ^ 2
    pct_var <- eig_sq / sum(eig_sq)
    subset_pct_var <- pct_var[1:est_dim]


    formatted_pcs <- rbind(subset_pct_var, sig_pcs)
    rownames(formatted_pcs)[1] <- "percent_variance"
    head(formatted_pcs)

    mp_var <- sum(pct_var[1:mp_res[[1]]])
    print(paste("Variance explained, GV:", round(sum(subset_pct_var), 2), "MP:", round(sum(mp_var), 2)))


    print("Saving files")

    if (include_schizo){
        savefile <- paste0(savedir, "schizo_int_celltype_pcs.rds")
    } else{
        savefile <- paste0(savedir, "int_celltype_pcs.rds")
        if (cleaned_mvals){
            savefile <- paste0(savedir, "int_cleaned_celltype_pcs.rds")
        }
    }
    savefile
    saveRDS(formatted_pcs, savefile)

    1
}


plot_umap <- function(cleaned_mvals){
    # plot the umap using pcs
    savefile <- paste0(savedir, "celltype_pcs.rds")
    savefile
    if (cleaned_mvals){
        savefile <- paste0(savedir, "int_cleaned_celltype_pcs.rds")
    }
    if (include_schizo){
        savefile <- paste0(savedir, "schizo_int_celltype_pcs.rds")
    }
    savefile
    pcs <- readRDS(savefile)
    pcs <- pcs[-1,]
    head(pcs)
    dim(pcs)



    # load clinical data for formatting the UMAP plot
    clinical <- readRDS("../../output/02_deconvolution/02_svd_batch_correction/restimate_proportions/c2_covariates/INT/remove_age/astro_formatted_clinical.rds")
    head(clinical)

    
    # run umap
    set.seed(1174117174)
    umap_data <- umap(
                     X=sub_pcs,
                      pca=min(100,ncol(pcs)),
                      n_neighbors=500,
                      spread=10,
                      min_dist=5,
                )

    # run umap for schizo data
    # using different parameters
    if (include_schizo){
        dim(sub_pcs)
        set.seed(1174117174)
        umap_data <- umap(
                         X=sub_pcs,
                          n_neighbors=1000,
                          spread=10,
                          min_dist=5,
                    )
    }


    # format the names for plotting
    colnames(umap_data) <- c('UMAP_1', 'UMAP_2')
    head(umap_data)

    # separate celltype from the column names
    formatted_umap <- umap_data %>%
        as.data.frame() %>%
        rownames_to_column('sample_celltype') %>%
        separate(sample_celltype, c('sample', 'celltype'), sep='_', extra='merge')
    head(formatted_umap)

    # save the results
    savefile <- paste0(savedir, "celltype_umap_results.rds")
    saveRDS(formatted_umap, savefile)
    formatted_umap <- readRDS(savefile)

    # format the clinical data to plot
    sub_clinical <- clinical[clinical$individual %in% unique(formatted_umap$sample),] %>%
        mutate(sample = individualID)
    head(sub_clinical)

    #formatted_umap <- formatted_umap %>%
        #left_join(sub_clinical) 
        #filter(!is.na(sex))

    # plot
    savefile <- paste0(savedir, "celltype_umap_plot.png")
    if (include_schizo){
        savefile <- paste0(savedir, "schizo_celltype_umap_plot.png")
    }
    savefile

    # plot the UMAP
    head(formatted_umap)
    cell_type_colors
    updated_ct_colors <- cell_type_colors
    updated_ct_colors['NeuN'] = 'deeppink'
    updated_ct_colors['Olig2'] = 'darkturquoise'
    p <- formatted_umap %>%
        mutate(formatted_ct = fct_recode(celltype,
                                         "Astrocyte"="astro",
                                         "Bulk"="bulk",
                                         "Endothelial"="endo",
                                         "Neuron"="neuron",
                                         "Oligo/OPC"="oligo_opc",
                                         "NeuN"="NeuN",
                                         Olig2="Olig2")) %>%
        mutate(formatted_ct = fct_relevel(formatted_ct,
                                          "Bulk", "Astrocyte", "Endothelial", "Neuron", "Oligo/OPC",
                                          "NeuN", "Olig2")) %>%
        ggplot() +
        geom_point(aes(x=UMAP_1, y=UMAP_2,
                       color=formatted_ct)) +
        scale_color_manual(values=updated_ct_colors, name="Cell type") +
        theme_bw() +
        theme(text=element_text(size=20)) +
        xlab("UMAP 1") +
        ylab("UMAP 2")
    ggsave(p, file=savefile, width=9, height=7)

    savefile <- paste0(savedir, "cleaned_celltype_umap_plot.rds")
    saveRDS(p, savefile)


    1
}

load_schizo <- function(){
    # load sorted data from the schizo study
    # containing nuclei-sorted methylation for neurons and oligos
    chroms <- 1:22

    chrom <- 22 # for testing
    load_chrom <- function(chrom){
        print(paste("Processing chromosome", chrom))
        datafile <- paste0("/home/eulalio/deconvolution/schizo_wgbs/output/03_format_methylation/chr",
                           chrom, "_betas.rds")
        schizo_betas <- readRDS(datafile)
        head(schizo_betas)

        schizo_betas
    }

    schizo_meth <- lapply(chroms, load_chrom)
    schizo_betas <- do.call(rbind, schizo_meth)
    head(schizo_betas)

    schizo_betas
}

main <- function(){
    cell_types <- c('neuron', 'oligo_opc', 'astro', 'endo', 'bulk')

    # load methylation data
    downsample=FALSE
    ct_meth <- lapply(cell_types, load_meth, cleaned_mvals, downsample)

    # combine the methylation together
    print(paste("Formatting methylation"))
    formatted_meth <- format_meth(ct_meth)
    head(formatted_meth)

    # load schizo data for this comparison
    if (include_schizo){
        schizo_betas <- load_schizo()
        head(schizo_betas)
        dim(schizo_betas)

        # intersect cpgs across studies
        common_cpgs <- intersect(rownames(formatted_meth), rownames(schizo_betas))
        length(common_cpgs)

        # subset to these common cpgs
        sub_schizo_betas <- schizo_betas[common_cpgs,] %>%
            na.omit()


        # remove zero variance cpgs
        var_betas <- sub_schizo_betas[apply(sub_schizo_betas, 1, var, na.rm=TRUE) != 0,] %>%
            na.omit()

        # remove values exactly equal to 0 or 1
        betas_matrix <- as.matrix(var_betas)
        betas_matrix[betas_matrix < 0] <- 0
        betas_matrix[betas_matrix > 1] <- 1
        betas_matrix <- (betas_matrix * (ncol(betas_matrix) - 1) + 0.5) / ncol(betas_matrix)
        sum(is.na(betas_matrix))

        # convert to mvals
        schizo_mvals <- beta2m(betas_matrix) %>%
            as.data.frame()
        head(schizo_mvals)

        # combine the schizo with the rosmap data
        head(formatted_meth)
        head(schizo_mvals)

        ordered_cpgs <- intersect(rownames(formatted_meth), rownames(schizo_mvals))
        head(ordered_cpgs)
        length(ordered_cpgs)

        ordered_rosmap <- formatted_meth[ordered_cpgs,]
        ordered_schizo <- schizo_mvals[ordered_cpgs,]

        if (!identical(rownames(ordered_rosmap), rownames(ordered_schizo))) stop("rownames do not match")

        joined_mvals <- cbind(ordered_schizo, ordered_rosmap)
        dim(joined_mvals)
        head(joined_mvals)

        formatted_meth <- joined_mvals
    }

    # run pca
    print(paste("Running PCA"))
    run_pca(formatted_meth, cleaned_mvals)

    plot_umap(cleaned_mvals)

    1
}

main()
