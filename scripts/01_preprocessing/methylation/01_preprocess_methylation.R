# this is code to 
# - load all data from .idat files
# - run quality control on samples + plot density and snp-heatmaps 
# - remove low quality probes and normalize data for each brain region separately
# all functions are below
# last updated:
# # 0/2020, Anna-Lena Lang
# 12/14/2021, Tiffany Eulalio

# wateRmelon annotation files - make sure these are installed, don't need to load yet
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
# library(IlluminaHumanMethylationEPICmanifest) 


#doParallel must be manually loaded to enable multithreading on servers 
library(doParallel) 
registerDoParallel(cores = 16) #This must be called or multithreading will not work

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#library(impute)
library(minfi)
#library(limma)
#library(RColorBrewer)
#library(missMethyl)
#library(Gviz)
#library(DMRcate)
#library(stringr)
#library(grid)
library(wateRmelon)
#library(methylumi)
library(ChAMP)
#library(MethylAid)
#library(DNAmArray)
#library(gemPlot)
#library(gplots)
library(ewastools)
#library(reshape)
#library(data.table)
library(maxprobes)
#library(randomcoloR)
#library(plotly)
#library(ggplot2)
#library(dplyr)
#library(viridis)
#library(circlize)
#library(gnn)
library(RnBeads.hg19)
library(RnBeads) # BiocManager::install('RnBeads.hg38') # make sure this is installed for RnBeads
#library(lumi)
library(data.table)
library(tidyverse)


# define important directories

# set to direct to your main folder
baseDir <- paste("/oak/stanford/groups/smontgom/eulalio/deconvolution/new_rosmap/")
list.files(baseDir)

datadir <- paste0(baseDir, "data/ROSMAP_data/Epigenetics/Epigenetics (DNA methylation array)/") # main data folder
list.files(datadir)

idatdir <- paste0(datadir, "IDAT Files") # folder containing idat files
list.files(idatdir)

savedir <- paste0(baseDir, "output/01_preprocessing/methylation/01_preprocessing_methylation/") # folder for output
dir.create(savedir)




# ---- functions -----
load_samples <- function(){
    # load the data files and combine some phenotypes together
    # save the combined df and return the df

    # NOTE: you might want to check all data files to see
    # if there are any other important features that we might
    # want to carry forward in the analyses


    # adjust the sample annotations provided
    # need a Barcode column
    sample_annotations <- paste0(datadir, "SYNAPSE_METADATA_MANIFEST.tsv")
    samples <- read_tsv(sample_annotations) %>%
        as.data.frame()
    head(samples)

    table(samples$analysisType)
    head(samples$path)

    # check the analysis type
    # filter based on this
    filtered_samples <- samples %>%
        filter(dataSubtype == 'raw',
               fileFormat == 'idat'
        )
    head(filtered_samples)

    table(filtered_samples$fileFormat)

    dim(samples)
    dim(filtered_samples)
    head(filtered_samples)

    # get the base filename
    formatted_annotations <- filtered_samples %>%
        mutate(Basename = str_remove(name, "_[a-zA-Z]{3}\\.idat")) %>%
        select(-path, -name) %>%
        unique()
    head(formatted_annotations)
    length(unique(formatted_annotations$Basename))
    length(unique(formatted_annotations$specimenID))
    formatted_annotations %>%
        select(specimenID, Basename) %>%
        count(specimenID) %>%
        arrange(-n) %>%
        head()

    dim(formatted_annotations)

    formatted_annotations
}



# -- PART 1 

# qc for samples filters samples that: with bisulfite conversion efficiency <80%, fail probe control metrics using ewastools,  are outliers in snp-heatmap (compared within each individual)
qualitycontrol_samples <- function(targets, rgset){

  # bisulfite conversion efficiency
  print(paste0("sample number : ", nrow(targets)))
  print("removing samples with bisulfite conversion efficiency < 80%")
  # get the bisulfite conversion effeciency
  # remove samples with < 80
  bsc <- bscon(rgset) %>% # checks bisulfite conversion efficiency
            as.data.frame() %>%
            rename(bsc=".")
  head(bsc)

  # check for failed samples
  bad_bsc <- bsc %>%
      filter(bsc < 80)
  head(bad_bsc)

  # get the samples to keep
  keep_bsc <- setdiff(colnames(rgset), rownames(bad_bsc))
  head(keep_bsc)
  length(keep_bsc)
  
  # subset the rgset
  rgset_bsc <- rgset[,keep_bsc]
  rgset_bsc
  rgset

  targets <- targets[targets$sample_id %in% colnames(rgset_bsc), ]
  print(paste0("sample number : ", ncol(rgset_bsc)))
  
  
  # use contol_metrics function from ewastools to remove failed samples
  print("use contol_metrics function from ewastools to remove failed samples")
  head(targets)
  
  # need to get filepath for each samples' idats
  # read in idats for the control_metrics function
  meth <- read_idats(targets$Basename, quiet = FALSE) 
  names(meth)
  head(meth$meta)

  # make sure targets and meth are in same order
  identical(meth$meta$sample_id, targets$sample_id)
  
  # get control metrics defined by illumina using ewastools
  ctrls <- control_metrics(meth)
  head(ctrls)

  targets$failed <- sample_failure(ctrls) # check for failed samples
  table(targets$failed) # 85 samples failed

  # filter to those that passed filtering
  targets_qcfilt <- targets %>%
      filter(!failed)
  head(targets_qcfilt)
  dim(targets_qcfilt)
  dim(targets)

  rgset_qcfilt <- rgset_bsc[, targets_qcfilt$sample_id]
  rgset_qcfilt # 84 samples removed

  print(paste0("sample number: ", nrow(targets_qcfilt)))

  
  # save files
  print("sample filtering done. saving new datasets:")
  
  savefile <- paste0(savedir, "rgset_qcsfilt_hg19.rds")
  saveRDS(rgset_qcfilt, savefile)

  rgset <- readRDS(savefile)
  
  savefile <- paste0(savedir, "targets_qcsfilt_hg19.rds")
  saveRDS(targets_qcfilt, savefile)

  targets <- readRDS(savefile)

  
  res <- list(rgset_qcfilt=rgset,
       targets_qcfilt=targets)
  res
}


rgset <- rgset_raw # set this for manual testing
preprocess_pipeline <- function(targets, rgset){
  # first running quality control on samples 
  # (checking bisuflite conversion, snp heatmaps and 
  # ewastools quality control to remove any outlier samples, 
  # this is done on full dataset
  res <- qualitycontrol_samples(targets, rgset) 
  names(res)
  
  # set these based on our newly filtered results
  targets <- res$targets_qcfilt
  rgset <- res$rgset_qcfilt
  

  # check quality of probes in each sample and remove probes with low quality 
  # (hight detection p-value detP > 0.01, beadcount, 
  # remove probes on sex chromosomes and remove cross-reactive probes. 
  # Finally remove do dye-bias-correction using noob and normalization 
  # using bmiq and remove all technical and biological replicates)
  qc_probes_norm_remove_reps(targets, rgset)
}


# -- PART 2
# bmiq normalization
bmiq_norm <- function(targets_champ, betas_champ){
  # BMIQ normalization 
  print(paste0("starting bmiq normalization"))
  betas_bmiq <<- champ.norm(beta = betas_champ, 
                            method = "BMIQ", 
                            plotBMIQ = FALSE, 
                            arraytype = "450K", cores = 1)
  
  savefile <- paste0(savedir, "betas_bmiq.rds")
  saveRDS(betas_bmiq, savefile)

  savefile <- paste0(savedir, "targets_bmiq.rds")
  saveRDS(targets_champ, savefile)

  
  print("bmiq done and results saved")
}

# remove biological and technical replicates 
remove_reps <- function(targets_ewasfilt_bmiq_incl_R_850k, beta_ewasfilt_bmiq_incl_R_850k){
  # remove technical and biological replicates
  print(paste0("removing tR and bR"))
  print(paste0("number of bR: " , sum(targets_ewasfilt_bmiq_incl_R_850k$biological_replicate == 1)))
  print(paste0("number of tR: " , sum(targets_ewasfilt_bmiq_incl_R_850k$technical_replicate == 1)))
  
  targets_ewasfilt_bmiq_orig_850k <<- targets[targets$original == 1,]
  beta_ewasfilt_bmiq_orig_850k <<- beta[ ,colnames(beta)%in% targets_ewasfilt_bmiq_orig_850k$sample_id]
  
  saveRDS(beta_ewasfilt_bmiq_orig_850k, file =  paste0(savedir, "beta_ewasfilt_bmiq_orig_850k.rds"))
  saveRDS(targets_ewasfilt_bmiq_orig_850k, file = paste0(savedir, "targets_ewasfilt_bmiq_orig_850k.rds"))
  print(paste0("tR and bR removed, new datasets saved. finished analysis for brain region, total count of samples: ", 
               nrow(targets_ewasfilt_bmiq_orig_850k)))
}


filter_snps <- function(betas, targets){
    # use the genotype VCFs from ROSMAP to get AF

    head(betas)
    head(targets)

    # need to map probes to positions
    head(ann450k)
    table(ann450k$Enhancer)

    # get just the relevant columns
    pos450k <- ann450k %>%
        as.data.frame() %>%
        select(chr, pos, strand, Name) %>%
        mutate(strand_char = ifelse(strand == '+', 'p', 'n')) %>%
        mutate(cpg = paste(chr, pos, sep='_'))
    head(pos450k)

    savefile <- paste0(savedir, "probe_position_map_hg37.rds")
    savefile
    saveRDS(pos450k, savefile)

    
    # get the probes that we have for betas
    ordered_pos <- pos450k[rownames(betas),]
    head(ordered_pos)

    # make sure the order matches
    identical(rownames(ordered_pos), rownames(betas))


    # read in the annoated vcfs to get AF
    chrom = 22

    process_chrom <- function(chrom){
        print(paste("Processing chrom", chrom, Sys.time()))
        # read in the 
        datafile <- paste0("../../../data/rosmap_wgs_harmonization_variant_calling/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_", chrom,
                           ".recalibrated_variants.annotated.vcf.gz")
        annots <- fread(datafile, skip='#CHROM')
        head(annots)

        # keep only the info with AF
        formatted_annots <- annots %>%
            mutate(AF = str_match(INFO, "AF=[0-9\\.]+")) %>%
            mutate(AF = str_remove(AF, "AF=")) %>%
            select(-INFO) %>%
            mutate(AF = as.numeric(AF))
        head(formatted_annots)

        summary(formatted_annots$AF)

        # subset to SNPS with AF > 0.01
        thresh=0
        remove_snps <- formatted_annots %>%
            filter(AF > thresh)
        head(remove_snps)

        print(paste("Number of SNPs on chromosome", chrom, "with AF >", thresh, "=", nrow(remove_snps)))
            
        remove_snps
    }

    remove_snps <- lapply(1:22, process_chrom)
    remove_snps_df <- do.call(rbind, remove_snps)

    head(remove_snps_df)
    
    savefile <- paste0(savedir, "remove_snps_af_", thresh, ".rds")
    savefile
    saveRDS(remove_snps_df, savefile)

    remove_snps_df <- readRDS(savefile) 


    # number of snps that we find in each chromosomes
    table(remove_snps_df['#CHROM'])


    # check how many of these snps overlap our probe positions
    head(ordered_pos)
    head(remove_snps_df)

    formatted_snps <- remove_snps_df %>%
        rename(chr='#CHROM',
               pos=POS
        ) %>%
        mutate(chr = paste0('chr', chr)) 
    head(formatted_snps)

    table(formatted_snps$chr)
        
    joined_snps_pos <- ordered_pos %>%
        left_join(formatted_snps)
    head(joined_snps_pos)

    # these are the SNPs to remove
    # do filtering here
    thresh = 0
    matched_snps <- joined_snps_pos %>%
        filter(!is.na(AF)) %>%
        filter(AF > thresh)

    head(matched_snps)
    nrow(matched_snps)
    length(unique(matched_snps$Name))


    # get the filtered probes with straight overlap
    head(ordered_pos)
    filtered_probes <- ordered_pos %>%
        filter(!cpg %in% matched_snps$cpg)
    head(filtered_probes)
    nrow(filtered_probes)


    table(matched_snps$chr)



    # filter betas to match
    head(betas)
    head(matched_snps)

    # decide which filtering method to use here
    betas_filtered <- betas %>%
        as.data.frame() %>%
        rownames_to_column('Name') %>%
        filter(!Name %in% matched_snps$Name)
    head(betas_filtered)

    # replace probe name wit position
    head(ordered_pos)
    formatted_betas <- betas_filtered %>%
        left_join(ordered_pos[c('Name', 'cpg')]) %>%
        column_to_rownames('cpg') %>%
        select(-Name)
    head(formatted_betas)

    dim(formatted_betas)

    savefile <- paste0(savedir, "snp", thresh, "_window0_filtered_betas.rds")
    savefile
    saveRDS(formatted_betas, savefile)


    savefile <- paste0(savedir, "snp", thresh, "_window0_filtered_targets.rds")
    savefile
    saveRDS(targets, savefile)

}


# qc for probes, normalize using noob and bmiq and finally remove biological/technical replicates 
qc_probes_norm_remove_reps <- function(targets, rgset) {
  print(paste0("starting quality control of probes"))
  print(paste0("probe number: ", nrow(rgset)))

  # use noob for dye bias Correction
  print(paste0("running noob dye bias correction"))
  mset_noob_qcsfilt_incl_R_850k <- preprocessNoob(rgset, dyeCorr = TRUE, dyeMethod = "single")
  RSet_noob_qcsfilt_incl_R_850k <- ratioConvert(mset_noob_qcsfilt_incl_R_850k, what = "both", keepCN = TRUE)


  # map the genome so we can filter out sex chromosomes
  GRset_noob_qcsfilt_incl_R_850k <- mapToGenome(RSet_noob_qcsfilt_incl_R_850k)
  print(paste0("probe number: ", nrow(GRset_noob_qcsfilt_incl_R_850k)))

  # save this early filtere data
  GRset_noob_qcsfilt_incl_R_850k
  betas <- getBeta(GRset_noob_qcsfilt_incl_R_850k)
  savefile <- paste0(savedir, "noob_preprocessed_meth.rds")
  saveRDS(betas, savefile)


  getLocations(GRset_noob_qcsfilt_incl_R_850k)

  # 1) data includes males and females --> remove probes on the sex chromosomes
  print("removing probes on sex-chromosomes")
  featureNames(GRset_noob_qcsfilt_incl_R_850k)
  unique(ann450k$chr)
  head(ann450k)

  # select probes that aren't on the sex chromosomes
  keep_probes <- ann450k %>%
      as.data.frame() %>%
      select(Name, chr) %>%
      filter(!chr %in% c('chrX', 'chrY')) 
  head(keep_probes)

  grset_sex_filtered <- GRset_noob_qcsfilt_incl_R_850k[keep_probes$Name,]
  grset_sex_filtered

  nrow(GRset_noob_qcsfilt_incl_R_850k) - nrow(grset_sex_filtered)


  # 2) remove probes with SNPs at CpG site
  thresh = 0
  print(paste("removing probes associated with SNPs of minor allele frequency (MAF) >= ", thresh))
  grset_snp_filtered <- dropLociWithSnps(grset_sex_filtered, maf=thresh)
  print(paste0("probe number: ", nrow(grset_snp_filtered)))

  nrow(grset_sex_filtered) - nrow(grset_snp_filtered)

  savefile <- paste0(savedir, "snp_filtered_thresh", thresh, ".rds")
  savefile
  saveRDS(grset_snp_filtered, savefile)

  grset_snp_filtered <- readRDS(savefile)
  
  # 3) finally exclude cross reactive probes (based on maxprobes github)
  print("removing cross reactive probes based on maxprobes::xreactive_probes")
  xreact <- xreactive_probes(array_type = "450K")
  head(xreact)
  keep2 <- !(featureNames(grset_snp_filtered) %in% xreact)
  grset_xreact_filtered <- grset_snp_filtered[keep2,]

  nrow(grset_snp_filtered) - nrow(grset_xreact_filtered)

  # convert to beta values
  betas <- getBeta(grset_xreact_filtered)
  print(paste0("probe number: ", nrow(betas)))
  
  
  # 4) filter samples by detection P value from ewastools, using champ.filter 
  # subset all datasets to match
  # need targets, detp, beta, and beadcount to match

  # get detection p -values 
  detP_ewas <- detectionP.minfi(rgset) # using detP value from ewastools - return detection p vals matrix
  detP_ewas <- detP_ewas[, order(colnames(detP_ewas))]
  
  savefile <- paste0(savedir, "detP_ewas_raw.rds")
  print(paste0("saving dataset with detection pvals", savefile))
  saveRDS(detP_ewas, file=savefile)

    detP_ewas <- readRDS(savefile)
    dim(detP_ewas)
    nrow(rgset) - nrow(detP_ewas)

  # get overlap of samples 
  ordered_samples <- intersect(targets$sample_id, colnames(detP_ewas)) %>%
      intersect(colnames(betas)) 
  head(ordered_samples)
  length(ordered_samples)

  # order the targets 
  targets <- targets[match(ordered_samples, targets$sample_id),] %>%
      mutate(Sample_Name = sample_id)
  identical(ordered_samples, targets$sample_id)

  # get overlap of probes
  ordered_probes <- intersect(rownames(detP_ewas), rownames(betas))
  length(ordered_probes)

  # order detp and betas
  ordered_detp <- detP_ewas[ordered_probes, ordered_samples]
  dim(ordered_detp)

  ordered_betas <- betas[ordered_probes, ordered_samples]
  dim(ordered_detp)


  # get beadcount for filtering later on
  beadcount <- beadcount(rgset) %>%
      as.data.frame()
  head(beadcount)

  # fix sample names for beadcount
  formatted_beadcount <- beadcount
  colnames(formatted_beadcount) <- gsub("^X", "", colnames(formatted_beadcount))
  head(formatted_beadcount)

  ordered_beadcount <- formatted_beadcount[ordered_probes, ordered_samples]
  dim(ordered_beadcount)


  print("removing: probes with detP > 0.01, probes with <3 beads in at least 5% of samples, samples where > 10% of probes have detection p value > 0.01,
        probes with NA using champ.filter")

  
  # filter based on beadcount and detP from ewastools using champ.filter
    identical(rownames(ordered_betas), rownames(ordered_detp))
    identical(colnames(ordered_betas), colnames(ordered_detp))

    
  champ <- champ.filter(beta = ordered_betas, 
                        pd = targets,
                        detP = ordered_detp,
                        beadcount = ordered_beadcount,
                        autoimpute = FALSE,
                        ProbeCutoff = 0,
                        SampleCutoff = 0.1,
                        detPcut = 0.01,
                        filterDetP = TRUE,
                        beadCutoff = 0.05,
                        filterBeads = TRUE,
                        filterNoCG = FALSE,
                        filterSNPs = FALSE,
                        population = NULL,
                        filterMultiHit = FALSE,
                        filterXY = FALSE,
                        fixOutlier = FALSE,
                        arraytype = "450K")
  betas_champ <- champ$beta
  targets_champ <- champ$pd
  print(paste0("probe number: ", nrow(betas_champ)))
  
  # save results
  savefile <- paste0(savedir, "betas_filtered.rds")
  saveRDS(betas_champ, savefile)

  betas_champ <- readRDS(savefile)
  
  savefile <- paste0(savedir, "targets_filtered.rds")
  saveRDS(targets_champ, savefile)

  targets_champ <- readRDS(savefile)


  print(paste0("Finished probe filtering, datasets incl replicates saved"))

    # bmiq normalization  
  bmiq_norm(targets_champ, betas_champ)
  
  # this removes technical/biological replicates - we don't have these so okay to skip
   #remove_reps(targets_ewasfilt_bmiq_incl_R_850k, beta_ewasfilt_bmiq_incl_R_850k)

  savefile <- paste0(savedir, "betas_bmiq.rds")
  betas <- readRDS(savefile)
  
  savefile <- paste0(savedir, "targets_bmiq.rds")
  targets <- readRDS(savefile)

  betas <- betas_var

  # filter probes that overlap SNPs with AF >= 0.01
  filter_snps(betas, targets)
} 

# get filtered rgset
get_filtered_rgset <- function(){
    # load the original rgset
  savefile <- paste0(savedir, "rgset_qcsfilt_hg19.rds")
  rgset <- readRDS(savefile)
  rgset
  
  thresh = 0
    savefile <- paste0(savedir, "snp", thresh, "_window0_filtered_betas.rds")
    savefile

    formatted_betas <- readRDS(savefile)
    head(formatted_betas)

    # filter the rgset down to the final set of cpgs
    dim(rgset)
    dim(formatted_betas)


    # map the cpg IDs to cpg positions to match betas to rgset
    head(rgset)
    head(formatted_betas)

  mapped_rgset <- mapToGenome(rgset)
  head(mapped_rgset)

  # check that row order is the same in mapped and rgset
  i = 1
  identical(mapped_rgset[i,], rgset[i,])

  str(mapped_rgset)
  head(mapped_rgset$Meth)

  str(rgset)
  head(rgset$Meth)


    savefile <- paste0(savedir, "probe_position_map_hg37.rds")
    savefile
    probe_map <- readRDS(savefile)
    head(probe_map)

    filtered_probes <- probe_map %>%
        filter(cpg %in% rownames(formatted_betas))
    head(filtered_probes)

    dim(filtered_probes)
    dim(probe_map)
    dim(formatted_betas)


    # filter the rgset to the final set
    head(mapped_rgset)
    filtered_grgset <- mapped_rgset[filtered_probes$Name,]
    filtered_grgset

    dim(filtered_rgset)
    dim(formatted_betas)

    filtered_rgset <- as(filtered_grgset, "RGChannelSet")

    props <- estimateCellCounts(rgset, 
                                 compositeCellType='DLPFC', 
                                 referencePlatform="IlluminaHumanMethylation450k",
                                 cellTypes=c('NeuN_neg', 'NeuN_pos'),
                                 processMethod='preprocessNoob'
                                )
    head(props)

    savefile <- paste0(savedir, "neun_estimated_cell_counts.rds")
    savefile
    saveRDS(props, savefile)

}




# ---- run analysis ----

# get annotation data for epic array - needed for mapping probes to genome
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

# read in sample sheet containing phenotype and array information
targets <- load_samples() 
head(targets)
dim(targets)

# read in raw idats - take a long time
rgset_raw <- read.metharray.exp(targets = targets, extended = TRUE)


# assign column names
colnames(rgset_raw)
head(rgset_raw)
targets$sample_id <- colnames(rgset_raw)


# run full analysis pipeline 
# this includes quality control of samples, probes, normalization and removal of replicates
preprocess_pipeline(targets_raw, rgset_raw)
