library(ggplot2)
library(tidyverse)

# use this for testing
chromdir = ""
#chromdir = "chrom18_only/"

# check and format eigenstrat output for qtl analysis
savedir <- paste0("../../output/05_qtl_analysis/08_format_eigenstrat/", chromdir)
dir.create(savedir, showWarnings=FALSE)


# -- look at variance explained for PCs
datafile <- paste0("../../output/05_qtl_analysis/07_geno_pca/", chromdir, "eigenstrat_smartpca_hg38.pca.evec")
vals <- read_table(datafile, n_max=1, col_names=FALSE) %>%
    t()
vals <- vals[-1,]
head(vals)

vals_df <- data.frame(eval=vals) %>%
    rownames_to_column('rn') %>%
    select(-rn) %>%
    na.omit() %>%
    mutate(eval=as.numeric(eval))
head(vals_df)
dim(vals_df)



# compute total variance
total_var <- sum(vals_df$eval)
total_var

vals_pct <- vals_df %>%
    mutate(pct_var = eval / total_var,
           cum_var = cumsum(pct_var)
    )
head(vals_pct)
vals_pct

savefile <- paste0(savedir, "pca_pct_var_hg38.tsv")
savefile
write_tsv(vals_pct, savefile)


# -- load the PCs
datafile <- paste0("../../output/05_qtl_analysis/07_geno_pca/", chromdir, "eigenstrat_smartpca_hg38.pca.evec")
pcs <- read_table(datafile, skip=1, c('individualID', paste0("PC", 1:10), 'race'))
head(pcs)
dim(pcs)

formatted_pcs <- pcs %>%
    separate(individualID, c('fid', 'specimenID'), sep=':') %>%
    select(-race, -fid) 
head(formatted_pcs)

savefile <- paste0(savedir, "genotype_pcs_hg38.rds")
savefile
saveRDS(formatted_pcs, savefile)

# identify samples that were removed
datafile <- paste0("../../output/05_qtl_analysis/05_filter_mishap/keep_samples.txt")
keep_samples <- read_tsv(datafile, col_names=c('fam', 'sampleID'))
head(keep_samples)
dim(keep_samples)

length(intersect(keep_samples$sampleID, formatted_pcs$specimenID))

removed_samps <- data.frame(outliers=setdiff(keep_samples$sampleID, formatted_pcs$specimenID))
head(removed_samps)
dim(removed_samps)

savefile <- paste0(savedir, "population_outliers_hg38.txt")
savefile
write_tsv(removed_samps, savefile, col_names=FALSE)


# remove samples that were filtered out in methylation data
#datafile <- paste0("/home/eulalio/deconvolution/new_rosmap/output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
#matched_clinical <- readRDS(datafile)
#head(matched_clinical)
#dim(matched_clinical)

#head(keep_samples)

#keep_matched_samples <- keep_samples %>%
    #filter(sampleID%in% matched_clinical$wgs.specimenID)
#head(keep_matched_samples)
#dim(keep_matched_samples)
#dim(keep_samples)



datafile <- paste0("../../output/01_preprocessing/01_subset_samples/matched_clincal_meta.rds")
meta <- readRDS(datafile)
head(meta)

head(formatted_pcs)
keep_meta <- meta %>%
    filter(wgs.specimenID %in% formatted_pcs$specimenID)
head(keep_meta)
dim(keep_meta)


# correlate PCs with phenotypes of interest
#datafile <- paste0("../../output/episcore/02_svd_batch_correction/astro_formatted_clinical.rds")
#clinical <- readRDS(datafile)

#head(clinical)
#head(formatted_pcs)

#length(intersect(keep_meta$individualID, clinical$individualID))

#matched_clinical <- clinical %>%
    #filter(individualID %in% keep_meta$individualID) 
#head(matched_clinical)

head(formatted_pcs)
ordered_pcs <- formatted_pcs %>%
    mutate(wgs.specimenID=specimenID) %>%
    left_join(keep_meta[c('wgs.specimenID', 'individualID')]) %>%
    select(-specimenID, -wgs.specimenID) %>%
    column_to_rownames('individualID')
head(ordered_pcs)

ordered_pcs <- ordered_pcs[meta$individualID,]

#identical(matched_clinical$individualID, rownames(ordered_pcs))

savefile <- paste0(savedir, "matched_meta.rds")
savefile
saveRDS(meta, savefile)


# traits that we'll test against 
names(meta)
traits <- c('individualID', 'Study', 'educ', 'age_at_visit_max', 'age_first_ad_dx', 
            'age_death', 'cts_mmse30_first_ad_dx', 'cts_mmse30_lv', 'pmi', 'braaksc', 'ceradsc',
            'cogdx', 'dcfdx_lv', 'diagnosis', 'sex', 'apoe4',
            'meth.batch'
            #'neuron', 'oligo_opc', 'astro', 'endo'
)

# subset the targets
target_traits <- meta[traits]
head(target_traits)



# put svs in dataframe
#svs_df <- as.data.frame(svs)

head(ordered_pcs)
pcs_df <- ordered_pcs
head(pcs_df)

identical(rownames(pcs_df), meta$individualID)

n.sv <- ncol(pcs_df)

n.traits <- length(traits) - 1

pvals_mat <- matrix(nrow=n.sv, ncol=n.traits)
cols <- rep(0,n.traits)

target_traits <- as.data.frame(target_traits)

# correlte the SVs with each of the traits
sv_num <- 1
trait_num <- 1
for (sv_num in 1:n.sv){
    for (trait_num in 1:n.traits){
        trait <- traits[trait_num+1]
        print(trait)

        if (class(target_traits[,trait]) != 'numeric'){
            pvals_mat[sv_num, trait_num] <- kruskal.test(pcs_df[,sv_num] ~ target_traits[,trait])$p.value
        } else{
            pvals_mat[sv_num, trait_num] <- cor.test(pcs_df[,sv_num], target_traits[,trait], method='spearman',
                                                     exact=FALSE)$p.value
        }

        cols[trait_num] <- trait
    }
}
colnames(pvals_mat) <- cols
pvals_mat

# store a corrected and uncorrected version of pvals
pvals_long <- pvals_mat %>%
    as.data.frame() %>%
    mutate(sv=paste("sv", row_number(), sep='_')) %>%
    gather(key='trait', value='pval', -sv) %>%
    group_by(trait) %>%
    mutate(bf_pval = p.adjust(pval, method='bonferroni'))
head(pvals_long)


# save the output
savefile <- paste0(savedir, "genotype_pc_trait_lm_pvals_hg38.rds")
savefile
saveRDS(pvals_long, savefile)

pvals_long <- readRDS(savefile)

head(pvals_long)
pvals_long %>%
    filter(bf_pval < 0.1)

function plot_pvals(){
    savefile <- paste0(savedir, "genotype_pc_trait_pvals_heatmap.png")
    savefile

    p <- ggplot(pvals_long, aes(x=sv, y=trait)) +
        geom_tile(aes(fill=bf_pval)) +
        scale_fill_gradient(limits=c(0,0.1), low='red', high='white') +
        geom_text(aes(label=round(bf_pval,2))) +
        theme_bw() +
        theme(text=element_text(size=18)) +
        ggtitle(paste("Genomic PCs, trait correlation p-values"))
    ggsave(plot=p, filename=savefile, width=10, height=7)
}
