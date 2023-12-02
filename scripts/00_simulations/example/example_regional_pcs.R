library(minfiData)
library(PCAtools)
library(tidyverse)

# download example methylation data
# obtain rPCS for a sample region


# Load the MethylSet data
data(MsetEx.sub)

# Display the first few rows of the dataset for a preliminary view
head(MsetEx.sub)

# Extract methylation M-values from the MethylSet
# M-values are logit-transformed beta-values and are often used in differential
# methylation analysis for improved statistical performance.
mvals <- getM(MsetEx.sub)
head(mvals)
dim(mvals)

# subset to a example region with N cpgs
N <- 50
region_mvals <- mvals[1:N,]
head(region_mvals)
dim(region_mvals)

# save the region methylation
savefile <- paste0("./example_region_methylation_df.rds")
saveRDS(region_mvals, savefile)

# apply whatever centering/scaling to the data
scaled_mvals <- scale(t(region_mvals), center=TRUE, scale=FALSE) %>%
    t()
head(scaled_mvals)
dim(scaled_mvals)

# run PCA - expects features as rows, samples as columns
# function seems to center the columns of data frame by default
pca_res <- PCAtools::pca(scaled_mvals)

names(pca_res)

# get the rotated data (regionalPCs)
rpcs <- pca_res$rotated
head(rpcs)

savefile <- paste0("./example_rpcs_df.rds")
saveRDS(rpcs, savefile)
