library(GenomicRanges)
library(tidyverse)

# check that the annotations for gwas and qtl SNPs match


## -- Wightman et al GWAS -- ##
# load gwas data
datafile <- paste0("/path/to/data/ad_gwas/wightman_etal_2021_results.txt.gz")
gwas <- read_tsv(datafile)
head(gwas)

# save the gwas data in a format for liftover
formatted_gwas <- gwas %>%
    mutate(Chromosome=paste0('chr', chr)) %>%
    select('#CHROM'=Chromosome, POS=PosGRCh37, REF=testedAllele,
           ALT=otherAllele
    ) %>%
    mutate(QUAL='.', FILTER='PASS', INFO='.') %>%
    filter(!REF %in% c('-'))
head(formatted_gwas)

formatted_gwas %>%
    filter(grepl('[^a-zA-Z]', ALT))

# change the ID name
formatted_gwas <- formatted_gwas %>%
    mutate(ID=paste(`#CHROM`, POS, REF, ALT, sep='_')) %>%
    select(`#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO)
head(formatted_gwas)

savefile <- paste0("/path/to/data/ad_gwas/wightman_etal_2021_results.vcf")

write_lines("##fileformat=VCFv4.2", savefile)
formatted_gwas |> colnames() |> paste0(collapse = "\t") |> write_lines(savefile, append = TRUE)
write_tsv(formatted_gwas, savefile, append=TRUE)
1
