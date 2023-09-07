#!/usr/bin/perl

# This script was adapted from the Eigensoft SmartPCA tool
# It is tailed for the ROSMAP analysis

use strict;
use warnings;

# Check if Eigensoft SmartPCA is in PATH
unless (`which smartpca.perl`){
    die "Error smartpca.perl not found. Please add it to your PATH.\n"
}

# MUST put smartpca bin directory in path for smartpca.perl to work
# You can add SmartPCA to your path by using these lines that follow
# and pointing to the correct location on your machine
#$ENV{'PATH'} = "/home/eulalio/programs/EIG-7.2.1/bin:$ENV{'PATH'}"; 
# add the executables to the path as well
#$ENV{'PATH'} = "/home/eulalio/programs/EIG-7.2.1/src/eigensrc:$ENV{'PATH'}"; 


# Directories and file paths
my $EXAMPLEDIR="/home/eulalio/programs/EIG-7.2.1/EIGENSTRAT/";
my $OUT="/home/eulalio/deconvolution/new_rosmap/output/05_qtl_analysis/07_geno_pca/eigenstrat_smartpca_hg38";

# Initialize the smarpca.perl command
$command = "smartpca.perl";


# input files
# geno file - PACKEDPED format
$command .= " -i ../../output/05_qtl_analysis/06_merge_bims/all_chroms_hg38.bed ";
# snp file - PACKEDPED format
$command .= " -a ../../output/05_qtl_analysis/06_merge_bims/all_chroms_ref_counts_hg38.map ";
# indiv file
$command .= " -b ../../output/05_qtl_analysis/06_merge_bims/all_chroms_ref_counts_hg38.ped " ;

# output files
$command .= " -o ${OUT}.pca ";
$command .= " -p ${OUT}.plot ";
$command .= " -e ${OUT}.eigenvalues ";
$command .= " -l ${OUT}.log ";

# Parameters
# number of PCs to output - default = 10
$command .= " -k 10 ";
# max iterations - default = 5
$command .= " -m 5 ";
# number of PCs to remove outliers during each iteration - default = 10
$command .= " -t 10 ";
# sigma, number of SDs to be removed as outlier - default = 6.0
$command .= " -s 6.0 ";

# Execute the command
print("$command\n");
system("$command");
