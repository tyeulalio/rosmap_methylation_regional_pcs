# create input files for ENLOC

import os
import gzip
import pandas as pd
from datetime import datetime
import sys
import math
import subprocess
import numpy as np

limit = -1

summary_types = ['avgs', 'pcs', 'cpgs']
cell_types = ['bulk', 'astro', 'endo', 'neuron', 'oligo_opc']
region_types = ['full_gene', 'promoters', 'gene_body', 'preTSS']

# get command line arguments
args = sys.argv

print(args)
summary_type = summary_types[int(args[1])]
cell_type = cell_types[int(args[2])]
region_type = region_types[int(args[3])]
job_idx = int(args[4])
savedir = args[5]
num_jobs = int(args[6])
priors_dir = args[7]
pheno_dir = args[8]
genotype_dir = args[9]

print("savedir:", savedir)


# use this from now on
runname = "{}_{}_{}".format(summary_type, cell_type, region_type)


def process_phenotype(savefile, gene, ordered_samples):
    # print("Processing phenotype data", datetime.now())
    # read in the phenotype data
    # output to savefile in appropriate format
    phenofile = "{}/{}_{}_{}_matched_pheno.bed.gz".format(pheno_dir, cell_type, region_type, summary_type)

    out = open(savefile, "a+")

    header = []
    round = 0
    with gzip.open(phenofile, mode='rb') as f:
        for line in f:
            # convert from binary to utf
            line = line.decode('utf-8')
            # split on tab
            line = line.strip().split('\t')

            # store the header line
            if round == 0:
                header = line

                sub_header = header[4:]
                # print("sub header:", len(sub_header))
                # print("ordered samples:", len(ordered_samples))
                # print(header[:10])
                # if not sub_header == ordered_samples:
                if not np.array_equal(sub_header, ordered_samples):
                    print("SAMPLES OUT OF ORDER")
                    for i in range(0,len(ordered_samples)):
                        print(sub_header[i], ordered_samples[i])
                        if i == 10:
                            break
                    exit()
                round += 1
                continue


            # get the fields for output
            field1 = "pheno"
            pheno_id = line[3]
            group_id = "group"

            if pheno_id != gene:
                continue

            # print("found pheno", pheno_id, gene)
            outlist = [field1, pheno_id, group_id] + line[4:]
            outline = " ".join(outlist)
            # print(outlist[1:10])

            out.write("{}\n".format(outline)) 

            break

            # stop for testing
            # if round == 1:
                # break
            # round += 1

    out.close()

    # use this to match order of samples for geno and controlled pheno
    ordered_samples = header[4:]
    return(ordered_samples)


def get_genes(savedir):
    # get a list of all the genes
    print("Getting gene data", datetime.now())
    phenofile = "{}/{}_{}_{}_matched_pheno.bed.gz".format(pheno_dir, cell_type, region_type, summary_type)
    print("Reading phenotypes from", phenofile)

    genesdf = pd.read_csv(phenofile, sep='\t')

    # get the genes
    geneslist = genesdf.gene_id.values

    # output genes list
    genesfile = "{}/genes_list.txt".format(savedir)
    genesonly = genesdf.gene_id
    genesonly.to_csv(genesfile, index=False, header=False)

    return(geneslist)
    

def process_genotype(sbamfile, priorsfile, ordered_samples, gene, missing_priors_file, genotype_dir):

    # need to adjust gene name for pcs to remove PC
    # print(gene)
    sub_gene = gene

    global summary_type
    if summary_type == 'pcs':
        sub_gene = gene.split('-')[0]

    # print(sub_gene)


    # file with counts for this gene
    geno_file = "{}/{}_genotype_counts.txt".format(genotype_dir, sub_gene)


    if not os.path.exists(geno_file):
        print("no geno file")
        with open(missing_priors_file, 'a+') as f:
            f.write("{}\n".format(sub_gene))
        exit()

    # out = open(sbamfile, "a+")

    cmd = "cat {} | sed 's/{}/{}/g' >> {}".format(geno_file, gene, sub_gene, sbamfile)
    # print("Running cmd", cmd)
    subprocess.run(cmd, shell=True)

    # out.close()

    # write a new snp pip file
    # originally used to remove SNPs but don't need to do that
    # need to rewrite file for some reason or else error in dap-g occurs
    # snpdf.index = snpdf.pheno.values
    # keepdf = snpdf.loc[keepsnps,:]
    # print("Saving SNP file")
    snpfile = "{}/{}/torus/{}.prior".format(priors_dir, runname, gene)

    if not os.path.exists(snpfile):
        print("no priors file")
        with open(missing_priors_file, 'a+') as f:
            f.write("{}\n".format(gene))
        return()

    snpdf = pd.read_csv(snpfile, delim_whitespace=True, names=['pheno', 'pip'])
    snpdf.to_csv(priorsfile, sep=' ', index=False, header=False)



def process_variables(savefile, ordered_samples):
    # print("Processing covariates", datetime.now())
    # read in the variable data
    # output to savefile in appropriate format
    varfile = "{}/{}_{}_{}_matched_covariates.tsv".format(pheno_dir, cell_type, region_type, summary_type)

    out = open(savefile, "a+")

    # need to transpose data so read in with pandas
    covs = pd.read_csv(varfile, sep='\t')

    # samples alraedy match ordering --
    
    # need to save each column to output file
    cols = covs.columns
    for i in range(1, covs.shape[1]):
        # grab the column
        col = covs.iloc[:, i].values
        colname = cols[i]

        # create the output
        outlist = ['controlled', colname, 'group'] + list(col)

        outline = " ".join([str(x) for x in outlist])

        out.write("{}\n".format(outline))

        # if i == limit:
            # break


    out.close()


def load_genotypes():
    # load the genotypes into a dictionary 
    # for quicker look-ups

    # read in the genotype data
    # output to savefile in appropriate format
    genofile = "/home/eulalio/deconvolution/new_rosmap/output/05_qtl_analysis/06_merge_bims/all_chroms_ref_counts_hg38.traw"

    header = []
    samples_idx = []

    genotypes_dict = {}

    round = 0
    with open(genofile, 'r') as f:
        for line in f:
            # split on tab
            line = line.strip().split('\t')

            # store the header line
            if round == 0:
                header = line

                # remove 0_ from the sample names
                header = [x.replace('0_', '') for x in header]

                # get the index of the ordered samples in header
                # samples_idx = [header.index(x) for x in ordered_samples]
                round += 1
                continue



            # get the fields for output
            field1 = "geno"
            pheno_id = line[1]
            group_id = "group"

            # replace NA values with 0 - DAP can handle missing values
            # line = ["2" if i=="NA" else i for i in line]

            # switch counts from ref to alt by subtracting from 2
            genos = line[6:]
            genos = [str(2-int(x)) if x != 'NA' else x for x in genos]


            # create output
            # outlist = [field1, pheno_id, group_id] + [line[x] for x in samples_idx]
            outlist = [field1, pheno_id, group_id] + genos
            outline = " ".join(outlist)

            genotypes_dict[pheno_id] = outline

            if round % 1000000 == 0:
                print("Loaded {} genotypes {}".format(round, datetime.now()))

            # stop for testing
            # if round == 1000000:
                # break
            round += 1

    return(genotypes_dict)

def get_ordered_samples(genotype_dir):
    datafile = "{}/ordered_samples.txt".format(genotype_dir)

    ordered_samples = pd.read_csv(datafile, header=None)

    ordered_samples = ordered_samples.iloc[:,0].values
    # print(ordered_samples)

    return(ordered_samples)


def main():
    # create the output directories
    global savedir
    if not os.path.exists(savedir):
        print("creating output directory", savedir)
        os.makedirs(savedir)
    savedir = "{}/{}".format(savedir, runname)
    if not os.path.exists(savedir):
        print("creating output directory", savedir)
        os.makedirs(savedir)

    priorsdir = "{}/priors".format(savedir)
    if not os.path.exists(priorsdir):
        print("creating output directory", priorsdir)
        os.makedirs(priorsdir)

    sbamdir = "{}/sbams".format(savedir)
    if not os.path.exists(sbamdir):
        print("creating output directory", sbamdir)
        os.makedirs(sbamdir)


    # create this file to keep track of gene with missing priors
    missing_priors_file = "{}/missing_priors_file.txt".format(savedir)
    with open(missing_priors_file, 'w+') as f:
        None

    # need to output each gene separately - get list of genes
    geneslist = get_genes(savedir)

    # subset the genes list for this job
    print("Starting job", job_idx)
    total_genes = len(geneslist)
    genes_per_job = math.ceil(total_genes  / num_jobs)
    start_idx = job_idx * genes_per_job
    end_idx = start_idx + genes_per_job
    # make sure we don't go past the end of array
    end_idx = min(end_idx, total_genes)

    print("total genes", total_genes)
    print("genes per job:", genes_per_job)
    print("start:", start_idx, "end", end_idx)

    if (end_idx < start_idx):
        print("No genes to process")
        return()

    geneslist = geneslist[start_idx:end_idx]
    print(len(geneslist))
    print(geneslist[0])
    print(geneslist[len(geneslist)-1])

    # print("Loading genotype {}".format(datetime.now()))
    # genotypes_dict = load_genotypes()


    # write an output file for each gene
    numgenes = 0
    for gene in geneslist:
        # gene = geneslist[gene_idx]
        print(gene)

        # print("Loading genotype {}".format(datetime.now()))
        # genotypes_dict = load_genotype(gene)


        # create the savefile for output
        sbamfile = "{}/{}_fastqtl_singl_snp_output.sbam".format(sbamdir, gene)
        # print("output file:", sbamfile)
        # if os.path.isfile(sbamfile):
            # continue

        # clear the file after testing
        with open(sbamfile, "w+"):
            None

        # create priors output file
        priorsfile = "{}/{}_priors.txt".format(priorsdir, gene)
        with open(priorsfile, "w+"):
            None

        # print("SBAM output file:", sbamfile)
        # print("PRIORS output file:", priorsfile)


        # process pheno
        # print("Processing pheno")
        global genotype_dir
        ordered_samples = get_ordered_samples(genotype_dir)

        _ = process_phenotype(sbamfile, gene, ordered_samples)

        # process geno
        # print("Processing geno")
        # process_genotype2(sbamfile, priorsfile, ordered_samples, gene, genotypes_dict, missing_priors_file)
        # print("Writing genotype {}".format(datetime.now()))
        process_genotype(sbamfile, priorsfile, ordered_samples, gene, missing_priors_file, genotype_dir)

        # process controlled variables
        # print("Processing controlled variables")
        process_variables(sbamfile, ordered_samples)


        if numgenes % 1000 == 0:
            print("Processed {} genes out of {} {}".format(numgenes, len(geneslist), datetime.now()))

        numgenes += 1




if __name__ == '__main__':
    main()
