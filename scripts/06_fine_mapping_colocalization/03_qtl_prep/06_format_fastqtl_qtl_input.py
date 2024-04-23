# create the gwas input for fastqtl

import os
from datetime import datetime
import gzip

summary_types = ['avgs', 'pcs', 'cpgs']
cell_types = ['bulk']
region_types = ['full_gene']


# get command line arguments
args = sys.argv
summary_type = summary_types[args[1]]
cell_type = cell_types[args[2]]
region_type = region_types[args[3]]

# use this from now on
runname = "{}_{}_{}".format(summary_type, cell_type, region_type)

# create output directory if it doesn't exit
savedir = "/path/to/output/06_fine_mapping_colocalization/03_qtl_prep/06_fastqtl_qtl_input"
if not os.path.exists(savedir):
    os.makedirs(savedir)
savedir = "{}/{}".format(savedir, runname)
if not os.path.exists(savedir):
    os.makedirs(savedir)

def format_output(out):
    # format the output for gwas fastqtl

    # read in the qtl file
    qtl_file = "/path/to/output/06_fine_mapping_colocalization/03_qtl_prep/05_fastqtl_qtl_annotations/{}/fastenloc.qtl.annotation.vcf.gz".format(runname)

    i = 0
    with gzip.open(qtl_file, 'rb') as f:
        for line in f:
            line = line.decode("utf-8").strip().split()

            # get the gene from the line
            gene = line[5].split(":")[0]
            # gene = gene.replace("-", ".")

            # add gene before snp 
            new_line = line
            new_line[2] = "{}_{}".format(gene, line[2])

            # save line to the output file
            out.write("{}\n".format("\t".join(new_line)))

            # if i == 3:
                # break
            i+=1


def main():

    # output file
    savefile = "{}/fastqtl_qtl_input.vcf".format(savedir)
    out = open(savefile, 'w+')

    format_output(out)

    out.close()

    # gzip the file
    zipfile = "{}.gz".format(savefile)
    with open(savefile, 'rb') as f_in, gzip.open(zipfile, 'wb') as f_out:
        f_out.writelines(f_in)


if __name__ == '__main__':
    main()
