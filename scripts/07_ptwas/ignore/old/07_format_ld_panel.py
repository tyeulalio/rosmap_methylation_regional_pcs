# format the LD panel using ROSMAP genotype data
import os
from datetime import datetime
import sys
import shutil

# get command line arg
args = sys.argv
if len(args) < 2:
    print("Enter chromosome as argument")
    exit()
CHROM = sys.argv[1]

TIMEFORMAT="%H:%M:%S"

def process_chromosomes(savedir):
    # read in lifted rosmap chromosome file
    rosmap_file = "/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf"

    # open the output file
    outfile = "{}/rosmap_wgs_normalized_chrom{}.vcf".format(savedir, CHROM)
    out = open(outfile, "w+")
    print("Writing output to", outfile)

    # curr_chrom = None
    print("Processing chromosome", CHROM)
    i = 0
    header = None
    with open(rosmap_file, 'r') as f:
        for line in f:
            # write comment lines directly to output
            if line.startswith("##"):
                out.write(line)
                continue

            # get the header line
            if line.startswith("#"):
                header = line.rstrip().split()
                out.write(line)
                continue

            i = i+1
            if i % 100000 == 0:
                print("Processing line", i, datetime.now().strftime(TIMEFORMAT))
            line = line.rstrip().split()

            # print(header[0:10])
            # print(line[0:10])

            linedict = dict(zip(header, line))
            
            # adjust chrom to match example
            linedict["#CHROM"] = linedict["#CHROM"].replace("chr", "")
            # if (not curr_chrom == linedict["#CHROM"]) | (not curr_chrom):
                # curr_chrom = linedict["#CHROM"]
                # print("Processing chromsome", curr_chrom, datetime.now().strftime("%H:%M:%S"))

            # skip line if we're not on this chromosome
            if int(linedict["#CHROM"]) != int(CHROM):
                continue
            # exit once we pass the chromosome 
            if int(linedict["#CHROM"]) > int(CHROM):
                break

            # need to adjust subject info
            subjects = header[9:] 
            # print(subjects[0:10])

            for subj in subjects:
                subj_info = linedict[subj]
                # print(subj_info)

                # take only the first part of info
                subj_geno = subj_info.split(':')[0]
                # print(subj_geno)

                formatted_geno = subj_geno.replace("/", "|")
                # print(formatted_geno)

                # replace info in linedict
                linedict[subj] = formatted_geno


            # print(linedict)
            outlist = [linedict[x] for x in header]
            outline = "\t".join(outlist)
            out.write("{}\n".format(outline))

    out.close()

def split_chromosomes(chromosomes):
    # split up the chromosomes from the main normalized wgs vcf

    # read in lifted rosmap chromosome file
    rosmap_file = "/home/eulalio/deconvolution/new_rosmap/data/rosmap_wgs_liftover_hg38/rosmap_wgs_all_chroms_normalized.vcf"

    tmpdir = "/tmp"

    # split header 
    print("Splitting header", datetime.now().strftime(TIMEFORMAT))
    header_file = "{}/header.txt".format(tmpdir)
    cmd = "grep '#' {} > {}".format(rosmap_file, header_file) 
    os.system(cmd)

    # split off chromosomes
    for chrom in chromosomes:
        print("Splitting chromosome" chrom)
        chrom_file = "{}/wgs_chrom{}.txt".format(tmpdir, chrom)
        cmd = "grep '^chr{}' {} > {}".format(rosmap_file, chrom_file)



def main():
    # create the output directory
    savedir = "/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/06_formatted_ld_panel"
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    chromosomes = range(20, 22)
    print("Processing chromosomes", chromosomes)

    # split up files
    print("Splitting VCF file")
    split_chromosomes(chromosomes)
    

    # process_chromosomes(savedir)


if __name__ == '__main__':
    main()
