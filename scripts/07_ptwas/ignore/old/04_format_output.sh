#!/home/eulalio/micromamba/bin/bash

# remove chr from the chromosome in QTL gambit output

GAMBITDIR="/home/eulalio/deconvolution/new_rosmap/output/07_ptwas/03_gambit_db"

GAMBITFILE="${GAMBITDIR}/all_gene.ptwas_weights.gambit.vcf.gz"
OUTFILE="${GAMBITDIR}/formatted_all_gene.ptwas_weights.gambit.vcf"

# get the header from gabmit file
echo "printing header"
zcat $GAMBITFILE | head | grep "^#" > $OUTFILE

# remove chr from chromosome and from ID
echo "printing remaining formatted lines"
zcat $GAMBITFILE | grep -v "^#" | sed -e "s/chr//g" >> $OUTFILE

# bgzip the output file
echo "zipping output file"
bgzip $OUTFILE

# tabix the output
echo "running tabix"
tabix -p 'vcf' "${OUTFILE}.gz"

