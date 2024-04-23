#!/home/eulalio/micromamba/bin/bash

# run the ptwas estimation step

GENEFILE="/home/eulalio/deconvolution/data/ptwas/sample.ptwas_est.in.txt"

PTWAS_est -d $GENEFILE \
    -t 0.5 \
    -n ptwas_test 
