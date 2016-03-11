#!/bin/ba sh
#MSUB -r testseamount    # Request na me
#MSUB -n 4               # Number of tasks to use
#MSUB -T 1800            # Elapsed time limit in seconds
#MSUB -o nh_%I.txt       # Standard output.%I is the job id
#MSUB -e nh_%I.txt       # Error output.%I is the job id
#MSUB -q standard
#MSUB -A gen9178         # IMPORTANT to specifuy here the correct project ID

# Project ID

set -x

module unload intel

module load intel/16.0.2.181 bullxmpi netcdf/4.2

cd ${BRIDGE_MSUB_PWD}

ccc_mprun -n 4 testseamount

