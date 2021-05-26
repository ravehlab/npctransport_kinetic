#!/usr/local/bin/bash

# Get folders correctly from USER_CONFIG.sh:
#MYDIR=`dirname $(readlink -f $0)`
source /cs/labs/ravehb/ravehb/Projects/npctransport_kinetic/slurm/USER_CONFIG.sh
TRAVERSE_CMD=${SLURM_FOLDER}/sbatch_traverse.sh

for MW in 27 34 41 47 54 67; do 
   sbatch ${TRAVERSE_CMD} ${MW} $1 $2 $3;
done
