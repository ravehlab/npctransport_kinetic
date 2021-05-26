#!/usr/local/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000M

# Get folders correctly from USER_CONFIG.sh:
MYDIR=`dirname $(readlink -f $0)`
source $MYDIR/USER_CONFIG.sh
#SCRIPTS_FOLDER="/cs/labs/ravehb/ravehb/Projects/npctransport_kinetic/"

echo Running simulation with MW=$1, no force NPC traverse rate = $2, force NPC traverse rate = $3

python3 $SCRIPTS_FOLDER/MW_stats_by_force.py $1 100 $2 $3
