#!/usr/local/bin/bash

#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2000M

# Get folders correctly from USER_CONFIG.sh:
#MYDIR=`dirname $(readlink -f $0)`
source /cs/labs/ravehb/ravehb/Projects/npctransport_kinetic/slurm/USER_CONFIG.sh

echo Running simulation with MW=$1, no force NPC traverse rate = $2, force NPC traverse rate = $3 RAN factor = $4

python3 $SCRIPTS_FOLDER/MW_stats_by_force.py $1 100 $2 $3 $4
