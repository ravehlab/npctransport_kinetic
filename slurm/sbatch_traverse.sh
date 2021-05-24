#!/bin/csh

#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2000M

python3 /cs/labs/ravehb/littlem/npctransport_kinetic/MW_stats_by_force.py $1 40 $2 $3 $4
