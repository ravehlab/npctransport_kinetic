#!/bin/env bash
#SBATCH --job-name=multiprocess
#SBATCH --output=logs/multiprocess_%j.out
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=122

# #SBATCH --nodes=32

source /cs/staff/ravehb/raveh_lab/External/py3/bin/activate
if [ -z ${SLURM_CPUS_PER_TASK+x} ] ; then
    SLURM_CPUS_PER_TASK=1
fi
echo CPUs per task ${SLURM_CPUS_PER_TASK}
python - << EOF
import subprocess
import stats_grid
import os

stats_grid.get_stats_on_grid(
    output="tmp",
    passive_range=(0.01,0.15),
    npc_traverse_range=(10,1000),
    k_on_range=(0.01,5),
    nx=11,
    ny=11,
    n_passive=10,
    cargo_concentration_M=0.1e-6,
    pickle_file="tmp.pkl",
    number_of_processors=${SLURM_CPUS_PER_TASK})
EOF
