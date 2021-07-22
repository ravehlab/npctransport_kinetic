#!/bin/env bash
#SBATCH --job-name=multiprocess
#SBATCH --output=logs/multiprocess_%j.out
#SBATCH --time=12:00:00
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
    output="run03",
    passive_range=(0.01,0.15),
    npc_traverse_range=(1,1000),
    k_on_range=(0.001,5),
    nx=33,
    ny=33,
    n_passive=22,
    cargo_concentration_M=0.1e-6,
    Ran_concentration_M=20e-6,
    v_N_L=2194e-15,
    v_C_L=627e-15,
    pickle_file="run03_reverse_N_C_volume.pkl",
    number_of_processors=${SLURM_CPUS_PER_TASK})
EOF
