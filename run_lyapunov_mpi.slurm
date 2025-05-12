#!/bin/bash
#SBATCH --job-name=lyap-hybrid
#SBATCH --account=beozcan
#SBATCH --partition=orfoz
#SBATCH -C weka
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=110
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --mail-user=berk@metu.edu.tr
#SBATCH --mail-type=ALL

echo "Job started on $(hostname) at $(date)"

module purge
module load comp/gcc/12.3.0
module load lib/openmpi/4.1.6
module load lib/gsl/2.7.0

export OMPI_MCA_pml=ucx
export OMPI_MCA_btl="^openib"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

cd /arf/scratch/beozcan/Lyapunov-Analysis || exit 1

mpirun -np ${SLURM_NTASKS} \
       --map-by ppr:1:node:PE=${SLURM_CPUS_PER_TASK} \
       --bind-to core \
       -x OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK} \
       ./builds/main

echo "Job finished at $(date)"
