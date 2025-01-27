#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=8GB
#SBATCH --time=04:00:00
#SBATCH --account=open
#SBATCH --partition=basic
#SBATCH --job-name=FTDE
export OMPI_MCA_mpi_cuda_support=0

echo "Job started on $(hostname) at $(date)"

source ~/.bashrc

module use /storage/icds/RISE/sw8/modules
module load vasp/vasp-6.3.1vtst

conda activate findtde

echo "Start: `date`"

$FTDE_DIR/findtde/find_tde -c midpoint -p vasp

echo "End: `date`"
