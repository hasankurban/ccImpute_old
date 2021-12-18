#!/bin/bash

#SBATCH -J ccimpute
#SBATCH -p gpu-debug
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mamalec@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:00:00
#SBATCH --mem=128G

export OMP_NUM_THREADS=16

#Load any modules that your program needs
module unload gcc
module load intel/19.0.5
module load gcc/9.3.0
module load boost/gnu/1.72.0
srun ~/R-4.1.2/build/bin/Rscript --vanilla ~/ccImpute/benchmark2.R

