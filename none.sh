#!/bin/bash

#SBATCH -J ccimpute
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mamalec@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --mem=128G

export OMP_NUM_THREADS=16

#Load any modules that your program needs
module load intel/19.0.5
module load r/4.1.1 

srun Rscript --vanilla ~/ccImpute/none.R

