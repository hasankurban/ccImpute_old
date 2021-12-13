#!/bin/bash

#SBATCH -J ccimpute
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mamalec@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:05:00
#SBATCH --mem=64G

#Load any modules that your program needs
module load intel/19.0.5
module load r/4.1.1 


srun R CMD BATCH ccImpute2.R
