#!/bin/bash

#SBATCH -J ccimpute
#SBATCH -p gpu-debug
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mamalec@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-node=15
#SBATCH --time=00:05:00
#SBATCH --mem=64G

#Load any modules that your program needs
module load cudatoolkit/11.2


#Run your program
srun ./benchmark.out 1 2048 80 1 2 10
#srun ./sample2.out

