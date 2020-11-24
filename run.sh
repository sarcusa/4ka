#!/bin/bash

#SBATCH --time=400:00 ##was 3000, then 1000
#SBATCH --mem=64000
#SBATCH --mail-type=ALL

##module load openmpi
module load R
##mpirun Rscript make.R
srun Rscript make.R

