#!/bin/bash

#SBATCH --time=400:00 ##was 3000, then 1000, then 400
#SBATCH --mem=64000
#SBATCH --mail-type=ALL
##SBATCH --partition=el6

##module load openmpi
module load R
##mpirun Rscript make.R
srun Rscript make.R

