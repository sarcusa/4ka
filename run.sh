#!/bin/bash

#SBATCH --time=3000:00
#SBATCH --mem=64000
#SBATCH --mail-type=ALL

module load openmpi
module load R
mpirun Rscript make.R


