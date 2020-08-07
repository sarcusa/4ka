#!/bin/bash

#SBATCH --time=1400:00
#SBATCH --mem=64000
#SBATCH --mail-type=ALL

module load R
srun Rscript make.R


