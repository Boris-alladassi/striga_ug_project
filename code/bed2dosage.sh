#!/bin/bash
#SBATCH --account=aces
#SBATCH --partition=secondary
#SBATCH --job-name=bed_to_dosage
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --output=bed_to_dosage.out
#SBATCH --error=bed_to_dosage.err
#SBATCH --mail-user=aboris@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda
module load anaconda3
source ~/.bashrc
conda activate popgen_env


bed_file="/projects/illinois/aces/shared/aboris/anal_stri_ga/ld_pruned/data_pruned_final"
output_file="/projects/illinois/aces/shared/aboris/code/dosage_dt2/dosage_output"

plink --bfile $bed_file --recode A --allow-extra-chr --out $output_file

conda deactivate
