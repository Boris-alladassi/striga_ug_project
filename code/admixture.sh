#!/bin/bash
#SBATCH --account=aces
#SBATCH --partition=secondary
#SBATCH --job-name=popgen_admixture
#SBATCH --time=03:30:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --output=admixture.out
#SBATCH --error=admixture.err
#SBATCH --mail-user=aboris@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda
module load anaconda3
source ~/.bashrc
conda activate popgen_env

bed_file="/projects/illinois/aces/shared/aboris/anal_stri_ga/ld_pruned/data_pruned_final_numeric.bed"
# ===== STEP 3: Run ADMIXTURE =====
# Choose K values (e.g., 1 to 10)
for K in {1..10}; do admixture --cv $bed_file $K | tee log${K}.out; done

# Collect cross-validation errors
grep -h CV log*.out > cv_error.txt

conda deactivate
