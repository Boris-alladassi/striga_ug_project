#!/bin/bash
#SBATCH --account=aces
#SBATCH --partition=secondary
#SBATCH --job-name=ld_pruning
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --output=ld_pruning.out
#SBATCH --error=ld_pruning.err
#SBATCH --mail-user=aboris@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda
module load anaconda3
source ~/.bashrc
conda activate popgen_env


bed_file="/projects/illinois/aces/shared/aboris/anal_stri_ga/geno_data_filtered"
out_file="/projects/illinois/aces/shared/aboris/anal_stri_ga/geno_data_filtered2"
#output_file="/projects/illinois/aces/shared/aboris/anal_stri_ga/ld_pruned/data_pruned"
#output_file2="/projects/illinois/aces/shared/aboris/anal_stri_ga/ld_pruned/data_pruned_final"

plink --bfile $bed_file --allow-extra-chr --set-missing-var-ids @:# --make-bed --out $out_file
plink --bfile $out_file --allow-extra-chr --indep-pairwise 50 10 0.2 --out data_pruned
plink --bfile $out_file --allow-extra-chr --extract data_pruned.prune.in --make-bed --out data_pruned_final

conda deactivate
