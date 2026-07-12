#!/bin/bash
#SBATCH --account=aces
#SBATCH --partition=secondary
#SBATCH --job-name=vcf2plink_filter
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --output=vcf2plink_filtering.out
#SBATCH --error=vcf2plink_filtering.err
#SBATCH --mail-user=aboris@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda
module load anaconda3
source ~/.bashrc
conda activate popgen_env

# Input and output paths
VCF_INPUT="/projects/illinois/aces/shared/aboris/anal_stri_ga/striga_filtered_84k.vcf"
VCF_PREFIX="/projects/illinois/aces/shared/aboris/anal_stri_ga/geno_data_filtered"

## vcftools already filtered 20% NAs and .05 MAF, leading to 84808 SNPs out of a possible 35339684 Sites
# Now, further filtering and Convert to PLINK format for R analysis
plink --vcf $VCF_INPUT --mind 0.40 --geno 0.10 --maf 0.10 --allow-extra-chr --make-bed --out $VCF_PREFIX

conda deactivate
