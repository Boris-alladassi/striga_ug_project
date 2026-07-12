#!/bin/bash
#SBATCH --account=aces
#SBATCH --partition=secondary
#SBATCH --job-name=vcf_filtering
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --mem=256000
#SBATCH --output=vcf_filtering.out
#SBATCH --error=vcf_filtering.err
#SBATCH --mail-user=aboris@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda
module load anaconda3
source ~/.bashrc
conda activate popgen_env

# Input and output paths
VCF_INPUT="/projects/illinois/aces/shared/aboris/anal_stri_ga/filtered_Raw_SNPs_biallelic_freqfilter.vcf.gz"
VCF_PREFIX="/projects/illinois/aces/shared/aboris/anal_stri_ga/geno_data_filtered"

# Step 1: Remove SNPs with >10% missing data and MAF < 0.10
#vcftools --gzvcf $VCF_INPUT --max-missing 0.90 --mac 1 --maf 0.10 --recode --recode-INFO-all --out step1
# Step 3: Convert to PLINK format for R analysis
plink --vcf step1.recode.vcf --allow-extra-chr --make-bed --out $VCF_PREFIX

conda deactivate
