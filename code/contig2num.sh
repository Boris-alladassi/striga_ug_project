#!/bin/bash
#SBATCH --account=aces
#SBATCH --partition=secondary
#SBATCH --job-name=contig_to_num
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --mem=64000
#SBATCH --output=contig_to_num.out
#SBATCH --error=contig_to_num.err
#SBATCH --mail-user=aboris@illinois.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load conda
module load anaconda3
source ~/.bashrc
conda activate popgen_env

# Input PLINK prefix (without extension)
PREFIX="/projects/illinois/aces/shared/aboris/anal_stri_ga/ld_pruned/data_pruned_final"

# Output prefix
OUTPREFIX="${PREFIX}_numeric"

echo ">>> Recoding $PREFIX.bim to integer chromosome codes..."

# Step 1. Extract unique contigs/scaffolds
cut -f1 ${PREFIX}.bim | sort | uniq > contigs.txt

# Step 2. Make mapping: integer -> contig
nl -w1 -s' ' contigs.txt > contig_map.txt

# Step 3. Replace contigs in .bim with integers
awk 'NR==FNR {map[$2]=$1; next} { $1=map[$1]; print }' contig_map.txt ${PREFIX}.bim > ${OUTPREFIX}.bim

# Step 4. Copy .bed and .fam to match new prefix
cp ${PREFIX}.bed ${OUTPREFIX}.bed
cp ${PREFIX}.fam ${OUTPREFIX}.fam

echo ">>> Done!"

conda deactivate
