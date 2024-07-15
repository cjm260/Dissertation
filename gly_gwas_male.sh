#!/bin/bash
#SBATCH --job-name=gwas_Gly_Male
#SBATCH --output=/home/cjm260/gwas_files/gwas_results/Gly/gwas_GlyMale_%A_%a.out
#SBATCH --error=/home/cjm260/gwas_files/gwas_results/Gly/gwas_GlyMale_%A_%a.err
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# Load the module environment
source /etc/profile.d/modules.sh

# Load necessary modules (adjust as needed)
module load plink/2.00-alpha

# Add plink2 directory to PATH
export PATH=$PATH:/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex

# Define file paths and other parameters
PHENO_FILE="/home/cjm260/gwas_files/glyPheno.txt"
PHENO_NAME="invRankNorm_resids"
OUT_DIR="/home/cjm260/gwas_files/gwas_results/Gly"
plink2="plink2"

# Array of chromosomes to process
CHROMOSOMES=(1 2 3 5 7 8 9 12 16 20)
CHR=${CHROMOSOMES[$SLURM_ARRAY_TASK_ID-1]}

# Define the input file prefix and output file prefix based on the chromosome
PFILE_PREFIX="/home/cjm260/gwas_files/filtINTERVALdata/gly/chr${CHR}_filtered"
OUT_PREFIX="${OUT_DIR}/gwas_male_chr${CHR}"

# Create output directory if it doesn't exist
#mkdir -p ${OUT_DIR}

# Run Plink for the current chromosome
${plink2} --pfile ${PFILE_PREFIX} \
          --pheno ${PHENO_FILE} --pheno-name ${PHENO_NAME} \
          --keep /home/cjm260/gwas_files/male_samples.txt \
          --mind 0.03 --geno 0.05 --maf 0.01 --hwe 0.00001 \
          --linear \
          --out ${OUT_PREFIX}
