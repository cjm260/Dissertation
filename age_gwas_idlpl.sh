#!/bin/bash
#SBATCH --job-name=gwas_idlpl_age
#SBATCH --error=/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_idlpl_%A_%a.err
#SBATCH --array=1-12
#SBATCH --mem=8G
#SBATCH --time=02:00:00

# Load the module environment
source /etc/profile.d/modules.sh

# Load necessary modules (adjust as needed)
module load plink/2.00-alpha

# Add plink2 directory to PATH
export PATH=$PATH:/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex

# Array of chromosomes to process
CHROMOSOMES=(1 2 5 7 8 9 11 15 17 19 20)
CHR=${CHROMOSOMES[$SLURM_ARRAY_TASK_ID-1]}

# Define the input file prefix based on the chromosome
PFILE_PREFIX="/home/cjm260/gwas_files/filtINTERVALdata/idlpl/chr${CHR}_filtered"

# Run GWAS 1
PHENO_FILE1="/home/cjm260/gwas_files/ageStratMetabPheno/idlplPheno.txt"
PHENO_NAME1="invRankNorm_resids"
OUT_DIR1="/home/cjm260/gwas_files/gwas_results/idlpl1"
OUT_PREFIX1="${OUT_DIR1}/gwas_chr${CHR}"

mkdir -p ${OUT_DIR1}

plink2 --pfile ${PFILE_PREFIX} \
       --pheno ${PHENO_FILE1} --pheno-name ${PHENO_NAME1} \
       --keep /home/cjm260/gwas_files/data_bin_18_29.txt \
       --mind 0.03 --geno 0.05 --maf 0.01 --hwe 0.00001 \
       --linear \
       --out ${OUT_PREFIX1}

# Run GWAS 2
PHENO_FILE2="/home/cjm260/gwas_files/ageStratMetabPheno/idlplPheno.txt"
PHENO_NAME2="invRankNorm_resids"
OUT_DIR2="/home/cjm260/gwas_files/gwas_results/idlpl2"
OUT_PREFIX2="${OUT_DIR2}/gwas_chr${CHR}"

mkdir -p ${OUT_DIR2}

plink2 --pfile ${PFILE_PREFIX} \
       --pheno ${PHENO_FILE2} --pheno-name ${PHENO_NAME2} \
       --keep /home/cjm260/gwas_files/data_bin_30_44.txt \
       --mind 0.03 --geno 0.05 --maf 0.01 --hwe 0.00001 \
       --linear \
       --out ${OUT_PREFIX2}

# Run GWAS 3
PHENO_FILE3="/home/cjm260/gwas_files/ageStratMetabPheno/idlplPheno.txt"
PHENO_NAME3="invRankNorm_resids"
OUT_DIR3="/home/cjm260/gwas_files/gwas_results/idlpl3"
OUT_PREFIX3="${OUT_DIR3}/gwas_chr${CHR}"

mkdir -p ${OUT_DIR3}

plink2 --pfile ${PFILE_PREFIX} \
       --pheno ${PHENO_FILE3} --pheno-name ${PHENO_NAME3} \
       --keep /home/cjm260/gwas_files/data_bin_45_59.txt \
       --mind 0.03 --geno 0.05 --maf 0.01 --hwe 0.00001 \
       --linear \
       --out ${OUT_PREFIX3}

# Run GWAS 4
PHENO_FILE4="/home/cjm260/gwas_files/ageStratMetabPheno/idlplPheno.txt"
PHENO_NAME4="invRankNorm_resids"
OUT_DIR4="/home/cjm260/gwas_files/gwas_results/idlpl4"
OUT_PREFIX4="${OUT_DIR4}/gwas_chr${CHR}"

mkdir -p ${OUT_DIR4}

plink2 --pfile ${PFILE_PREFIX} \
       --pheno ${PHENO_FILE4} --pheno-name ${PHENO_NAME4} \
       --keep /home/cjm260/gwas_files/data_bin_60_76.txt \
       --mind 0.03 --geno 0.05 --maf 0.01 --hwe 0.00001 \
       --linear \
       --out ${OUT_PREFIX4}
