
import sys
import subprocess
import pandas as pd
import time

#command to create score file:
#echo "/home/cjm260/gwas_files/gwas_results/idlpl1/OPGS003434_model_filt.txt" > /rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/modelfiles/score_file_idlpl.txt

if __name__ == "__main__":
    #PGS_file = str(sys.argv[1])
    #PGS_id = PGS_file.split('.')[0]

    prs_cal_tool = "/home/cjm260/rds/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/GRS_resources/calc_PS_lvls.sh"
    genotype_prefix = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr"

    
    #"/rds/project/asb38/rds-asb38-ceu-ukbiobank/genetics/P7439/post_qc_data/imputed/HRC_UK10K/plink_format/GRCh37/pgen/ukb_imp_v3_dedup_chr"


    #work_dir ="/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/yx322/UKB_omics_PGS/Somalogic_full"
    #out_path = "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/yx322/UKB_omics_PGS/Somalogic_full"
    work_dir = "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/workingdirectory/idlpl"
    out_path = "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/workingdirectory/idlpl"

    score_file = "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/modelfiles/score_file_idlpl.txt"


    
    #"/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/PGS_UKB/selected_scores/selectedScores_v1-2_Somalogic_trans_5e8_cis_5e8.txt"


     #"36000"
    subprocess.call([prs_cal_tool, '--time','12:0:0', 
                     "--partition", "icelake-himem", 
                     "--score-file", score_file,
                     "--type",'l', 
                     "--score-rsid", 'rsid', "--score-chr", "chr","--score-pos",'pos',"--score-EA", "effect_allele", "--score-OA", 'other_allele', '--score-weight', 'effect',
                     '--account', 'INOUYE-SL3-CPU',
                     "--cohort-name",'UKB',
                     "--single-out", 'UKB_Nightingale_idlpl',
                     "--genotype-prefix", genotype_prefix,
                     "--work", work_dir,
                     "--mem-per-chr", "100000", "--out",out_path])