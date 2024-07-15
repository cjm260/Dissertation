library(dplyr)
library(tidyr)
library(data.table)
library(survival)
#Coronary atherosclerosis
#################
#CORONARY ATHEROSCLEROSIS
#Example Sex Stratified Cox Reg for idlplcine
#Covariate file
ukb_qc_file <- "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/others/sampleQC_fromUKB_withHeaders.txt"
df_ukb_qc <- fread(ukb_qc_file, skip = 31, data.table = FALSE)
df_ukb_qc <- df_ukb_qc[, c("#UKB_ID1", "genotyping.array", paste0("PC", 1:10))]
mapData <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/reference_files/genetic/reference_files/full_release/sampleID_map.txt")
df_ukb_qc <- merge(df_ukb_qc, mapData, by.x = '#UKB_ID1', by.y = 'UKB_sample_ID')

#Coronary atherosclerosis
df_pheno <- fread("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/PGSCatalog/PheWAS_202310/20240119_phenotyped/Phecode_411.4.csv.gz")
df_pheno <- df_pheno %>%
  select(eid, sex, age_at_recruitment, PHECODE_AgeAsTimescale, PHECODE_AgeAsTimescale_Years) %>%
  rename(idno = eid, age = age_at_recruitment, PHENOTYPE = PHECODE_AgeAsTimescale, CENSOR_AGE = PHECODE_AgeAsTimescale_Years)
# Incorporate PCs and array info
df_pheno_all <- merge(df_pheno, df_ukb_qc, by.x = 'idno', by.y = 'Adiposity_sample_ID')

#Full Model PGS File
pgs_file <- "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/omics_PGS_scores/UKB_Nightingale_5e8/UKB_Nightingale.sscore.gz" 
pgs_full <- fread(pgs_file, sep = '\t')
# Rename columns for the PGS file of all the UKB samples
colnames(pgs_full) <- c('IID', sapply(colnames(pgs_full)[-1], function(col) gsub("\\.pct", "_", strsplit(col, "_")[[1]][2])))

# Only analyzing idlpl
pgs_full_idlpl <- pgs_full %>%
  select(IID, idlpl)

# Reduced model PGS File for idlpl
pgs_reduced <- fread("/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/workingdirectory/idlpl/UKB_Nightingale_idlpl.sscore.gz")
pgs_reduced <- pgs_reduced %>%
  rename(idlpl = score_sum)

# Data frame to save results
df_save <- data.frame(Regression = character(), Trait = character(), HR = numeric(), HR_low = numeric(), HR_high = numeric(), pvalue = numeric(), zscore = numeric(), cindex = numeric(), N_Case = integer(), N_All = integer(), stringsAsFactors = FALSE)

run_cox <- function(df, model_name) {
  cox_model <- coxph(Surv(CENSOR_AGE, PHENOTYPE) ~ adjusted_score + sex, data = df)
  summary_cox <- summary(cox_model)
  c_index <- summary_cox$concordance[1]
  hr <- summary_cox$coefficients[1, 2]
  hr_low <- summary_cox$conf.int[,"lower .95"]
  hr_high <- summary_cox$conf.int[,"upper .95"]
  pvalue <- summary_cox$coefficients[1,5]
  zscore <- summary_cox$coefficients[1,4]
  n_cases <- sum(df$PHENOTYPE == 1)
  n_all <- nrow(df)
  
  df_result <- data.frame(Regression = model_name, Trait = "idlpl", HR = hr, HR_low = hr_low, HR_high = hr_high, pvalue = pvalue, zscore = zscore, cindex = c_index, N_Case = n_cases, N_All = n_all)
  return(df_result)
}

# Function to prepare data, standardize, and adjust PGS for PCs
# Function to prepare data, standardize, and adjust PGS for PCs
prepare_data <- function(df_pgs, df_pheno_all) {
  # Check for missing values in key columns before merge
  if (any(is.na(df_pgs$IID)) | any(is.na(df_pheno_all$idno))) {
    stop("Missing values found in key columns 'IID' or 'idno'. Please remove or impute missing values.")
  }
  
  # Ensure key columns are of the same type
  df_pgs$IID <- as.character(df_pgs$IID)
  df_pheno_all$idno <- as.character(df_pheno_all$idno)
  
  # Merge the data frames
  df_pheno_test <- merge(df_pgs, df_pheno_all, by.x = 'IID', by.y = 'idno')
  score <- colnames(df_pgs)[2]
  
  df_pheno_test <- df_pheno_test %>%
    select(IID, all_of(score), sex, PHENOTYPE, CENSOR_AGE, age, genotyping.array, paste0('PC', 1:10))
  
  #remove rows with NAs
  df_pheno_test = na.omit(df_pheno_test)
  
  print(df_pheno_test)
  # Standardize continuous variables
  df_pheno_test[,2] = scale(df_pheno_test[,2])
  df_pheno_test$PC1 = scale(df_pheno_test$PC1)
  df_pheno_test$PC2 = scale(df_pheno_test$PC2)
  df_pheno_test$PC3 = scale(df_pheno_test$PC3)
  df_pheno_test$PC4 = scale(df_pheno_test$PC4)
  df_pheno_test$PC5 = scale(df_pheno_test$PC5)
  df_pheno_test$PC6 = scale(df_pheno_test$PC6)
  df_pheno_test$PC7 = scale(df_pheno_test$PC7)
  df_pheno_test$PC8 = scale(df_pheno_test$PC8)
  df_pheno_test$PC9 = scale(df_pheno_test$PC9)
  df_pheno_test$PC10 = scale(df_pheno_test$PC10)
  
  # Adjust PGS for PCs and genotyping array
  reg_trait <- as.formula(paste0(score, ' ~ ', paste0('PC', 1:10, collapse = ' + '), ' + as.factor(genotyping.array)'))
  adj_result <- lm(reg_trait, data = df_pheno_test)
  df_pheno_test$adjusted_score <- residuals(adj_result)
  
  
  return(df_pheno_test)
}

# Prepare data for full and reduced models
df_pheno_test_full <- prepare_data(pgs_full_idlpl, df_pheno_all)
df_pheno_test_reduced <- prepare_data(pgs_reduced, df_pheno_all)


# Split data into age bins
# Subset data based on age bins for full dataset
df_bin1_full <- subset(df_pheno_test_full, age >= 37 & age <= 45)
df_bin2_full <- subset(df_pheno_test_full, age >= 45 & age <= 55)
df_bin3_full <- subset(df_pheno_test_full, age >= 55 & age <= 65)
df_bin4_full <- subset(df_pheno_test_full, age >= 65 & age <= 75)

# Subset data based on age bins for reduced dataset
df_bin1_reduced <- subset(df_pheno_test_reduced, age >= 37 & age <= 45)
df_bin2_reduced <- subset(df_pheno_test_reduced, age >= 45 & age <= 55)
df_bin3_reduced <- subset(df_pheno_test_reduced, age >= 55 & age <= 65)
df_bin4_reduced <- subset(df_pheno_test_reduced, age >= 65 & age <= 75)

# Get the number of rows for each subset in the full dataset
nrow_bin1_full <- nrow(df_bin1_full)
nrow_bin2_full <- nrow(df_bin2_full)
nrow_bin3_full <- nrow(df_bin3_full)
nrow_bin4_full <- nrow(df_bin4_full)

# Get the number of rows for each subset in the reduced dataset
nrow_bin1_reduced <- nrow(df_bin1_reduced)
nrow_bin2_reduced <- nrow(df_bin2_reduced)
nrow_bin3_reduced <- nrow(df_bin3_reduced)
nrow_bin4_reduced <- nrow(df_bin4_reduced)

# Run Cox regressions and save results
df_save <- rbind(df_save, run_cox(df_bin1_full, "Age 37-45 Full"))
df_save <- rbind(df_save, run_cox(df_bin2_full, "Age 45-55 Full"))
df_save <- rbind(df_save, run_cox(df_bin3_full, "Age 55-65 Full"))
df_save <- rbind(df_save, run_cox(df_bin4_full, "Age 65-75 Full"))

df_save <- rbind(df_save, run_cox(df_bin1_reduced, "Age 37-45 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin2_reduced, "Age 45-55 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin3_reduced, "Age 55-65 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin4_reduced, "Age 65-75 Reduced"))


cox_model <- coxph(Surv(CENSOR_AGE, PHENOTYPE) ~ adjusted_score + sex, data = df_bin1_full)

# Write results to file
write.table(df_save, file = "/home/cjm260/pt3_Cox_results/age/idlpl_PGS_411.4.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#Unstable Angina

#Unstable Angina
################
library(dplyr)
library(tidyr)
library(data.table)
library(survival)
#UNSTABLE ANGINA
#Example Sex Stratified Cox Reg for idlplcine
#Covariate file
ukb_qc_file <- "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/others/sampleQC_fromUKB_withHeaders.txt"
df_ukb_qc <- fread(ukb_qc_file, skip = 31, data.table = FALSE)
df_ukb_qc <- df_ukb_qc[, c("#UKB_ID1", "genotyping.array", paste0("PC", 1:10))]
mapData <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/reference_files/genetic/reference_files/full_release/sampleID_map.txt")
df_ukb_qc <- merge(df_ukb_qc, mapData, by.x = '#UKB_ID1', by.y = 'UKB_sample_ID')

#Coronary atherosclerosis
df_pheno <- fread("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/PGSCatalog/PheWAS_202310/20240119_phenotyped/Phecode_411.1.csv.gz")
df_pheno <- df_pheno %>%
  select(eid, sex, age_at_recruitment, PHECODE_AgeAsTimescale, PHECODE_AgeAsTimescale_Years) %>%
  rename(idno = eid, age = age_at_recruitment, PHENOTYPE = PHECODE_AgeAsTimescale, CENSOR_AGE = PHECODE_AgeAsTimescale_Years)
# Incorporate PCs and array info
df_pheno_all <- merge(df_pheno, df_ukb_qc, by.x = 'idno', by.y = 'Adiposity_sample_ID')

#Full Model PGS File
pgs_file <- "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/omics_PGS_scores/UKB_Nightingale_5e8/UKB_Nightingale.sscore.gz" 
pgs_full <- fread(pgs_file, sep = '\t')
# Rename columns for the PGS file of all the UKB samples
colnames(pgs_full) <- c('IID', sapply(colnames(pgs_full)[-1], function(col) gsub("\\.pct", "_", strsplit(col, "_")[[1]][2])))

# Only analyzing idlpl
pgs_full_idlpl <- pgs_full %>%
  select(IID, idlpl)

# Reduced model PGS File for idlpl
pgs_reduced <- fread("/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/workingdirectory/idlpl/UKB_Nightingale_idlpl.sscore.gz")
pgs_reduced <- pgs_reduced %>%
  rename(idlpl = score_sum)

# Data frame to save results
df_save <- data.frame(Regression = character(), Trait = character(), HR = numeric(), HR_low = numeric(), HR_high = numeric(), pvalue = numeric(), zscore = numeric(), cindex = numeric(), N_Case = integer(), N_All = integer(), stringsAsFactors = FALSE)

run_cox <- function(df, model_name) {
  cox_model <- coxph(Surv(CENSOR_AGE, PHENOTYPE) ~ adjusted_score + sex, data = df)
  summary_cox <- summary(cox_model)
  c_index <- summary_cox$concordance[1]
  hr <- summary_cox$coefficients[1, 2]
  hr_low <- summary_cox$conf.int[,"lower .95"]
  hr_high <- summary_cox$conf.int[,"upper .95"]
  pvalue <- summary_cox$coefficients[1,5]
  zscore <- summary_cox$coefficients[1,4]
  n_cases <- sum(df$PHENOTYPE == 1)
  n_all <- nrow(df)
  
  df_result <- data.frame(Regression = model_name, Trait = "idlpl", HR = hr, HR_low = hr_low, HR_high = hr_high, pvalue = pvalue, zscore = zscore, cindex = c_index, N_Case = n_cases, N_All = n_all)
  return(df_result)
}

# Function to prepare data, standardize, and adjust PGS for PCs
# Function to prepare data, standardize, and adjust PGS for PCs
prepare_data <- function(df_pgs, df_pheno_all) {
  # Check for missing values in key columns before merge
  if (any(is.na(df_pgs$IID)) | any(is.na(df_pheno_all$idno))) {
    stop("Missing values found in key columns 'IID' or 'idno'. Please remove or impute missing values.")
  }
  
  # Ensure key columns are of the same type
  df_pgs$IID <- as.character(df_pgs$IID)
  df_pheno_all$idno <- as.character(df_pheno_all$idno)
  
  # Merge the data frames
  df_pheno_test <- merge(df_pgs, df_pheno_all, by.x = 'IID', by.y = 'idno')
  score <- colnames(df_pgs)[2]
  
  df_pheno_test <- df_pheno_test %>%
    select(IID, all_of(score), sex, PHENOTYPE, CENSOR_AGE, age, genotyping.array, paste0('PC', 1:10))
  
  #remove rows with NAs
  df_pheno_test = na.omit(df_pheno_test)
  
  print(df_pheno_test)
  # Standardize continuous variables
  df_pheno_test[,2] = scale(df_pheno_test[,2])
  df_pheno_test$PC1 = scale(df_pheno_test$PC1)
  df_pheno_test$PC2 = scale(df_pheno_test$PC2)
  df_pheno_test$PC3 = scale(df_pheno_test$PC3)
  df_pheno_test$PC4 = scale(df_pheno_test$PC4)
  df_pheno_test$PC5 = scale(df_pheno_test$PC5)
  df_pheno_test$PC6 = scale(df_pheno_test$PC6)
  df_pheno_test$PC7 = scale(df_pheno_test$PC7)
  df_pheno_test$PC8 = scale(df_pheno_test$PC8)
  df_pheno_test$PC9 = scale(df_pheno_test$PC9)
  df_pheno_test$PC10 = scale(df_pheno_test$PC10)
  
  # Adjust PGS for PCs and genotyping array
  reg_trait <- as.formula(paste0(score, ' ~ ', paste0('PC', 1:10, collapse = ' + '), ' + as.factor(genotyping.array)'))
  adj_result <- lm(reg_trait, data = df_pheno_test)
  df_pheno_test$adjusted_score <- residuals(adj_result)
  
  
  return(df_pheno_test)
}

# Prepare data for full and reduced models
df_pheno_test_full <- prepare_data(pgs_full_idlpl, df_pheno_all)
df_pheno_test_reduced <- prepare_data(pgs_reduced, df_pheno_all)


# Split data into age bins
# Subset data based on age bins for full dataset
df_bin1_full <- subset(df_pheno_test_full, age >= 37 & age <= 45)
df_bin2_full <- subset(df_pheno_test_full, age >= 45 & age <= 55)
df_bin3_full <- subset(df_pheno_test_full, age >= 55 & age <= 65)
df_bin4_full <- subset(df_pheno_test_full, age >= 65 & age <= 75)

# Subset data based on age bins for reduced dataset
df_bin1_reduced <- subset(df_pheno_test_reduced, age >= 37 & age <= 45)
df_bin2_reduced <- subset(df_pheno_test_reduced, age >= 45 & age <= 55)
df_bin3_reduced <- subset(df_pheno_test_reduced, age >= 55 & age <= 65)
df_bin4_reduced <- subset(df_pheno_test_reduced, age >= 65 & age <= 75)

# Get the number of rows for each subset in the full dataset
nrow_bin1_full <- nrow(df_bin1_full)
nrow_bin2_full <- nrow(df_bin2_full)
nrow_bin3_full <- nrow(df_bin3_full)
nrow_bin4_full <- nrow(df_bin4_full)

# Get the number of rows for each subset in the reduced dataset
nrow_bin1_reduced <- nrow(df_bin1_reduced)
nrow_bin2_reduced <- nrow(df_bin2_reduced)
nrow_bin3_reduced <- nrow(df_bin3_reduced)
nrow_bin4_reduced <- nrow(df_bin4_reduced)

# Run Cox regressions and save results
df_save <- rbind(df_save, run_cox(df_bin1_full, "Age 37-45 Full"))
df_save <- rbind(df_save, run_cox(df_bin2_full, "Age 45-55 Full"))
df_save <- rbind(df_save, run_cox(df_bin3_full, "Age 55-65 Full"))
df_save <- rbind(df_save, run_cox(df_bin4_full, "Age 65-75 Full"))

df_save <- rbind(df_save, run_cox(df_bin1_reduced, "Age 37-45 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin2_reduced, "Age 45-55 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin3_reduced, "Age 55-65 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin4_reduced, "Age 65-75 Reduced"))



# Write results to file
write.table(df_save, file = "/home/cjm260/pt3_Cox_results/age/idlpl_PGS_411.1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



#Essential Hypertension

#Essential Hypertension


############## 
library(dplyr)
library(tidyr)
library(data.table)
library(survival)
#ESSENTIAL HYPERTENSION
#Example Sex Stratified Cox Reg for idlplcine
#Covariate file
ukb_qc_file <- "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/others/sampleQC_fromUKB_withHeaders.txt"
df_ukb_qc <- fread(ukb_qc_file, skip = 31, data.table = FALSE)
df_ukb_qc <- df_ukb_qc[, c("#UKB_ID1", "genotyping.array", paste0("PC", 1:10))]
mapData <- fread("/rds/project/jmmh2/rds-jmmh2-post_qc_data/uk_biobank/reference_files/genetic/reference_files/full_release/sampleID_map.txt")
df_ukb_qc <- merge(df_ukb_qc, mapData, by.x = '#UKB_ID1', by.y = 'UKB_sample_ID')

#Coronary atherosclerosis
df_pheno <- fread("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/PGSCatalog/PheWAS_202310/20240119_phenotyped/Phecode_401.1.csv.gz")
df_pheno <- df_pheno %>%
  select(eid, sex, age_at_recruitment, PHECODE_AgeAsTimescale, PHECODE_AgeAsTimescale_Years) %>%
  rename(idno = eid, age = age_at_recruitment, PHENOTYPE = PHECODE_AgeAsTimescale, CENSOR_AGE = PHECODE_AgeAsTimescale_Years)
# Incorporate PCs and array info
df_pheno_all <- merge(df_pheno, df_ukb_qc, by.x = 'idno', by.y = 'Adiposity_sample_ID')

#Full Model PGS File
pgs_file <- "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/omics_PGS_scores/UKB_Nightingale_5e8/UKB_Nightingale.sscore.gz" 
pgs_full <- fread(pgs_file, sep = '\t')
# Rename columns for the PGS file of all the UKB samples
colnames(pgs_full) <- c('IID', sapply(colnames(pgs_full)[-1], function(col) gsub("\\.pct", "_", strsplit(col, "_")[[1]][2])))

# Only analyzing idlpl
pgs_full_idlpl <- pgs_full %>%
  select(IID, idlpl)

# Reduced model PGS File for idlpl
pgs_reduced <- fread("/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/workingdirectory/idlpl/UKB_Nightingale_idlpl.sscore.gz")
pgs_reduced <- pgs_reduced %>%
  rename(idlpl = score_sum)

# Data frame to save results
df_save <- data.frame(Regression = character(), Trait = character(), HR = numeric(), HR_low = numeric(), HR_high = numeric(), pvalue = numeric(), zscore = numeric(), cindex = numeric(), N_Case = integer(), N_All = integer(), stringsAsFactors = FALSE)

run_cox <- function(df, model_name) {
  cox_model <- coxph(Surv(CENSOR_AGE, PHENOTYPE) ~ adjusted_score + sex, data = df)
  summary_cox <- summary(cox_model)
  c_index <- summary_cox$concordance[1]
  hr <- summary_cox$coefficients[1, 2]
  hr_low <- summary_cox$conf.int[,"lower .95"]
  hr_high <- summary_cox$conf.int[,"upper .95"]
  pvalue <- summary_cox$coefficients[1,5]
  zscore <- summary_cox$coefficients[1,4]
  n_cases <- sum(df$PHENOTYPE == 1)
  n_all <- nrow(df)
  
  df_result <- data.frame(Regression = model_name, Trait = "idlpl", HR = hr, HR_low = hr_low, HR_high = hr_high, pvalue = pvalue, zscore = zscore, cindex = c_index, N_Case = n_cases, N_All = n_all)
  return(df_result)
}

# Function to prepare data, standardize, and adjust PGS for PCs
# Function to prepare data, standardize, and adjust PGS for PCs
prepare_data <- function(df_pgs, df_pheno_all) {
  # Check for missing values in key columns before merge
  if (any(is.na(df_pgs$IID)) | any(is.na(df_pheno_all$idno))) {
    stop("Missing values found in key columns 'IID' or 'idno'. Please remove or impute missing values.")
  }
  
  # Ensure key columns are of the same type
  df_pgs$IID <- as.character(df_pgs$IID)
  df_pheno_all$idno <- as.character(df_pheno_all$idno)
  
  # Merge the data frames
  df_pheno_test <- merge(df_pgs, df_pheno_all, by.x = 'IID', by.y = 'idno')
  score <- colnames(df_pgs)[2]
  
  df_pheno_test <- df_pheno_test %>%
    select(IID, all_of(score), sex, PHENOTYPE, CENSOR_AGE, age, genotyping.array, paste0('PC', 1:10))
  
  #remove rows with NAs
  df_pheno_test = na.omit(df_pheno_test)
  
  print(df_pheno_test)
  # Standardize continuous variables
  df_pheno_test[,2] = scale(df_pheno_test[,2])
  df_pheno_test$PC1 = scale(df_pheno_test$PC1)
  df_pheno_test$PC2 = scale(df_pheno_test$PC2)
  df_pheno_test$PC3 = scale(df_pheno_test$PC3)
  df_pheno_test$PC4 = scale(df_pheno_test$PC4)
  df_pheno_test$PC5 = scale(df_pheno_test$PC5)
  df_pheno_test$PC6 = scale(df_pheno_test$PC6)
  df_pheno_test$PC7 = scale(df_pheno_test$PC7)
  df_pheno_test$PC8 = scale(df_pheno_test$PC8)
  df_pheno_test$PC9 = scale(df_pheno_test$PC9)
  df_pheno_test$PC10 = scale(df_pheno_test$PC10)
  
  # Adjust PGS for PCs and genotyping array
  reg_trait <- as.formula(paste0(score, ' ~ ', paste0('PC', 1:10, collapse = ' + '), ' + as.factor(genotyping.array)'))
  adj_result <- lm(reg_trait, data = df_pheno_test)
  df_pheno_test$adjusted_score <- residuals(adj_result)
  
  
  return(df_pheno_test)
}

# Prepare data for full and reduced models
df_pheno_test_full <- prepare_data(pgs_full_idlpl, df_pheno_all)
df_pheno_test_reduced <- prepare_data(pgs_reduced, df_pheno_all)


# Split data into age bins
# Subset data based on age bins for full dataset
df_bin1_full <- subset(df_pheno_test_full, age >= 37 & age <= 45)
df_bin2_full <- subset(df_pheno_test_full, age >= 45 & age <= 55)
df_bin3_full <- subset(df_pheno_test_full, age >= 55 & age <= 65)
df_bin4_full <- subset(df_pheno_test_full, age >= 65 & age <= 75)

# Subset data based on age bins for reduced dataset
df_bin1_reduced <- subset(df_pheno_test_reduced, age >= 37 & age <= 45)
df_bin2_reduced <- subset(df_pheno_test_reduced, age >= 45 & age <= 55)
df_bin3_reduced <- subset(df_pheno_test_reduced, age >= 55 & age <= 65)
df_bin4_reduced <- subset(df_pheno_test_reduced, age >= 65 & age <= 75)

# Get the number of rows for each subset in the full dataset
nrow_bin1_full <- nrow(df_bin1_full)
nrow_bin2_full <- nrow(df_bin2_full)
nrow_bin3_full <- nrow(df_bin3_full)
nrow_bin4_full <- nrow(df_bin4_full)

# Get the number of rows for each subset in the reduced dataset
nrow_bin1_reduced <- nrow(df_bin1_reduced)
nrow_bin2_reduced <- nrow(df_bin2_reduced)
nrow_bin3_reduced <- nrow(df_bin3_reduced)
nrow_bin4_reduced <- nrow(df_bin4_reduced)

# Run Cox regressions and save results
df_save <- rbind(df_save, run_cox(df_bin1_full, "Age 37-45 Full"))
df_save <- rbind(df_save, run_cox(df_bin2_full, "Age 45-55 Full"))
df_save <- rbind(df_save, run_cox(df_bin3_full, "Age 55-65 Full"))
df_save <- rbind(df_save, run_cox(df_bin4_full, "Age 65-75 Full"))

df_save <- rbind(df_save, run_cox(df_bin1_reduced, "Age 37-45 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin2_reduced, "Age 45-55 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin3_reduced, "Age 55-65 Reduced"))
df_save <- rbind(df_save, run_cox(df_bin4_reduced, "Age 65-75 Reduced"))


cox_model <- coxph(Surv(CENSOR_AGE, PHENOTYPE) ~ adjusted_score + sex, data = df_bin1_full)

# Write results to file
write.table(df_save, file = "/home/cjm260/pt3_Cox_results/age/idlpl_PGS_401.1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
