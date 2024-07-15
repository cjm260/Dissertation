#Age Stratified Metabolites
library(dplyr)
library(tidyr)
library(R.utils)
library(data.table)
#library(broom)
#calculated PGS gile of UKB Samples
pgs_file = "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/omics_PGS_scores/UKB_Nightingale_5e8/UKB_Nightingale.sscore.gz" 
df_ukb_pgs = fread(pgs_file, sep = '\t')

#rename columns for the PGS file of all the UKB samples
colnames(df_ukb_pgs) <- c('IID', sapply(colnames(df_ukb_pgs)[-1], function(col) gsub("\\.pct", "_", strsplit(col, "_")[[1]][2])))
traits_all <- colnames(df_ukb_pgs)[-1]

#old_trait_level_file <- "/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/yx322/Nightingale_prerelease_postqc_from_Scott/postqc_adjusted_age_sex_by_yu.txt"
old_trait_level_file = "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/others/postqc_adjusted_age_sex_by_yu.txt"
trait_level_file <- "/rds/project/asb38/rds-asb38-ceu-ukbiobank/phenotype/P7439/post_qc_data/NMR_metabolite_biomarkers/nmr_techadj.txt"

df_trait_level <- fread(trait_level_file, sep = '\t')
df_trait_level_old <- fread(old_trait_level_file, sep = '\t')

df_trait_level <- df_trait_level[df_trait_level$visit_index == 0, ]
#removes visit index
col_names <- c('eid', colnames(df_trait_level_old)[-c(1:2)])
#replaces df_trait_level with older colnames
df_trait_level <- df_trait_level[, ..col_names]
#merges pgs and trait level files together.
df_merged <- merge(df_ukb_pgs, df_trait_level, by.x = 'IID', by.y = 'eid')

#demographic data
demodata = fread("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/anthropometrics/output/anthropometrics.txt")
demodata = demodata %>%
  rename(IID = eid)


#Bin the data into age groups
bin1 = subset(demodata, age >= 37 & age < 46)
table(bin1$sex)
Fems1 = bin1[bin1$sex == "Female",]
#Randomly sample females up to the number of males in bin1
sampsFems1 = Fems1[sample(nrow(Fems1), 29728, replace = FALSE), ]
males1 = bin1[bin1$sex == "Male", ]
#bin now has the same number of males and females
bin1 = rbind(sampsFems1, males1)

bin2 = subset(demodata, age >= 46 & age < 55)
table(bin2$sex)
Fems2 = bin2[bin2$sex == "Female",]
#Randomly sample females up to the number of males in bin2
sampsFems2 = Fems2[sample(nrow(Fems2), 61595, replace = FALSE), ]
males2 = bin2[bin2$sex == "Male", ]
#bin now has the same number of males and females
bin2 = rbind(sampsFems2, males2)


bin3 = subset(demodata, age >= 55 & age < 64)
table(bin3$sex)
Fems3 = bin3[bin3$sex == "Female",]
#Randomly sample females up to the number of males in bin3
sampsFems3 = Fems3[sample(nrow(Fems3), 95363, replace = FALSE), ]
males3 = bin3[bin3$sex == "Male", ]
#bin now has the same number of males and females
bin3 = rbind(sampsFems3, males3)


bin4 = subset(demodata, age >= 64 & age < 73)
table(bin4$sex)
Fems4 = bin4[bin4$sex == "Female",]
#Randomly sample females up to the number of males in bin4
sampsFems4 = Fems4[sample(nrow(Fems4), 73033, replace = FALSE), ]
males4 = bin4[bin4$sex == "Male", ]
#bin now has the same number of males and females
bin4 = rbind(sampsFems4, males4)


#bin5 = subset(demodata, age >= 73 & age < 83)
#table(bin5$sex)
#males5 = bin5[bin5$sex == "Male",]
#Randomly sample females up to the number of males in bin5
#sampsMales5 = males5[sample(nrow(males5), 3796, replace = FALSE), ]
#fems5 = bin5[bin5$sex == "Female", ]
#bin now has the same number of males and females
#bin5 = rbind(sampsMales5, fems5)


#Split df_merged into age groups
#Stratify the trait data by age bins
#only merge eid, age to the trait data
bin1cols = bin1[,c("IID", "age")]
#Merge the trait data corresponding to the individuals in each bin
df_merged_1 = merge(df_merged, bin1cols, by = "IID")  #nrow = 33480

bin2cols = bin2[,c("IID", "age")]
#Merge the trait data corresponding to the individuals in each bin
df_merged_2 = merge(df_merged, bin2cols, by = "IID")

bin3cols = bin3[,c("IID", "age")]
#Merge the trait data corresponding to the individuals in each bin
df_merged_3 = merge(df_merged, bin3cols, by = "IID") #nrow = 122838

bin4cols = bin4[,c("IID", "age")]
#Merge the trait data corresponding to the individuals in each bin
df_merged_4 = merge(df_merged, bin4cols, by = "IID") #nrow = 90188

#bin5cols = bin5[,c("IID", "age")]
#Merge the trait data corresponding to the individuals in each bin
#df_merged_5 = merge(df_merged, bin5cols, by = "IID") #only 6570


# read ukb vs interval trait name mapping
name_map_file <- "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_prs/others/interval_ukb_nmr_map.csv"
#df_name_map <- fread(name_map_file, sep = ',', encoding = "ISO-8859-1")
df_name_map <- fread(name_map_file, sep = ',', encoding = "Latin-1")
df_name_map <- df_name_map[, .(interval_name = `INTERVAL column name`, ukb_name = `UKB column name`)]
df_name_map <- df_name_map[!is.na(interval_name) & !is.na(ukb_name), ]

# with significant qtls in INTERVAL and measured in UKB
df_name_map_overlap <- df_name_map[interval_name %in% traits_all, ]
# get UKB sample pcs data
#ukb_qc_file <- "/home/yx322/rds/rds-jmmh2-post_qc_data/uk_biobank/reference_files/genetic/reference_files/full_release/QC_documents/sampleQC_fromUKB_withHeaders.txt"
ukb_qc_file = "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/others/sampleQC_fromUKB_withHeaders.txt"
# UKB sample ID mapings
#ukb_id_map_file = "/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/others/sampleID_map.txt "
# read ukb samples qc data to get array and PCs info
df_ukb_qc <- fread(ukb_qc_file, skip = 31, sep = ' ')


df_ukb_id_map = fread("/home/cjm260/sampleID_map.txt", sep = '\t')
# Merge to get Adiposity sample ID
df_ukb_qc <- merge(df_ukb_id_map, df_ukb_qc, by.x = 'UKB_sample_ID', by.y = '#UKB_ID1')
df_ukb_qc <- df_ukb_qc[df_ukb_qc$in.white.British.ancestry.subset == 1, ]
df_ukb_qc <- df_ukb_qc[, c('Adiposity_sample_ID', paste0('PC', 1:10))]

df_ukb_qc2 = df_ukb_qc %>%
  rename(IID = Adiposity_sample_ID)

# Initialize results dataframe for BIN 1
df_results_1 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                             R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                             stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:220) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$interval_name
  col2 <- row$ukb_name
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_1 %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  df_one <- df_one %>% rename(trait_value = all_of(col2))
  
  ids <- unique(c(df_one$IID, ids))
  
  # Adjust for PCs
  reg_trait <- as.formula(paste('trait_value ~', paste(paste0('PC', 1:10), collapse = ' + ')))
  result <- lm(reg_trait, data = df_one)
  resids = residuals(result)
  ranks <- rank(resids, ties.method = "average")
  # Step 4: Transform ranks to a normal distribution
  transformed_residuals <- qnorm((ranks - 0.5) / length(ranks))
  #attach inverse rank normalized residuals back to original df
  df_one$trait_value_adjusted = transformed_residuals
  
  r2 <- cor(df_one[[col1]], df_one$trait_value_adjusted)^2
  r2_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted)$p.value
  
  sp_score <- cor(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman")
  sp_score_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman", exact = FALSE)$p.value
  
  df_results_1 <- df_results_1 %>% add_row(INTERVAL_Name = col1, UKB_Name = col2, N_val_samples = nrow(df_one),
                                           R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_1 <- df_results_1 %>% arrange(desc(R2))

#BIN 2
#BIN 2
# Initialize results dataframe for BIN 2
df_results_2 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                           R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                           stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:220) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$interval_name
  col2 <- row$ukb_name
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_2 %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  df_one <- df_one %>% rename(trait_value = all_of(col2))
  
  ids <- unique(c(df_one$IID, ids))
  
  # Adjust for PCs
  reg_trait <- as.formula(paste('trait_value ~', paste(paste0('PC', 1:10), collapse = ' + ')))
  result <- lm(reg_trait, data = df_one)
  resids = residuals(result)
  ranks <- rank(resids, ties.method = "average")
  # Step 4: Transform ranks to a normal distribution
  transformed_residuals <- qnorm((ranks - 0.5) / length(ranks))
  #attach inverse rank normalized residuals back to original df
  df_one$trait_value_adjusted = transformed_residuals
  
  r2 <- cor(df_one[[col1]], df_one$trait_value_adjusted)^2
  r2_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted)$p.value
  
  sp_score <- cor(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman")
  sp_score_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman", exact = FALSE)$p.value
  
  df_results_2 <- df_results_2 %>% add_row(INTERVAL_Name = col1, UKB_Name = col2, N_val_samples = nrow(df_one),
                                           R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_2 <- df_results_2 %>% arrange(desc(R2))


#BIN 3
# Initialize results dataframe for BIN 3
df_results_3 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                           R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                           stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:220) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$interval_name
  col2 <- row$ukb_name
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_3 %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  df_one <- df_one %>% rename(trait_value = all_of(col2))
  
  ids <- unique(c(df_one$IID, ids))
  
  # Adjust for PCs
  reg_trait <- as.formula(paste('trait_value ~', paste(paste0('PC', 1:10), collapse = ' + ')))
  result <- lm(reg_trait, data = df_one)
  resids = residuals(result)
  ranks <- rank(resids, ties.method = "average")
  # Step 4: Transform ranks to a normal distribution
  transformed_residuals <- qnorm((ranks - 0.5) / length(ranks))
  #attach inverse rank normalized residuals back to original df
  df_one$trait_value_adjusted = transformed_residuals
  
  r2 <- cor(df_one[[col1]], df_one$trait_value_adjusted)^2
  r2_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted)$p.value
  
  sp_score <- cor(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman")
  sp_score_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman", exact = FALSE)$p.value
  
  df_results_3 <- df_results_3 %>% add_row(INTERVAL_Name = col1, UKB_Name = col2, N_val_samples = nrow(df_one),
                                           R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_3 <- df_results_3 %>% arrange(desc(R2))


#BIN 4
# Initialize results dataframe for BIN 4
df_results_4 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                           R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                           stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:220) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$interval_name
  col2 <- row$ukb_name
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_4 %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  df_one <- df_one %>% rename(trait_value = all_of(col2))
  
  ids <- unique(c(df_one$IID, ids))
  
  # Adjust for PCs
  reg_trait <- as.formula(paste('trait_value ~', paste(paste0('PC', 1:10), collapse = ' + ')))
  result <- lm(reg_trait, data = df_one)
  resids = residuals(result)
  ranks <- rank(resids, ties.method = "average")
  # Step 4: Transform ranks to a normal distribution
  transformed_residuals <- qnorm((ranks - 0.5) / length(ranks))
  #attach inverse rank normalized residuals back to original df
  df_one$trait_value_adjusted = transformed_residuals
  
  r2 <- cor(df_one[[col1]], df_one$trait_value_adjusted)^2
  r2_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted)$p.value
  
  sp_score <- cor(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman")
  sp_score_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman", exact = FALSE)$p.value
  
  df_results_4 <- df_results_4 %>% add_row(INTERVAL_Name = col1, UKB_Name = col2, N_val_samples = nrow(df_one),
                                           R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_4 <- df_results_4 %>% arrange(desc(R2))

#BIN 5
# Initialize results dataframe for BIN 5
df_results_5 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                           R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                           stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:220) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$interval_name
  col2 <- row$ukb_name
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_5 %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  df_one <- df_one %>% rename(trait_value = all_of(col2))
  
  ids <- unique(c(df_one$IID, ids))
  
  # Adjust for PCs
  reg_trait <- as.formula(paste('trait_value ~', paste(paste0('PC', 1:10), collapse = ' + ')))
  result <- lm(reg_trait, data = df_one)
  df_one$trait_value_adjusted <- residuals(result)
  
  r2 <- cor(df_one[[col1]], df_one$trait_value_adjusted)^2
  r2_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted)$p.value
  
  sp_score <- cor(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman")
  sp_score_p <- cor.test(df_one[[col1]], df_one$trait_value_adjusted, method = "spearman", exact = FALSE)$p.value
  
  df_results_5 <- df_results_5 %>% add_row(INTERVAL_Name = col1, UKB_Name = col2, N_val_samples = nrow(df_one),
                                           R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_5 <- df_results_5 %>% arrange(desc(R2))

#Save age bin results
path1 = "/home/cjm260/B1_MetabAgeStrat.csv"
write.csv(df_results_1, path1, row.names = FALSE)
path2 = "/home/cjm260/B2_MetabAgeStrat.csv"
write.csv(df_results_2, path2, row.names = FALSE)
path3 = "/home/cjm260/B3_MetabAgeStrat.csv"
write.csv(df_results_3, path3, row.names = FALSE)
path4 = "/home/cjm260/B4_MetabAgeStrat.csv"
write.csv(df_results_4, path4, row.names = FALSE)
#path5 = "/home/cjm260/B5_MetabAgeStrat.csv"
#write.csv(df_results_5, path5, row.names = FALSE)

#SLOPE CALCULATION SECTION
B1 = fread("/home/cjm260/B1_MetabAgeStrat.csv")
B2 = fread("/home/cjm260/B2_MetabAgeStrat.csv")
B3 = fread("/home/cjm260/B3_MetabAgeStrat.csv")
B4 = fread("/home/cjm260/B4_MetabAgeStrat.csv")
#B5 = fread("/home/cjm260/B5_MetabAgeStrat.csv")

#Prepare bin dfs to be merged
B1Prep = B1 %>%
  rename(R2_B1 = R2) %>%
  rename(R2_P_B1 = R2_P) %>%
  rename(SR_B1 = SR) %>%
  rename(SR_P_B1 = SR_P) %>%
  rename(N_val_samples_B1 = N_val_samples) 
B2Prep = B2 %>%
  rename(R2_B2 = R2) %>%
  rename(R2_P_B2 = R2_P) %>%
  rename(SR_B2 = SR) %>%
  rename(SR_P_B2 = SR_P) %>%
  rename(N_val_samples_B2 = N_val_samples) %>%
  select(-UKB_Name)
B3Prep = B3 %>%
  rename(R2_B3 = R2) %>%
  rename(R2_P_B3 = R2_P) %>%
  rename(SR_B3 = SR) %>%
  rename(SR_P_B3 = SR_P) %>%
  rename(N_val_samples_B3 = N_val_samples) %>%
  select(-UKB_Name)
B4Prep = B4 %>%
  rename(R2_B4 = R2) %>%
  rename(R2_P_B4 = R2_P) %>%
  rename(SR_B4 = SR) %>%
  rename(SR_P_B4 = SR_P) %>%
  rename(N_val_samples_B4 = N_val_samples) %>%
  select(-UKB_Name)
#B5Prep = B5 %>%
  #rename(R2_B5 = R2) %>%
  #rename(R2_P_B5 = R2_P) %>%
  #rename(SR_B5 = SR) %>%
  #rename(SR_P_B5 = SR_P) %>%
  #rename(N_val_samples_B5 = N_val_samples) %>%
  #select(-UKB_Name)


#merge all the dfs together
#OR TO MAKE IT EASIER TO READ SHOULD I HAVE SEPARATE DFS SO I ONLY HAVE ONE SLOPE COLUMN?
m1 = merge(B1Prep, B2Prep, by = "INTERVAL_Name")
m2 = merge(m1, B3Prep, by = "INTERVAL_Name")
merged = merge(m2, B4Prep, by = "INTERVAL_Name")
#merged = merge(m3, B5Prep, by = "INTERVAL_Name")
ageSlopeDf = merged %>%
  select(INTERVAL_Name, UKB_Name, SR_B1, SR_P_B1, N_val_samples_B1, SR_B2, SR_P_B2, N_val_samples_B2, SR_B3, SR_P_B3, N_val_samples_B3, SR_B4, SR_P_B4, N_val_samples_B4) %>%
  arrange(desc(SR_B1))

#Initialize dfs to hold slopes
s1_2 = numeric(nrow(ageSlopeDf))
s2_3 = numeric(nrow(ageSlopeDf))
s3_4 = numeric(nrow(ageSlopeDf))
#s4_5 = numeric(nrow(ageSlopeDf))
avg_slope = numeric(nrow(ageSlopeDf))

for (i in 1:nrow(ageSlopeDf)) {
  s1_2[i] = abs(ageSlopeDf$SR_B2[i] - ageSlopeDf$SR_B1[i])
  s2_3[i] = abs(ageSlopeDf$SR_B3[i] - ageSlopeDf$SR_B2[i])
  s3_4[i] = abs(ageSlopeDf$SR_B4[i] - ageSlopeDf$SR_B3[i])
  #4_5[i] = abs(ageSlopeDf$SR_B5[i] - ageSlopeDf$SR_B4[i])
  avg_slope[i] = (s1_2[i] + s2_3[i] + s3_4[i])/3
}

ageSlopeDf$Slope_B1_B2 = s1_2
ageSlopeDf$Slope_B2_B3 = s2_3
ageSlopeDf$Slope_B3_B4 = s3_4
#ageSlopeDf$Slope_B4_B5 = s4_5
ageSlopeDf$Avg_Slope = avg_slope

ageSlopeDf = ageSlopeDf %>%
  arrange(desc(Avg_Slope))

#slopePath = "/home/cjm260/SlopeDf_MetabAgeStrat.csv"
#write.csv(ageSlopeDf, slopePath, row.names = FALSE)

#Remove pct entries

rows_to_remove <- grepl("pct$", ageSlopeDf$UKB_Name)
df_filtered =  ageSlopeDf[!rows_to_remove, ]

#THIS FULL DF WILL BE IN THE APPENDIX
pathnew = "/home/cjm260/SlopeDfFiltered_MetabAgeStrat.csv"
write.csv(df_filtered, pathnew, row.names = FALSE)


dfJustSlope = df_filtered %>%
  select(INTERVAL_Name, UKB_Name, Slope_B1_B2, Slope_B2_B3, Slope_B3_B4, Avg_Slope)
pathJust = "/home/cjm260/SlopeDfUSEFIG_MetabAgeStrat.csv"
write.csv(dfJustSlope, pathJust, row.names = FALSE)




