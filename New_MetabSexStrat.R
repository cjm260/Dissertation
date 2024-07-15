#Sex Stratified Metabolites
library(dplyr)
library(R.utils)
library(data.table)
library(tidyr)
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
#Bin the data into sex groups
age_bins <- list(
  bin1 = c(37, 45),
  bin2 = c(46, 54),
  bin3 = c(55, 63),
  bin4 = c(64, 72),
  bin5 = c(73, 83)
)

assign_age_bin <- function(age) {
  for (bin_name in names(age_bins)) {
    bin <- age_bins[[bin_name]]
    if (age >= bin[1] && age <= bin[2]) {
      return(bin_name)
    }
  }
  return(NA) # If age doesn't fall into any bin
}
demodata$age_bin = sapply(demodata$age, assign_age_bin)
#Sex bins
#bin 1 is female
bin1 = subset(demodata, sex == "Female")
#bin 2 is male
bin2 = subset(demodata, sex == "Male")

table(bin1$age_bin)
#3796 females in age bin 5. Resample the remaining four age bins to size 3796 each.
#"bin_1" refers to age bin 1
ab1 = bin1[bin1$age_bin == "bin1"]
sampsAB1 = ab1[sample(nrow(ab1), 3796, replace = FALSE)]

ab2 = bin1[bin1$age_bin == "bin2"]
sampsAB2 = ab2[sample(nrow(ab2), 3796, replace = FALSE)]

ab3 = bin1[bin1$age_bin == "bin3"]
sampsAB3 = ab3[sample(nrow(ab3), 3796, replace = FALSE)]

ab4 = bin1[bin1$age_bin == "bin4"]
sampsAB4 = ab4[sample(nrow(ab4), 3796, replace = FALSE)]

ab5 = bin1[bin1$age_bin == "bin5"]
#Bin1 now contains 3696 females from each age group
bin1 = rbind(sampsAB1, sampsAB2, sampsAB3, sampsAB4, ab5) #18980 rows total

#Repeat the same for the males bin2
table(bin2$age_bin)
#5221 males in bin 5 is the minimum. Resample the remaining four age bins to size 5221.
ab1 = bin2[bin2$age_bin == "bin1"]
sampsAB1 = ab1[sample(nrow(ab1), 5221, replace = FALSE)]

ab2 = bin2[bin2$age_bin == "bin2"]
sampsAB2 = ab2[sample(nrow(ab2), 5221, replace = FALSE)]

ab3 = bin2[bin2$age_bin == "bin3"]
sampsAB3 = ab3[sample(nrow(ab3), 5221, replace = FALSE)]

ab4 = bin2[bin2$age_bin == "bin4"]
sampsAB4 = ab4[sample(nrow(ab4), 5221, replace = FALSE)]

ab5 = bin2[bin2$age_bin == "bin5"]
bin2 = rbind(sampsAB1, sampsAB2, sampsAB3, sampsAB4, ab5) #26105 rows total

#Split df_merge into females and males
bin1 = bin1 %>%
  select(IID, sex, age)
bin2 = bin2 %>%
  select(IID, sex, age)
#df_merged data with only females
df_merged_fem = merge(df_merged, bin1, by = "IID")
df_merged_male = merge(df_merged, bin2, by = "IID")

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

# Initialize results dataframe for females
df_results_fem <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                         R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                         stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:220) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$interval_name
  col2 <- row$ukb_name
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_fem %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  df_one <- df_one %>% rename(trait_value = all_of(col2))
  
  ids <- unique(c(df_one$IID, ids))
  
  # Adjust for PCs
  reg_trait <- as.formula(paste('trait_value ~', paste(paste0('PC', 1:10), collapse = ' + ')))
  result <- lm(reg_trait, data = df_one)
  #df_one$trait_value_adjusted <- residuals(result)
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
  
  df_results_fem <- df_results_fem %>% add_row(INTERVAL_Name = col1, UKB_Name = col2, N_val_samples = nrow(df_one),
                                       R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_fem <- df_results_fem %>% arrange(desc(R2))


# Initialize results dataframe for males
df_results_male <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                             R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                             stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:220) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$interval_name
  col2 <- row$ukb_name
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_male %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  df_one <- df_one %>% rename(trait_value = all_of(col2))
  
  ids <- unique(c(df_one$IID, ids))
  
  # Adjust for PCs
  reg_trait <- as.formula(paste('trait_value ~', paste(paste0('PC', 1:10), collapse = ' + ')))
  result <- lm(reg_trait, data = df_one)
  #df_one$trait_value_adjusted <- residuals(result)
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
  
  df_results_male <- df_results_male %>% add_row(INTERVAL_Name = col1, UKB_Name = col2, N_val_samples = nrow(df_one),
                                           R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_male <- df_results_male %>% arrange(desc(R2))

#Save female results and male results 
fpath = "/home/cjm260/Fem_MetabSexStrat.csv"
write.csv(df_results_fem, fpath, row.names = FALSE)
mpath = "/home/cjm260/Male_MetabSexStrat.csv"
write.csv(df_results_male, mpath, row.names = FALSE)

#Slope calculation
females = fread("/home/cjm260/Fem_MetabSexStrat.csv")
rows_to_remove <- grepl("pct$", females$UKB_Name)
females =  slopeDf[!rows_to_remove, ]

males = fread("/home/cjm260/Male_MetabSexStrat.csv")

females = df_results_fem
males = df_results_male

#Prepare female and male dfs to be merged
femPrep = females %>%
  rename(R2_F = R2) %>%
  rename(R2_P_F = R2_P) %>%
  rename(SR_F = SR) %>%
  rename(SR_P_F = SR_P) %>%
  rename(N_val_samples_F = N_val_samples)
malePrep = males %>%
  rename(R2_M = R2) %>%
  rename(R2_P_M = R2_P) %>%
  rename(SR_M = SR) %>%
  rename(SR_P_M = SR_P) %>%
  rename(N_val_samples_M = N_val_samples)

sexSlopeDf = merge(femPrep, malePrep, by = "INTERVAL_Name")
sexSlopeDf = sexSlopeDf %>%
  arrange(desc(R2_F))

#Calculate the slopes of each trait and store in a vector
slopes = numeric(nrow(sexSlopeDf))
for (i in 1:nrow(sexSlopeDf)) {
  slopes[i] = abs(sexSlopeDf$SR_F[i] - sexSlopeDf$SR_M[i])
}
sexSlopeDf$slopes = slopes
sexSlopeDf = sexSlopeDf %>%
  arrange(desc(slopes))

slopePath = "/home/cjm260/SlopeDf_MetabSexStrat.csv"
write.csv(sexSlopeDf, slopePath, row.names = FALSE)

#Remove pct entries
slopeDf = fread("/home/cjm260/SlopeDf_MetabSexStrat.csv")
rows_to_remove <- grepl("pct$", slopeDf$UKB_Name.x)
df_filtered =  slopeDf[!rows_to_remove, ]

pathnew = "/home/cjm260/SlopeDfFiltered_MetabSexStrat.csv"
write.csv(df_filtered, pathnew, row.names = FALSE)

df_filtered = fread("/home/cjm260/SlopeDfFiltered_MetabSexStrat.csv")
df_filtered = df_filtered %>%
  rename(UKB_Name = UKB_Name.x)
dfJustSlope = df_filtered %>%
  select(INTERVAL_Name, UKB_Name, slopes)
pathJust = "/home/cjm260/SlopeDfUSEFIG_MetabSexStrat.csv"
write.csv(dfJustSlope, pathJust, row.names = FALSE)

