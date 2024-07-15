#Sex Stratified Olink proteins
library(dplyr)
library(tidyr)
library(data.table)


pgscores = fread("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/yx322/UKB_omics_PGS/Olink_full/UKB_Olink.sscore.gz")
colnames(pgscores) <- gsub("^olink_|_br_effects_[1-5]$", "", colnames(pgscores)) #col names in the Uniprot ID format
#all traits in the omicspred modle
traits_all = colnames(pgscores)[-1]

trait_levels = fread("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/protein_olink/EUR_wide_olink_data_2923.pheno")
trait_levels = trait_levels%>%
  select(-FID) #traits are in the pid format
#PG score and trait data
df_merged = merge(pgscores, trait_levels, by = "IID")

df_name_map = fread("/rds/project/asb38/rds-asb38-ceu-ukbiobank/projects/P7439/inouyelab/share_space/curated_ukb_data/protein_olink/olink_protein_info.txt")
df_name_map = df_name_map %>%
  select(pid, UniProt)
df_name_map <- df_name_map[!is.na(pid) & !is.na(UniProt), ]
# with significant qtls in INTERVAL and measured in UKB
df_name_map_overlap <- df_name_map[UniProt %in% traits_all, ]
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

#Split the data into sexes
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

#FEMALE CORRELATION CALC
# Initialize results dataframe
df_results_fem <- data.frame(pid = character(), UniProt = character(), R2 = numeric(),
                         R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                         stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap) = 298 protiens
for (i in 1:nrow(df_name_map_overlap)) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$pid
  col2 <- row$UniProt
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_fem %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  #In this case, pid/col1 are the trait values
  df_one <- df_one %>% rename(trait_value = all_of(col1))
  
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
  
  r2 <- cor(df_one[[col2]], df_one$trait_value_adjusted)^2
  r2_p <- cor.test(df_one[[col2]], df_one$trait_value_adjusted)$p.value
  
  sp_score <- cor(df_one[[col2]], df_one$trait_value_adjusted, method = "spearman")
  sp_score_p <- cor.test(df_one[[col2]], df_one$trait_value_adjusted, method = "pearson")$p.value
  
  df_results_fem <- df_results_fem %>% add_row(pid = col1, UniProt = col2, N_val_samples = nrow(df_one),
                                       R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_fem <- df_results_fem %>% arrange(desc(R2))


#MALE CORRELATION CALC
# Initialize results dataframe
df_results_male <- data.frame(pid = character(), UniProt = character(), R2 = numeric(),
                             R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                             stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)
for (i in 1:nrow(df_name_map_overlap)) {
  row <- df_name_map_overlap[i, ] #test on row 1 (ace, acetate)
  col1 <- row$pid
  col2 <- row$UniProt
  
  #ace col contains all the PGS, acetate col contains all the trait levels for each IDI
  df_one <- df_merged_male %>% select(IID, all_of(col1), all_of(col2)) %>% drop_na()
  #df_one <- merge(df_one, df_ukb_qc, by.x = 'IID', by.y = 'Adiposity_sample_ID')
  #alternative after renameing aidposity_sample_id to IID
  df_one = merge(df_one, df_ukb_qc2, by = 'IID')
  #In this case, pid/col1 are the trait values
  df_one <- df_one %>% rename(trait_value = all_of(col1))
  
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
  
  r2 <- cor(df_one[[col2]], df_one$trait_value_adjusted)^2
  r2_p <- cor.test(df_one[[col2]], df_one$trait_value_adjusted)$p.value
  
  sp_score <- cor(df_one[[col2]], df_one$trait_value_adjusted, method = "spearman")
  sp_score_p <- cor.test(df_one[[col2]], df_one$trait_value_adjusted, method = "pearson")$p.value
  
  df_results_male <- df_results_male %>% add_row(pid = col1, UniProt = col2, N_val_samples = nrow(df_one),
                                               R2 = r2, R2_P = r2_p, SR = sp_score, SR_P = sp_score_p)
}

df_results_male <- df_results_male %>% arrange(desc(R2))

#Save female results and male results 
fpath = "/home/cjm260/Fem_OlinkSexStrat.csv"
write.csv(df_results_fem, fpath, row.names = FALSE)
mpath = "/home/cjm260/Male_OlinkSexStrat.csv"
write.csv(df_results_male, mpath, row.names = FALSE)

#Slope calculation
females = fread("/home/cjm260/Fem_OlinkSexStrat.csv")
males = fread("/home/cjm260/Male_OlinkSexStrat.csv")

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
  rename(N_val_samples_M = N_val_samples) %>%
  select(-pid)

sexSlopeDf = merge(femPrep, malePrep, by = "UniProt")
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

slopePath = "/home/cjm260/SlopeDf_OlinkSexStrat.csv"
write.csv(sexSlopeDf, slopePath, row.names = FALSE)

#just slopes for figure
dfJustSlope = sexSlopeDf %>%
  select(UniProt, pid, slopes)

justSlopePath = "/home/cjm260/SlopeDfUSEFIG_OlinkSexStrat.csv"
write.csv(dfJustSlope, justSlopePath, row.names = FALSE)

