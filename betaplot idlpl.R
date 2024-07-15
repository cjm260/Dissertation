library(data.table)
library(dplyr)
library(qqman)

#bin2,3 = 339 rows, bin1 = 337 rows

#idlpl bin 1

chr1 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr1.invRankNorm_resids.glm.linear")
chr2 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr2.invRankNorm_resids.glm.linear")
#chr3 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr3.invRankNorm_resids.glm.linear")
#chr4 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr4.invRankNorm_resids.glm.linear")
chr5 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr5.invRankNorm_resids.glm.linear")
#chr6 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr6.invRankNorm_resids.glm.linear")
chr7 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr7.invRankNorm_resids.glm.linear")
chr8 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr8.invRankNorm_resids.glm.linear")
chr9 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr9.invRankNorm_resids.glm.linear")
#chr10 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr10.invRankNorm_resids.glm.linear")
chr11 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr11.invRankNorm_resids.glm.linear")
#chr12 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr12.invRankNorm_resids.glm.linear")
#chr13 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr13.invRankNorm_resids.glm.linear")
#chr14 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr14.invRankNorm_resids.glm.linear")
chr15 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr15.invRankNorm_resids.glm.linear")
#chr16 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr16.invRankNorm_resids.glm.linear")
chr17 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr17.invRankNorm_resids.glm.linear")
#chr18 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr18.invRankNorm_resids.glm.linear")
chr19 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr19.invRankNorm_resids.glm.linear")
chr20 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr20.invRankNorm_resids.glm.linear")
#chr21 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr21.invRankNorm_resids.glm.linear")
#chr22 = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/gwas_chr22.invRankNorm_resids.glm.linear")

#Rbinding them all together
chr_list <- list(chr1, chr2, chr5, chr7, chr8, chr9, chr11, chr15, chr17, chr19, chr20)


# Use do.call with rbind to bind all data frames together
bin1 <- do.call(rbind, chr_list) #260 SNPs. All look significant

#open path to print the plot to
path = "/home/cjm260/gwas_files/gwas_results/idlpl1_ManhatPlot.png"
png(path)

manhattan(bin1, chr="#CHROM", bp="POS", snp="ID", p="P", annotatePval = 0.01, main = "SNPs Associated with Phospholipids in IDL in Age Bin 1 INTERVAL Samples")

#fastman(bin1, chr="#CHROM", bp="POS", snp="ID", p="P",speedup=TRUE, logp = TRUE, scattermore = FALSE,annotatePval = 0.01, main = "SNPs Associated with Free Cholesterol in Very Large HDL in Male INTERVAL Samples")

dev.off()


#############idlpl bin 2


chr1 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr1.invRankNorm_resids.glm.linear")
chr2 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr2.invRankNorm_resids.glm.linear")
#chr3 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr3.invRankNorm_resids.glm.linear")
#chr4 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr4.invRankNorm_resids.glm.linear")
chr5 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr5.invRankNorm_resids.glm.linear")
#chr6 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr6.invRankNorm_resids.glm.linear")
chr7 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr7.invRankNorm_resids.glm.linear")
chr8 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr8.invRankNorm_resids.glm.linear")
chr9 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr9.invRankNorm_resids.glm.linear")
#chr10 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr10.invRankNorm_resids.glm.linear")
chr11 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr11.invRankNorm_resids.glm.linear")
#chr12 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr12.invRankNorm_resids.glm.linear")
#chr13 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr13.invRankNorm_resids.glm.linear")
#chr14 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr14.invRankNorm_resids.glm.linear")
chr15 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr15.invRankNorm_resids.glm.linear")
#chr16 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr16.invRankNorm_resids.glm.linear")
chr17 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr17.invRankNorm_resids.glm.linear")
#chr18 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr18.invRankNorm_resids.glm.linear")
chr19 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr19.invRankNorm_resids.glm.linear")
chr20 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr20.invRankNorm_resids.glm.linear")
#chr21 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr21.invRankNorm_resids.glm.linear")
#chr22 = fread("/home/cjm260/gwas_files/gwas_results/idlpl2/gwas_chr22.invRankNorm_resids.glm.linear")

#Rbinding them all together
chr_list <- list(chr1, chr2, chr5, chr7, chr8, chr9, chr11, chr15, chr17, chr19, chr20)


# Use do.call with rbind to bind all data frames together
bin2 <- do.call(rbind, chr_list) #260 SNPs. All look significant

#open path to print the plot to
path = "/home/cjm260/gwas_files/gwas_results/idlpl2_ManhatFemPlot.png"
png(path)

manhattan(bin2, chr="#CHROM", bp="POS", snp="ID", p="P", annotatePval = 0.01, main = "SNPs Associated with Phospholipids in IDL in Age Bin 2 INTERVAL Samples")

#fastman(bin1, chr="#CHROM", bp="POS", snp="ID", p="P",speedup=TRUE, logp = TRUE, scattermore = FALSE,annotatePval = 0.01, main = "SNPs Associated with Free Cholesterol in Very Large HDL in Male INTERVAL Samples")

dev.off()


#############idlpl bin 3


chr1 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr1.invRankNorm_resids.glm.linear")
chr2 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr2.invRankNorm_resids.glm.linear")
#chr3 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr3.invRankNorm_resids.glm.linear")
#chr4 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr4.invRankNorm_resids.glm.linear")
chr5 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr5.invRankNorm_resids.glm.linear")
#chr6 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr6.invRankNorm_resids.glm.linear")
chr7 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr7.invRankNorm_resids.glm.linear")
chr8 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr8.invRankNorm_resids.glm.linear")
chr9 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr9.invRankNorm_resids.glm.linear")
#chr10 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr10.invRankNorm_resids.glm.linear")
chr11 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr11.invRankNorm_resids.glm.linear")
#chr12 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr12.invRankNorm_resids.glm.linear")
#chr13 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr13.invRankNorm_resids.glm.linear")
#chr14 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr14.invRankNorm_resids.glm.linear")
chr15 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr15.invRankNorm_resids.glm.linear")
#chr16 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr16.invRankNorm_resids.glm.linear")
chr17 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr17.invRankNorm_resids.glm.linear")
#chr18 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr18.invRankNorm_resids.glm.linear")
chr19 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr19.invRankNorm_resids.glm.linear")
chr20 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr20.invRankNorm_resids.glm.linear")
#chr21 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr21.invRankNorm_resids.glm.linear")
#chr22 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr22.invRankNorm_resids.glm.linear")

#Rbinding them all together
chr_list <- list(chr1, chr2, chr5, chr7, chr8, chr9, chr11, chr15, chr17, chr19, chr20)


# Use do.call with rbind to bind all data frames together
bin3 <- do.call(rbind, chr_list) #260 SNPs. All look significant

#open path to print the plot to
path = "/home/cjm260/gwas_files/gwas_results/idlpl3_ManhatFemPlot.png"
png(path)

manhattan(bin3, chr="#CHROM", bp="POS", snp="ID", p="P", annotatePval = 0.01, main = "SNPs Associated with Phospholipids in IDL in Age Bin 3 INTERVAL Samples")

#fastman(bin1, chr="#CHROM", bp="POS", snp="ID", p="P",speedup=TRUE, logp = TRUE, scattermore = FALSE,annotatePval = 0.01, main = "SNPs Associated with Free Cholesterol in Very Large HDL in Male INTERVAL Samples")

dev.off()



#############idlpl bin 4


chr1 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr1.invRankNorm_resids.glm.linear")
chr2 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr2.invRankNorm_resids.glm.linear")
#chr3 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr3.invRankNorm_resids.glm.linear")
#chr4 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr4.invRankNorm_resids.glm.linear")
chr5 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr5.invRankNorm_resids.glm.linear")
#chr6 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr6.invRankNorm_resids.glm.linear")
chr7 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr7.invRankNorm_resids.glm.linear")
chr8 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr8.invRankNorm_resids.glm.linear")
chr9 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr9.invRankNorm_resids.glm.linear")
#chr10 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr10.invRankNorm_resids.glm.linear")
chr11 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr11.invRankNorm_resids.glm.linear")
#chr12 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr12.invRankNorm_resids.glm.linear")
#chr13 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr13.invRankNorm_resids.glm.linear")
#chr14 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr14.invRankNorm_resids.glm.linear")
chr15 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr15.invRankNorm_resids.glm.linear")
#chr16 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr16.invRankNorm_resids.glm.linear")
chr17 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr17.invRankNorm_resids.glm.linear")
#chr18 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr18.invRankNorm_resids.glm.linear")
chr19 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr19.invRankNorm_resids.glm.linear")
chr20 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr20.invRankNorm_resids.glm.linear")
#chr21 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr21.invRankNorm_resids.glm.linear")
#chr22 = fread("/home/cjm260/gwas_files/gwas_results/idlpl3/gwas_chr22.invRankNorm_resids.glm.linear")

#Rbinding them all together
chr_list <- list(chr1, chr2, chr5, chr7, chr8, chr9, chr11, chr15, chr17, chr19, chr20)


# Use do.call with rbind to bind all data frames together
bin4 <- do.call(rbind, chr_list) #260 SNPs. All look significant

#open path to print the plot to
path = "/home/cjm260/gwas_files/gwas_results/idlpl3_ManhatFemPlot.png"
png(path)

manhattan(bin4, chr="#CHROM", bp="POS", snp="ID", p="P", annotatePval = 0.01, main = "SNPs Associated with Phospholipids in IDL in Age Bin 4 INTERVAL Samples")

#fastman(bin1, chr="#CHROM", bp="POS", snp="ID", p="P",speedup=TRUE, logp = TRUE, scattermore = FALSE,annotatePval = 0.01, main = "SNPs Associated with Free Cholesterol in Very Large HDL in Male INTERVAL Samples")
dev.off()

######################################Create a scatter plot of the effect sizes for each SNP
library(ggplot2)

# Merge dataframes by ID
# Merge dataframes by ID
library(tidyr)
library(purrr)
# Renaming columns with suffixes
rename_columns <- function(df, suffix) {
  df %>%
    rename_with(~ paste0(., suffix), c("BETA", "SE", "T_STAT", "P"))
}

# Renaming columns in each data frame
bin1 <- rename_columns(bin1, "_1")
bin2 <- rename_columns(bin2, "_2")
bin3 <- rename_columns(bin3, "_3")
bin4 <- rename_columns(bin4, "_4")

m1 = merge(bin1, bin2, by = "ID")
m2 = merge(m1, bin3, by = "ID")
merged_df = merge(m2, bin4, by = "ID")

# Display the combined data frame
print(merged_df)

pathmerge = "/home/cjm260/gwas_files/gwas_results/idlpl1/fullidlpldata.txt"
write.table(merged_df, pathmerge, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# Extract the BETA values and plot  bin 2 vs. bin 1
x <- merged_df$BETA_1
y <- merged_df$BETA_2
labels <- merged_df$ID
path2 = "/home/cjm260/gwas_files/gwas_results/newNightPlots/idlplScatPlot1.png"
png(path2)
# Create scatter plot
plot(x, y, pch=19, col="black", main="Phospholipids in IDL\nAge Bin 2 vs. Age Bin 1 Effect Size",
     xlab="Age Bin 1 Samples Effect Size (Beta)", ylab="Age Bin 2 Samples Effect Size (Beta)")

# Add labels
text(x, y, labels=labels, pos=4, cex=0.7, col="blue")  # pos=4 places the text to the right of the points
# Add diagonal line y = x
abline(a=0, b=1, col="red", lty=2)  # a=0 is the intercept, b=1 is the slope


dev.off()

# Extract the BETA values and plot  bin 3 vs. bin 2
x <- merged_df$BETA_2
y <- merged_df$BETA_3
labels <- merged_df$ID
path3 = "/home/cjm260/gwas_files/gwas_results/newNightPlots/idlplScatPlot2.png"
png(path3)
# Create scatter plot
plot(x, y, pch=19, col="black", main="Phospholipids in IDL\nAge Bin 3 vs. Age Bin 2 Effect Size",
     xlab="Age Bin 2 Samples Effect Size (Beta)", ylab="Age Bin 3 Samples Effect Size (Beta)")

# Add labels
text(x, y, labels=labels, pos=4, cex=0.7, col="blue")  # pos=4 places the text to the right of the points
# Add diagonal line y = x
abline(a=0, b=1, col="red", lty=2)  # a=0 is the intercept, b=1 is the slope


dev.off()


# Extract the BETA values and plot  bin 4 vs. bin 3
x <- merged_df$BETA_3
y <- merged_df$BETA_4
labels <- merged_df$ID
path4 = "/home/cjm260/gwas_files/gwas_results/newNightPlots/idlplScatPlot3.png"
png(path4)
# Create scatter plot
plot(x, y, pch=19, col="black", main="Phospholipids in IDL\nAge Bin 4 vs. Age Bin 3 Effect Size",
     xlab="Age Bin 3 Samples Effect Size (Beta)", ylab="Age Bin 4 Samples Effect Size (Beta)")

# Add labels
text(x, y, labels=labels, pos=4, cex=0.7, col="blue")  # pos=4 places the text to the right of the points
# Add diagonal line y = x
abline(a=0, b=1, col="red", lty=2)  # a=0 is the intercept, b=1 is the slope


dev.off()



################### Extract age specific SNPs
pathMerged = "/home/cjm260/gwas_files/gwas_results/idlpl1/fullidlpldata.txt"
idlplDf = fread(pathMerged) #nrow 395
nrow(idlplDf)

idlplDf$diff1_2 = abs(idlplDf$BETA_2 - idlplDf$BETA_1)
idlplDf$diff3_2 = abs(idlplDf$BETA_3 - idlplDf$BETA_2)
idlplDf$diff4_3 = abs(idlplDf$BETA_4- idlplDf$BETA_3)
#remove rsIDs where diff > 0.1
idlplDfSub = subset(idlplDf, diff1_2 > 0.1 | diff3_2 > 0.1 | diff4_3 > 0.1) #nrow 93
nrow(idlplDfSub) #93 SNPS to REMOVE

rsIDs = idlplDfSub$ID

path = "/home/cjm260/gwas_files/gwas_results/idlpl1/removeidlplrsIDs.txt"
write.table(rsIDs, path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#load in the Free Cholesterol in Very Large HDL omcispred file and remove the data from removeidlplrsIDs.txt
idlplFile = fread("/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/modelfiles/Metabolomics/OPGS003434_model.txt") #330 SNPs
rsIDsRemove = fread("/home/cjm260/gwas_files/gwas_results/idlpl1/removeidlplrsIDs.txt")

# Step 1: Extract the IDs from rsIDs dataframe
ids_to_remove <- rsIDsRemove$x #18 IDs to remove

# Debugging output: Check for common values
common_ids <- intersect(idlplFile$rsid, ids_to_remove)
print(head(common_ids))
print(length(common_ids)) #93 rsIDs

# Step 2: Filter the idlplFile dataframe to remove rows with these IDs
idlplFile_filtered <- idlplFile[!idlplFile$rsid %in% rsIDsRemove$x, ] #237 SNPs remaining
#idlplFile_filtered contains SNPS which were not noticeably sex-specific. Will now see if taking out the sex specific
#SNPs impacts the Spearman correlation difference between men and women.

path2 = "/home/cjm260/gwas_files/gwas_results/idlpl1/OPGS003434_model_filt.txt"
write.table(idlplFile_filtered, path2, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

########################################Reperform sex stratified analysis
#Part 3 
#Redoing the sex/age specific analysis only for the traits with those differences.

#idlplcine
library(dplyr)
library(data.table)
library(tidyr)

#calculated PGS gile of UKB Samples
df_ukb_pgs = fread("/rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/workingdirectory/idlpl/UKB_Nightingale_idlpl.sscore.gz")

df_ukb_pgs = df_ukb_pgs %>%
  rename(idlpl = score_sum)

#rename columns for the PGS file of all the UKB samples
#colnames(df_ukb_pgs) <- c('IID', sapply(colnames(df_ukb_pgs)[-1], function(col) gsub("\\.pct", "_", strsplit(col, "_")[[1]][2])))
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

row <- df_name_map_overlap[1, ] #test on row 1 (ace, acetate)
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


df_results_1 <- df_results_1 %>% arrange(desc(R2))

#### bin 2
# Initialize results dataframe for BIN 1
df_results_2 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                           R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                           stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)

row <- df_name_map_overlap[1, ] #test on row 1 (ace, acetate)
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


df_results_2 <- df_results_2 %>% arrange(desc(R2))


#### bin 3
# Initialize results dataframe for BIN 1
df_results_3 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                           R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                           stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)

row <- df_name_map_overlap[1, ] #test on row 1 (ace, acetate)
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


df_results_3 <- df_results_3 %>% arrange(desc(R2))


#### bin 4
# Initialize results dataframe for BIN 1
df_results_4 <- data.frame(INTERVAL_Name = character(), UKB_Name = character(), R2 = numeric(),
                           R2_P = numeric(), SR = numeric(), SR_P = numeric(), N_val_samples = integer(),
                           stringsAsFactors = FALSE)

ids <- c()
#nrow(df_name_map_overlap)

row <- df_name_map_overlap[1, ] #test on row 1 (ace, acetate)
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


df_results_4 <- df_results_4 %>% arrange(desc(R2))


#combine age bin results, one on top of other
fulldf = rbind(df_results_1, df_results_2, df_results_3, df_results_4)
fulldf$Age_Bin = c("Bin 1", "Bin 2", "Bin 3", "Bin 4")
fulldf = fulldf %>%
  select(Age_Bin, everything())
pathfulldf = "/home/cjm260/pt3_results/NightAgeSpec/idlplFullResults.txt"
write.table(fulldf, pathfulldf, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

#DF of just SLOPE FOR PAPER
# Assuming df_results_male is already defined and slopes is a vector of the same length
# Create the new dataframe
df_new = df_results_1 %>%
  select(INTERVAL_Name, UKB_Name)

df_new$Slope_B1_B2 = abs(df_results_2$SR - df_results_1$SR)
df_new$Slope_B2_B3 = abs(df_results_3$SR - df_results_2$SR)
df_new$Slope_B3_B4 = abs(df_results_4$SR - df_results_3$SR)
df_new$Avg_Slope = (df_new$Slope_B1_B2 + df_new$Slope_B2_B3+df_new$Slope_B3_B4)/3
# Display the new dataframe
print(df_new)
pathnewdf = "/home/cjm260/pt3_results/NightAgeSpec/idlplSlope.txt"
write.table(df_new, pathnewdf, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)





