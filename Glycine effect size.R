library(dplyr)
library(data.table)

merged_data = fread("/home/cjm260/gwas_files/gwas_results/Gly/fullGlydata.txt")
x <- merged_data$BETA_male
y <- merged_data$BETA_female
labels <- merged_data$ID
path2 = "/home/cjm260/gwas_files/gwas_results/newNightPlots/glyScatPlotMF.png"



png(path2)
# Create scatter plot
plot(x, y, pch=19, col="black", main="Glycine\nFemale vs. Male SNP Effect Size",
     xlab="Male Samples Effect Size (Beta)", ylab="Female Samples Effect Size (Beta)")

# Add labels
text(x, y, labels=labels, pos=4, cex=0.7, col="blue")  # pos=4 places the text to the right of the points
# Add diagonal line y = x
abline(a=0, b=1, col="red", lty=2)  # a=0 is the intercept, b=1 is the slope


dev.off()

#####idlp
library(dplyr)
library(data.table)

merged_data = fread("/home/cjm260/gwas_files/gwas_results/idlp/fullidlpdata.txt")
x <- merged_data$BETA_male
y <- merged_data$BETA_female
labels <- merged_data$ID
path2 = "/home/cjm260/gwas_files/gwas_results/newNightPlots/idlpScatPlotMF.png"



png(path2)
# Create scatter plot
plot(x, y, pch=19, col="black", main="Concentration of IDL Particles\nFemale vs. Male SNP Effect Size",
     xlab="Male Samples Effect Size (Beta)", ylab="Female Samples Effect Size (Beta)")

# Add labels
text(x, y, labels=labels, pos=4, cex=0.7, col="blue")  # pos=4 places the text to the right of the points
# Add diagonal line y = x
abline(a=0, b=1, col="red", lty=2)  # a=0 is the intercept, b=1 is the slope


dev.off()



#### mldlfc
library(dplyr)
library(data.table)

merged_data = fread("/home/cjm260/gwas_files/gwas_results/mldlfc/fullmldlfcdata.txt")
x <- merged_data$BETA_male
y <- merged_data$BETA_female
labels <- merged_data$ID
path2 = "/home/cjm260/gwas_files/gwas_results/newNightPlots/mldlfcScatPlotMF.png"



png(path2)
# Create scatter plot
plot(x, y, pch=19, col="black", main="Free cholesterol in medium LDL\nFemale vs. Male SNP Effect Size",
     xlab="Male Samples Effect Size (Beta)", ylab="Female Samples Effect Size (Beta)")

# Add labels
text(x, y, labels=labels, pos=4, cex=0.7, col="blue")  # pos=4 places the text to the right of the points
# Add diagonal line y = x
abline(a=0, b=1, col="red", lty=2)  # a=0 is the intercept, b=1 is the slope


dev.off()



#mvldlce
library(dplyr)
library(data.table)

merged_data = fread("/home/cjm260/gwas_files/gwas_results/mvldlce/fullmvldlcedata.txt")
x <- merged_data$BETA_male
y <- merged_data$BETA_female
labels <- merged_data$ID
path2 = "/home/cjm260/gwas_files/gwas_results/newNightPlots/mvldlceScatPlotMF.png"



png(path2)
# Create scatter plot
plot(x, y, pch=19, col="black", main="Cholesteryl esters in medium VLDL\nFemale vs. Male SNP Effect Size",
     xlab="Male Samples Effect Size (Beta)", ylab="Female Samples Effect Size (Beta)")

# Add labels
text(x, y, labels=labels, pos=4, cex=0.7, col="blue")  # pos=4 places the text to the right of the points
# Add diagonal line y = x
abline(a=0, b=1, col="red", lty=2)  # a=0 is the intercept, b=1 is the slope


dev.off()
