#Subsetting rsIDs in .pgen files to only 

#code successful in extracting rsIDs
# Initialize an output file
output_file="/home/cjm260/gwas_files/trait_rsIDs/gly_rsIDs.txt"
> $output_file  # This creates an empty file or clears the existing file

# Extract the rsID column from 
#glycine is OPGS003445
awk '{print $1}' /rds/project/jmmh2/rds-jmmh2-projects/inouye_lab_other/impute_genomics/omics_age_sex/modelfiles/Metabolomics/OPGS003445_model.txt >> $output_file


#Subset and merge the .pgen files
#!/bin/bash

# Step 1: Set up paths and filenames
rsids_file="/home/cjm260/gwas_files/trait_rsIDs/gly_rsIDs.txt"
output_dir="/home/cjm260/gwas_files/filtINTERVALdata/gly"
filtered_files_list="${output_dir}/filtered_files.txt"
merged_output_prefix="${output_dir}/merged_filtered"


# Step 3: Initialize the filtered files list
> $filtered_files_list

# Step 4: Loop through each chromosome and subset the .pgen files
for chr in {1..22}
do
    input_prefix="/rds/project/rds-pNR2rM6BWWA/interval/imputed/uk10k_1000g_b37/imputed/plink_format/pgen/impute_dedup_${chr}_interval"  # Change this to the actual path of your input files
    output_prefix="${output_dir}/chr${chr}_filtered"
    
    # Subset the .pgen, .pvar, and .psam files
    plink2 --pfile $input_prefix --extract $rsids_file --make-pgen --out $output_prefix
    
    # Append the output prefix to the filtered files list for merging
    echo "${output_prefix}" >> $filtered_files_list
done

# Step 5: Merge the filtered .pgen files
plink2 --pmerge-list $filtered_files_list --make-pgen --out $merged_output_prefix

# Step 6: Clean up
rm $filtered_files_list

