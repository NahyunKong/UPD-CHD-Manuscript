install.packages("data.table")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
library(data.table)
library(dplyr)
suppressMessages(library(tidyverse))
suppressMessages(library(grid))
suppressMessages(library(VariantAnnotation))
if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr")
}
install.packages("writexl")
library(writexl)
library(vcfR)
library(plyr)
library(tidyr)  
library(dplyr)  
##sample info
# Define the base path
base_path <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/final"

# Create the samples data frame with sample IDs, parent, and chromosome regions
samples <- data.frame(
  sample_id = c("1-14785", "1-14523", "1-06216", "1-08240", "1-03820", "1-15162", "1-03573", "1-15783","1-00180","1-02917","1-04975"),
  parent = c("M", "M", "M", "M", "M", "F", "M", "M", "M", "M", "M"),
  chromosome_regions = c(
    "chr16:15561-1656434, chr16:61655365-82858291, chr16:86541212-900605855",
    "chr16:254515-13920136, chr16:84762997-90095778",
    "chr16:227459-8746241, chr16:57728489-89765102",
    "chr4:85671-7717984, chr4:25003349-55890365, chr4:173333596-177339880",
    "chr15:75337984-96757410",
    "chr19:43204093-43259192",
    "chr8:473541-23566156",
    "chr12:14838-133219299",
    "chr14:16133531-20617678",
    "chr16:35629-5300048, 75909788-90032062",
    "chr9:52745-1197232"
  ),
  genome_build = "hg38"
)

# Generate the VCF file paths dynamically based on the sample IDs
samples$txt_file <- paste0(
  base_path, "/", 
  #samples$sample_id, "/", 
  samples$sample_id, 
  "_annot.hg38_multianno.txt"
)
samples$vcf_file <- paste0(
  base_path, "/", 
  #samples$sample_id, "/", 
  samples$sample_id, 
  "_annot.hg38_multianno.vcf"
)
# Function to parse chromosome regions
parse_regions <- function(region_str) {
  regions_split <- strsplit(region_str, ",")[[1]]
  region_df <- do.call(rbind, lapply(regions_split, function(x) {
    chr_region <- strsplit(x, ":|-")[[1]]
    data.frame(
      chromosome = chr_region[1],
      start = as.integer(chr_region[2]),
      end = as.integer(chr_region[3])
    )
  }))
  return(region_df)
}

# Function to read and filter a VCF file based on chromosome, regions, and inheritance pattern
filter_vcf <- function(txt_file, vcf_file, chromosome, regions, inheritance_pattern, sample_id) {
  # Read the txt file using vcfR
  #txt_file="/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/final/1-14785_annot.hg38_multianno.txt"
  txt_data <- fread(txt_file, sep = "\t", header = TRUE)
  # info fields
  info_cols <- txt_data %>% dplyr::rename(ID = Otherinfo6) %>%
    dplyr::select(-c(Chr, Start, End, Ref, Alt)) %>%
    dplyr::select(-starts_with("Otherinfo"))
  
  # Separate Otherinfo13 (Proband), Otherinfo14 (Mother), and Otherinfo15 (Father) into GT, DP, AD, and GQ
  new_data <- txt_data %>%
    tidyr::separate(Otherinfo13, into = c("GT_Proband",  "AD_Proband", "DP_Proband","GQ_Proband", "PL_Proband"), sep = ":", fill = "right") %>%
    dplyr::select(-PL_Proband) %>%  # Keep only relevant columns
    tidyr::separate(Otherinfo14, into = c("GT_Mother", "AD_Mother", "DP_Mother",  "GQ_Mother", "PL_Mother"), sep = ":", fill = "right") %>%
    dplyr::select(-PL_Mother) %>%
    tidyr::separate(Otherinfo15, into = c("GT_Father", "AD_Father", "DP_Father", "GQ_Father", "PL_Father"), sep = ":", fill = "right") %>%
    dplyr::select( -PL_Father)
  
  
  # Extract fixed fields (Chr, Start, End, Ref, Alt) and separate INFO field
  fixed_data <- new_data %>%
    dplyr::select(Chr, Start, End, Ref, Alt)
  
  # Combine fixed fields with genotype data (GT, DP, AD, GQ) and INFO
  combined_data <- cbind(fixed_data, new_data[, c("GT_Proband", "DP_Proband", "AD_Proband", "GQ_Proband", "GT_Mother", "DP_Mother", "AD_Mother", "GQ_Mother", "GT_Father", "DP_Father", "AD_Father", "GQ_Father")], info_cols )
  
  
  # Filter based on the chromosome and regions
  filtered_df <- combined_data %>%
    rowwise() %>%  # Process each row individually
    filter(
      Chr == chromosome,  # Ensure chromosome matches
      any(map_lgl(1:nrow(regions), ~ (Start >= regions[.x, "start"] & Start <= regions[.x, "end"])))
    ) 
  # parents and proband DP >= 8 GQ>=20 
  filtered_df <- filtered_df %>%
    filter(!is.na(suppressWarnings(as.numeric(DP_Proband))) & suppressWarnings(as.numeric(DP_Proband)) >= 8) %>%
    filter(!is.na(suppressWarnings(as.numeric(DP_Mother))) & suppressWarnings(as.numeric(DP_Mother)) >= 4) %>%
    filter(!is.na(suppressWarnings(as.numeric(DP_Father))) & suppressWarnings(as.numeric(DP_Father)) >= 4)
  
  filtered_df <- filtered_df %>%
    filter(!is.na(suppressWarnings(as.numeric(GQ_Proband))) & suppressWarnings(as.numeric(GQ_Proband)) >= 20) %>%
    filter(!is.na(suppressWarnings(as.numeric(GQ_Mother))) & suppressWarnings(as.numeric(GQ_Mother)) >= 10) %>%
    filter(!is.na(suppressWarnings(as.numeric(GQ_Father))) & suppressWarnings(as.numeric(GQ_Father)) >= 10)

  # Apply the inheritance pattern filtering
  if (inheritance_pattern == "M") {
    recessive <- filtered_df %>%
      filter(!GT_Father %in% c("./1", "0/1", "1/2", "1/3","1/4","1/5", "1/1", ".|1", "0|1", "1|2", "1|3","1|4","1|5", "1|1")) %>%
      filter(GT_Mother %in% c("./1", "0/1", "1/2", "1/3", "1/4","1/5",".|1", "0|1", "1|2", "1|3", "1|4","1|5")) %>%
      filter(GT_Proband == "1/1"|GT_Proband == "1|1")
  } else if (inheritance_pattern == "F") {
    recessive <- filtered_df %>%
      filter(!GT_Mother %in% c("./1", "0/1", "1/2", "1/3","1/4","1/5", "1/1", ".|1", "0|1", "1|2", "1|3","1|4","1|5", "1|1")) %>%
      filter(GT_Father %in% c("./1", "0/1", "1/2", "1/3", "1/4","1/5",".|1", "0|1", "1|2", "1|3", "1|4","1|5")) %>%
      filter(GT_Proband == "1/1"|GT_Proband == "1|1")
  } else {
    stop("Invalid inheritance pattern. Use 'M' for maternal or 'P' for paternal.")
  }
  print(head(recessive))
  # Extract recessive IDs
  recessive_ids <- recessive$ID
  
  # Read the VCF file using vcfR
  vcf <- read.vcfR(vcf_file)
  
  # Extract the ID column from the VCF
  vcf_data <- getFIX(vcf)
  vcf_ids <- vcf_data[, "ID"]
  
  # Filter VCF by matching IDs
  matching_ids <- vcf_ids %in% recessive_ids
  print(head(vcf_ids))
  print(head(recessive_ids))
  filtered_vcf <- vcf[matching_ids, ]
  
  # Define the output VCF file path
  output_vcf <- paste0(base_path, "/",  sample_id, "_isodisomy_recessive_dp_mq.vcf.gz")
  
  # Write the filtered VCF back to a file
  write.vcf(filtered_vcf, file = output_vcf)
  
  return(recessive)
}

# Function to transform the combined VCF data frame
transform_vcf_data <- function(txt_df, sample_id) {
  # Create the transformed data frame by combining the original columns with the new columns
  transformed_df <- cbind(
    txt_df[, c(1:17)],
    ProbandID = sample_id, 
    FatherID = paste(sample_id, "02", sep = "-"),
    MotherID = paste(sample_id, "01", sep = "-"),
    txt_df[, c(18:ncol(txt_df))]
  )
  
  return(transformed_df)
}
# Initialize an empty list to store transformed data frames
transformed_dfs <- data.frame()

# Loop over each sample
for (i in 1:nrow(samples)) {
  # Get sample details
  sample_id <- samples$sample_id[i]
#  if (sample_id %in% c("1-02917")) {
  txt_file <- samples$txt_file[i]
  genome_build <- samples$genome_build[i]
  chromosome_regions <- samples$chromosome_regions[i]
  parent <- samples$parent[i]
  vcf_file <- samples$vcf_file[i]
  # Parse the regions
  regions <- parse_regions(chromosome_regions)
  
  # Choose the inheritance pattern based on parent
  inheritance_pattern <- parent
  # Filter the VCF # Load necessary libraries
  print(c(sample_id,regions$chromosome[1], regions, inheritance_pattern,txt_file) )
  # Filter the VCF & write vcf
  filtered_txt <- filter_vcf(txt_file, vcf_file, regions$chromosome[1], regions, parent, sample_id)
  #print("filtered_txt")
  #head(filtered_txt )
  
  # Transform the combined VCF data frame
  transformed_df <- transform_vcf_data(filtered_txt, sample_id)
  #print(transformed_df)  
  # Store the transformed data frame in the list
  transformed_dfs <- rbind(transformed_dfs, transformed_df)
#  }
}


# View the combined data frame
head(transformed_dfs)

# Step 4: Rare damaging filtering
new_data <- transformed_dfs
new_data$gnomAD_exome_ALL[which(new_data$gnomAD_exome_ALL == "." )] <- 0
new_data$gnomAD_genome_ALL[which(new_data$gnomAD_genome_ALL == "." )] <- 0
new_data$bravo[which(new_data$bravo == ".")] <- 0
#new_data$CADD_phred[which(new_data$CADD_phred == ".")] <- 0


# Filter for rare damaging variants without gnomAD genome
damaging_RV <- new_data %>%
  filter(as.numeric(gnomAD_exome_ALL) < 0.1) %>%
  filter(as.numeric(bravo) < 0.1) %>%
  filter(as.numeric(gnomAD_genome_ALL) < 0.1) %>%
  filter((is.na(CADD_phred) | suppressWarnings((as.numeric(CADD_phred)) >= 10 )& ExonicFunc.refGene == "nonsynonymous SNV") |
           (ExonicFunc.refGene == "stopgain" |
              ExonicFunc.refGene == "frameshift insertion" |
              ExonicFunc.refGene == "frameshift deletion" | 
              Func.refGene == "splicing"))
# Step 5: Write the final data to a new file
file_name <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/recessive_calling/hg38_damaging_RV_cadd10.txt"
fwrite(damaging_RV , file_name, sep = "\t")

