#!/usr/bin/env Rscript
install.packages("writexl")
install.packages("openxlsx")
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(writexl)
  library(openxlsx)
  
})

# -------------------------
# Input / Output
# -------------------------
in_xlsx  <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/sv_analysis/SV_calling/result/6_IGV/all_proband_filtered_1206_remove_dup.xlsx"
out_xlsx <- sub("\\.xlsx$", "_parsed.xlsx", in_xlsx)

# ---------- INFO keys ----------
info_keys_ordered <- c(
  "Genotype","Alt_reads","Ref_reads","Total_reads",
  "Sample_support","Allele_Freq_ALL",
  "Pop_Count_AFR","Pop_Freq_AFR",
  "Pop_Count_AMR","Pop_Freq_AMR",
  "Pop_Count_EAS","Pop_Freq_EAS",
  "Pop_Count_EUR","Pop_Freq_EUR",
  "Pop_Count_SAS","Pop_Freq_SAS",
  "Pop_Count_ALL","Pop_Freq_ALL",
  "Allele_Count_AFR","Allele_Freq_AFR",
  "Allele_Count_AMR","Allele_Freq_AMR",
  "Allele_Count_EAS","Allele_Freq_EAS",
  "Allele_Count_EUR","Allele_Freq_EUR",
  "Allele_Count_SAS","Allele_Freq_SAS",
  "Allele_Count_ALL","Allele_Freq_ALL2"
)

# ---------- Core columns ----------
core_cols <- c(
  "Samples_ID",
  "SV_chrom","SV_start","SV_end","SV_type","SV_length",
  "Gene_name","Closest_left","Closest_right","Gene_count","RE_gene",
  "clinical_significance","clinical_source",
  "AnnotSV_ranking_score","AnnotSV_ranking_criteria","ACMG_class"
)

# ---------- Helper ----------
extract_info_value <- function(info_string, key) {
  if (is.na(info_string) || !nzchar(info_string)) return(NA_character_)
  m <- str_match(info_string, paste0("(?:(?:^)|;)", key, "=([^;]*)"))
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

# ---------- Read ----------
df <- read_excel(in_xlsx)

if (!("INFO" %in% names(df))) {
  stop("ERROR: Input file does not contain an 'INFO' column.")
}

# ---------- Parse INFO ----------
parsed_cols <- map_dfc(
  info_keys_ordered,
  ~ tibble(!!.x := vapply(df$INFO, extract_info_value,
                          FUN.VALUE = character(1), key = .x))
)

df2 <- bind_cols(df, parsed_cols)

# ---------- Convert numerics ----------
numeric_keys <- setdiff(info_keys_ordered, c("Sample_support","Genotype"))
df2 <- df2 %>%
  mutate(across(all_of(numeric_keys), ~ suppressWarnings(as.numeric(.x))))

# ---------- Rename + recode ACMG ----------
df2 <- df2 %>%
  rename(ACMG_class_from_AnnotSV = ACMG_class) %>%
  mutate(
    ACMG_class_from_AnnotSV = case_when(
      ACMG_class_from_AnnotSV == 1 ~ "Benign",
      ACMG_class_from_AnnotSV == 2 ~ "Likely_Benign",
      ACMG_class_from_AnnotSV == 3 ~ "VUS",
      ACMG_class_from_AnnotSV == 4 ~ "Likely_Pathogenic",
      ACMG_class_from_AnnotSV == 5 ~ "Pathogenic",
      TRUE ~ as.character(ACMG_class_from_AnnotSV)
    )
  )

core_cols <- sub("^ACMG_class$", "ACMG_class_from_AnnotSV", core_cols)

# ---------- Final column order ----------
re_gene_idx <- match("RE_gene", core_cols)
if (is.na(re_gene_idx)) stop("ERROR: 'RE_gene' not found in core_cols.")

final_cols <- c(
  core_cols[1:re_gene_idx],
  info_keys_ordered,
  core_cols[(re_gene_idx + 1):length(core_cols)]
)

df_final <- df2 %>%
  select(any_of(final_cols)) %>%
  
  # ---------- ROW SORT ----------
arrange(
  Samples_ID,
  SV_chrom,
  SV_start,
  SV_end
)

# ---------- Write ----------
openxlsx::write.xlsx(df_final, out_xlsx, overwrite = TRUE)

message("Wrote compact supplementary table: ", out_xlsx)
