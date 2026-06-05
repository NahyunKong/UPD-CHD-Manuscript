suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(forestplot)
  library(grid)
})

vcf_wgs <- "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/all_chr_for_upd_and_methly_gene_analysis/burden_test/RESULT/pcgc-wgs-n6430.UPDregions.509trio_annotated.hg38_multianno.vcf"
vcf_exome <- "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/all_chr_for_upd_and_methly_gene_analysis/burden_test/data/fromKen/PCGC-exomes-UPDregions-3869trio_annotated.hg38_multianno.vcf"
vcf_joint_106216 <- "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/all_chr_for_upd_and_methly_gene_analysis/results/joint_call/1-06216/1-06216_annot.hg38_multianno.vcf"
trio_file <- "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/all_chr_for_upd_and_methly_gene_analysis/burden_test/data/burden_test_trio_blindID_final3739.txt"
outdir <- "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/all_chr_for_upd_and_methly_gene_analysis/burden_test/RESULT"

filtered_cache_file <- file.path(outdir, "all_variants_ALT_DP8_cadd_aachange.tsv")

case_ids <- c(
  "1-00660","1-08240","1-03573","1-04975",
  "1-15783","1-00180","1-03820","1-02917","1-06216",
  "1-14523","1-14785","1-15162"
)

af_threshold <- 0.001
dp_min <- 8

manual_groups <- tribble(
  ~gene, ~pathway_group,
  "TEX12", "Synaptonemal complex",
  "C14orf39", "Synaptonemal complex",
  "SYCE1", "Synaptonemal complex",
  "SYCE2", "Synaptonemal complex",
  "SYCE3", "Synaptonemal complex",
  "SYCP2", "Synaptonemal complex",
  "SYCP3", "Synaptonemal complex",
  "SYCP1", "Synaptonemal complex",
  "SPO11", "Meiotic recombination",
  "DMC1", "Meiotic recombination",
  "RAD51", "Meiotic recombination",
  "RAD51C", "Meiotic recombination",
  "MSH4", "Meiotic recombination",
  "MSH5", "Meiotic recombination",
  "MLH1", "Meiotic recombination",
  "MLH3", "Meiotic recombination",
  "MUS81", "Meiotic recombination",
  "SMC1B", "Cohesin complex",
  "SMC1A", "Cohesin complex",
  "SMC3", "Cohesin complex",
  "REC8", "Cohesin complex",
  "RAD21", "Cohesin complex",
  "RAD21L", "Cohesin complex",
  "STAG3", "Cohesin complex",
  "STAG1", "Cohesin complex",
  "STAG2", "Cohesin complex",
  "WAPL", "Cohesin complex",
  "PDS5A", "Cohesin complex",
  "PDS5B", "Cohesin complex",
  "NIPBL", "Cohesin complex",
  "MAU2", "Cohesin complex",
  "SGO1", "Cohesin complex",
  "SGO2", "Cohesin complex",
  "BUB1", "Spindle assembly checkpoint",
  "BUB1B", "Spindle assembly checkpoint",
  "BUB3", "Spindle assembly checkpoint",
  "MAD1L1", "Spindle assembly checkpoint",
  "MAD2L1", "Spindle assembly checkpoint",
  "CDC20", "Spindle assembly checkpoint",
  "PLK4", "Centrosome biogenesis",
  "STIL", "Centrosome biogenesis",
  "CEP152", "Centrosome biogenesis",
  "SAS6", "Centrosome biogenesis",
  "CENPJ", "Centrosome biogenesis",
  "WDR90", "Centrosome biogenesis",
  "DYNC1H1", "Spindle positioning",
  "NUMA1", "Spindle positioning",
  "KIF11", "Spindle positioning",
  "KIF2C", "Spindle positioning"
)

target_genes <- unique(manual_groups$gene)

get_vcf_colnames <- function(vcf_file) {
  header_line <- system(paste("grep -m 1 '^#CHROM' ", shQuote(vcf_file)), intern = TRUE)
  strsplit(header_line, "\t")[[1]]
}

is_alt_gt <- function(gt) {
  !is.na(gt) &
    !gt %in% c("0/0", "0|0", "./.", ".|.") &
    grepl("[1-9]", gt)
}

classify_variant <- function(func, exonic_func) {
  case_when(
    func == "splicing" ~ "LoF",
    exonic_func %in% c(
      "stopgain", "stoploss", "startloss",
      "frameshift_insertion",
      "frameshift_deletion",
      "frameshift_substitution"
    ) ~ "LoF",
    exonic_func == "nonsynonymous_SNV" ~ "Missense",
    exonic_func == "synonymous_SNV" ~ "Synonymous",
    TRUE ~ "Other"
  )
}

parse_info <- function(info, key) {
  pattern <- paste0("(^|;)", key, "=([^;]*)")
  out <- str_match(info, pattern)[, 3]
  out[out %in% c(".", "", "NA")] <- NA
  out
}

parse_info_num <- function(info, key) {
  suppressWarnings(as.numeric(parse_info(info, key)))
}

parse_sample_fast <- function(fmt, sample_str) {
  sample_parts <- tstrsplit(sample_str, ":", fixed = TRUE)
  fmt_parts <- strsplit(fmt, ":", fixed = TRUE)
  
  get_field <- function(field) {
    idx <- vapply(fmt_parts, function(x) match(field, x), integer(1))
    out <- rep(NA_character_, length(sample_str))
    
    for (k in unique(idx[!is.na(idx)])) {
      rows <- which(idx == k)
      if (k <= length(sample_parts)) {
        out[rows] <- sample_parts[[k]][rows]
      }
    }
    
    out[out %in% c(".", "", "NA")] <- NA
    out
  }
  
  data.table(
    GT = get_field("GT"),
    DP = suppressWarnings(as.numeric(get_field("DP")))
  )
}

run_fisher <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2)
  ft <- fisher.test(mat)
  
  mat_or <- mat
  if (any(mat_or == 0)) mat_or <- mat_or + 0.5
  
  a2 <- mat_or[1, 1]
  b2 <- mat_or[2, 1]
  c2 <- mat_or[1, 2]
  d2 <- mat_or[2, 2]
  
  or <- (a2 * d2) / (b2 * c2)
  se <- sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2)
  
  data.frame(
    odds_ratio = or,
    ci_low = exp(log(or) - 1.96 * se),
    ci_high = exp(log(or) + 1.96 * se),
    p_value = ft$p.value
  )
}

trio_map <- fread(trio_file, header = FALSE, select = 1:3)
colnames(trio_map) <- c("child_id", "mother_id", "father_id")

sample_info <- bind_rows(
  trio_map %>% transmute(sample_id = child_id, family_id = child_id, role = "child", type = "child"),
  trio_map %>% transmute(sample_id = mother_id, family_id = child_id, role = "mother", type = "parent"),
  trio_map %>% transmute(sample_id = father_id, family_id = child_id, role = "father", type = "parent")
) %>%
  mutate(group = ifelse(family_id %in% case_ids, "case", "control"))

all_trio_samples <- unique(sample_info$sample_id)

read_vcf_samples <- function(vcf_file, sample_ids, source_name) {
  
  message("Reading header for ", source_name, "...")
  vcf_colnames <- get_vcf_colnames(vcf_file)
  available_samples <- intersect(sample_ids, vcf_colnames)
  
  message(source_name, " samples requested: ", length(sample_ids))
  message(source_name, " samples available: ", length(available_samples))
  
  if (length(available_samples) == 0) return(data.table())
  
  message("Reading VCF body for ", source_name, "...")
  
  dt <- fread(
    cmd = paste("grep -v '^#' ", shQuote(vcf_file)),
    header = FALSE,
    sep = "\t",
    quote = ""
  )
  
  setnames(dt, vcf_colnames)
  
  dt[, Func.refGene := parse_info(INFO, "Func.refGene")]
  dt[, Gene.refGene := parse_info(INFO, "Gene.refGene")]
  dt[, ExonicFunc.refGene := parse_info(INFO, "ExonicFunc.refGene")]
  dt[, AAChange.refGene := parse_info(INFO, "AAChange.refGene")]
  
  dt[, Gene.refGene := str_replace_all(Gene.refGene, "\\\\x3b", ";")]
  dt[, AAChange.refGene := str_replace_all(AAChange.refGene, "\\\\x3b", ";")]
  
  dt[, gnomAD_exome_ALL := parse_info_num(INFO, "gnomAD_exome_ALL")]
  dt[, gnomAD_genome_ALL := parse_info_num(INFO, "gnomAD_genome_ALL")]
  dt[, bravo := parse_info_num(INFO, "bravo_freeze8")]
  dt[, CADD_raw := parse_info_num(INFO, "CADD_raw")]
  dt[, CADD_phred := parse_info_num(INFO, "CADD_phred")]
  
  fmt <- dt$FORMAT
  
  out_list <- lapply(available_samples, function(sid) {
    parsed <- parse_sample_fast(fmt, dt[[sid]])
    
    x <- data.table(
      source = source_name,
      sample_id = sid,
      Chr = dt$`#CHROM`,
      Start = dt$POS,
      Ref = dt$REF,
      Alt = dt$ALT,
      Gene.refGene = dt$Gene.refGene,
      Func.refGene = dt$Func.refGene,
      ExonicFunc.refGene = dt$ExonicFunc.refGene,
      AAChange.refGene = dt$AAChange.refGene,
      GT = parsed$GT,
      DP = parsed$DP,
      gnomAD_exome_ALL = dt$gnomAD_exome_ALL,
      gnomAD_genome_ALL = dt$gnomAD_genome_ALL,
      bravo = dt$bravo,
      CADD_phred = dt$CADD_phred
    )
    
    x[
      is_alt_gt(GT) &
        !is.na(DP) &
        DP >= dp_min
    ]
  })
  
  out <- rbindlist(out_list, fill = TRUE)
  
  out <- out %>%
    mutate(
      variant_class = classify_variant(Func.refGene, ExonicFunc.refGene),
      max_pop_af = pmax(
        gnomAD_exome_ALL,
        gnomAD_genome_ALL,
        bravo,
        na.rm = TRUE
      ),
      max_pop_af = ifelse(is.infinite(max_pop_af), NA, max_pop_af)
    )
  
  message(source_name, " ALT + DP>=", dp_min, " variants kept: ", nrow(out))
  out
}

wgs_colnames <- get_vcf_colnames(vcf_wgs)
exome_colnames <- get_vcf_colnames(vcf_exome)
joint_colnames <- get_vcf_colnames(vcf_joint_106216)

wgs_samples <- intersect(all_trio_samples, wgs_colnames)
remaining_after_wgs <- setdiff(all_trio_samples, wgs_samples)
exome_fallback_samples <- intersect(remaining_after_wgs, exome_colnames)
remaining_after_exome <- setdiff(remaining_after_wgs, exome_fallback_samples)
joint_rescue_samples <- intersect(remaining_after_exome, joint_colnames)
still_missing <- setdiff(remaining_after_exome, joint_rescue_samples)

message("Total trios: ", length(all_trio_samples) / 3)
message("Total samples: ", length(all_trio_samples))
message("Samples loaded from WGS: ", length(wgs_samples))
message("Samples loaded from exome fallback: ", length(exome_fallback_samples))
message("Samples rescued from joint-call VCF: ", length(joint_rescue_samples))
message("Samples still missing: ", length(still_missing))

if (length(still_missing) > 0) {
  message("Missing sample IDs:")
  print(still_missing)
}

if (file.exists(filtered_cache_file) && file.info(filtered_cache_file)$size > 0) {
  message("Loading filtered cached variants from: ", filtered_cache_file)
  all_variants_raw <- fread(filtered_cache_file)
  
  if (!"AAChange.refGene" %in% colnames(all_variants_raw)) {
    stop("Cached file does not contain AAChange.refGene. Delete cache or use updated filtered_cache_file.")
  }
  
} else {
  message("No filtered cache found. Reading VCFs from scratch...")
  
  wgs_variants <- read_vcf_samples(vcf_wgs, wgs_samples, "WGS")
  exome_variants <- read_vcf_samples(vcf_exome, exome_fallback_samples, "EXOME_fallback")
  joint_variants <- read_vcf_samples(vcf_joint_106216, joint_rescue_samples, "JOINT_RESCUE")
  
  all_variants_raw <- rbindlist(
    list(wgs_variants, exome_variants, joint_variants),
    fill = TRUE
  )
  
  fwrite(all_variants_raw, file = paste0(filtered_cache_file, ".tmp"), sep = "\t")
  file.rename(paste0(filtered_cache_file, ".tmp"), filtered_cache_file)
}

all_variants_raw <- all_variants_raw %>%
  mutate(
    Gene.refGene = str_replace_all(Gene.refGene, "\\\\x3b", ";"),
    AAChange.refGene = str_replace_all(AAChange.refGene, "\\\\x3b", ";")
  )

all_variants <- all_variants_raw %>%
  inner_join(sample_info, by = "sample_id")

child_df <- all_variants %>% filter(type == "child")
parent_df <- all_variants %>% filter(type == "parent")

sample_info_vcf <- sample_info %>%
  filter(sample_id %in% c(wgs_samples, exome_fallback_samples, joint_rescue_samples))

n_case_child <- sample_info_vcf %>% filter(type == "child", group == "case") %>% distinct(sample_id) %>% nrow()
n_control_child <- sample_info_vcf %>% filter(type == "child", group == "control") %>% distinct(sample_id) %>% nrow()
n_case_parent <- sample_info_vcf %>% filter(type == "parent", group == "case") %>% distinct(sample_id) %>% nrow()
n_control_parent <- sample_info_vcf %>% filter(type == "parent", group == "control") %>% distinct(sample_id) %>% nrow()

run_burden <- function(df, n_case, n_ctrl) {
  
  rare <- df %>% filter(is.na(max_pop_af) | max_pop_af < af_threshold)
  
  rare_burden <- rare %>%
    filter(variant_class %in% c("LoF", "Missense")) %>%
    mutate(class = "Rare LoF or Missense")
  
  gene_df <- rare_burden %>%
    mutate(Gene.refGene = str_replace_all(Gene.refGene, "\\\\x3b", ";")) %>%
    separate_rows(Gene.refGene, sep = ";|,") %>%
    mutate(Gene.refGene = str_trim(Gene.refGene)) %>%
    filter(!is.na(Gene.refGene), Gene.refGene != ".", Gene.refGene != "") %>%
    rename(gene = Gene.refGene) %>%
    filter(gene %in% target_genes)
  
  carrier <- gene_df %>%
    distinct(sample_id, group, gene, class) %>%
    count(gene, class, group, name = "n")
  
  wide <- carrier %>%
    pivot_wider(names_from = group, values_from = n, values_fill = 0)
  
  if (!"case" %in% colnames(wide)) wide$case <- 0
  if (!"control" %in% colnames(wide)) wide$control <- 0
  
  wide <- wide %>%
    mutate(
      case_carriers = case,
      control_carriers = control,
      case_noncarriers = n_case - case,
      control_noncarriers = n_ctrl - control
    )
  
  wide %>%
    rowwise() %>%
    mutate(
      fisher = list(run_fisher(case_carriers, case_noncarriers, control_carriers, control_noncarriers))
    ) %>%
    unnest(fisher) %>%
    ungroup() %>%
    group_by(class) %>%
    mutate(p_adj_bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
    ungroup() %>%
    left_join(manual_groups, by = "gene") %>%
    arrange(p_value)
}

child_res <- run_burden(child_df, n_case_child, n_control_child)
parent_res <- run_burden(parent_df, n_case_parent, n_control_parent)

make_forestplot <- function(burden_results_forest,
                            n_case_plot,
                            n_control_plot,
                            plot_title,
                            top_n = 15) {
  
  plot_min <- 0.05
  plot_max <- 20
  
  forest_df <- burden_results_forest %>%
    filter(class == "Rare LoF or Missense") %>%
    filter(case_carriers > 0) %>%
    slice_min(order_by = p_value, n = top_n, with_ties = FALSE) %>%
    mutate(
      OR_text = ifelse(is.infinite(odds_ratio), "Inf", sprintf("%.2f", odds_ratio)),
      CI_text = paste0("[", sprintf("%.2f", ci_low), "; ", ifelse(is.infinite(ci_high), "Inf", sprintf("%.2f", ci_high)), "]"),
      P_raw_text = ifelse(p_value < 0.001, formatC(p_value, format = "e", digits = 2), sprintf("%.3f", p_value)),
      Bonferroni_text = case_when(
        p_adj_bonferroni < 0.001 ~ paste0(formatC(p_adj_bonferroni, format = "e", digits = 2), "*"),
        p_adj_bonferroni < 0.05 ~ paste0(sprintf("%.3f", p_adj_bonferroni), "*"),
        TRUE ~ sprintf("%.3f", p_adj_bonferroni)
      ),
      Case_text = paste0(case_carriers, "/", n_case_plot),
      Control_text = paste0(control_carriers, "/", n_control_plot),
      odds_ratio_plot = pmin(pmax(ifelse(is.infinite(odds_ratio), plot_max, odds_ratio), plot_min), plot_max),
      ci_low_plot = pmin(pmax(ifelse(is.na(ci_low), plot_min, ci_low), plot_min), plot_max),
      ci_high_plot = pmin(pmax(ifelse(is.infinite(ci_high), plot_max, ci_high), plot_min), plot_max)
    ) %>%
    arrange(p_value)
  
  if (nrow(forest_df) == 0) return(NULL)
  
  tabletext <- rbind(
    c("Gene", "Category", "Class", "Case", "Control", "OR", "95% CI", "P value", "Adj.P value"),
    as.matrix(forest_df[, c("gene", "pathway_group", "class", "Case_text", "Control_text", "OR_text", "CI_text", "P_raw_text", "Bonferroni_text")])
  )
  
  p <- forestplot(
    labeltext = tabletext,
    mean = c(NA, forest_df$odds_ratio_plot),
    lower = c(NA, forest_df$ci_low_plot),
    upper = c(NA, forest_df$ci_high_plot),
    zero = 1,
    xlog = TRUE,
    clip = c(plot_min, plot_max),
    xticks = c(0.05, 0.1, 0.5, 1, 2, 10, 20),
    boxsize = 0.18,
    lineheight = unit(0.45, "cm"),
    colgap = unit(4, "mm"),
    graph.pos = 6,
    xlab = "Odds Ratio",
    title = plot_title,
    is.summary = c(TRUE, rep(FALSE, nrow(forest_df))),
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.75),
      ticks = gpar(cex = 0.8),
      xlab = gpar(cex = 0.9),
      title = gpar(cex = 1.1, fontface = "bold")
    ),
    col = fpColors(box = "gray40", line = "gray30", summary = "gray20", zero = "black")
  )
  
  print(p)
  invisible(p)
}

run_group_burden <- function(df, n_case, n_ctrl) {
  
  rare <- df %>% filter(is.na(max_pop_af) | max_pop_af < af_threshold)
  
  rare_burden <- rare %>%
    filter(variant_class %in% c("LoF", "Missense")) %>%
    mutate(class = "Rare LoF or Missense")
  
  group_df <- rare_burden %>%
    mutate(Gene.refGene = str_replace_all(Gene.refGene, "\\\\x3b", ";")) %>%
    separate_rows(Gene.refGene, sep = ";|,") %>%
    mutate(Gene.refGene = str_trim(Gene.refGene)) %>%
    filter(!is.na(Gene.refGene), Gene.refGene != ".", Gene.refGene != "") %>%
    rename(gene = Gene.refGene) %>%
    inner_join(manual_groups, by = "gene")
  
  carrier <- group_df %>%
    distinct(sample_id, group, pathway_group, class) %>%
    count(pathway_group, class, group, name = "n")
  
  wide <- carrier %>%
    pivot_wider(names_from = group, values_from = n, values_fill = 0)
  
  if (!"case" %in% colnames(wide)) wide$case <- 0
  if (!"control" %in% colnames(wide)) wide$control <- 0
  
  wide <- wide %>%
    mutate(
      case_carriers = case,
      control_carriers = control,
      case_noncarriers = n_case - case,
      control_noncarriers = n_ctrl - control
    )
  
  wide %>%
    rowwise() %>%
    mutate(
      fisher = list(run_fisher(case_carriers, case_noncarriers, control_carriers, control_noncarriers))
    ) %>%
    unnest(fisher) %>%
    ungroup() %>%
    group_by(class) %>%
    mutate(p_adj_bonferroni = p.adjust(p_value, method = "bonferroni")) %>%
    ungroup() %>%
    arrange(p_value)
}

child_group_res <- run_group_burden(child_df, n_case_child, n_control_child)
parent_group_res <- run_group_burden(parent_df, n_case_parent, n_control_parent)

make_group_forestplot <- function(group_results,
                                  n_case_plot,
                                  n_control_plot,
                                  plot_title,
                                  top_n = 20) {
  
  plot_min <- 0.05
  plot_max <- 50
  
  forest_df <- group_results %>%
    filter(class == "Rare LoF or Missense") %>%
    arrange(p_value) %>%
    slice_head(n = top_n) %>%
    mutate(
      OR_text = ifelse(is.infinite(odds_ratio), "Inf", sprintf("%.2f", odds_ratio)),
      CI_text = paste0("[", sprintf("%.2f", ci_low), "; ", ifelse(is.infinite(ci_high), "Inf", sprintf("%.2f", ci_high)), "]"),
      RawP_text = ifelse(p_value < 0.001, formatC(p_value, format = "e", digits = 2), sprintf("%.3f", p_value)),
      Bonferroni_text = case_when(
        p_adj_bonferroni < 0.001 ~ paste0(formatC(p_adj_bonferroni, format = "e", digits = 2), "*"),
        p_adj_bonferroni < 0.05 ~ paste0(sprintf("%.3f", p_adj_bonferroni), "*"),
        TRUE ~ sprintf("%.3f", p_adj_bonferroni)
      ),
      Case_text = paste0(case_carriers, "/", n_case_plot),
      Control_text = paste0(control_carriers, "/", n_control_plot),
      odds_ratio_plot = pmin(pmax(ifelse(is.infinite(odds_ratio), plot_max, odds_ratio), plot_min), plot_max),
      ci_low_plot = pmin(pmax(ifelse(is.na(ci_low), plot_min, ci_low), plot_min), plot_max),
      ci_high_plot = pmin(pmax(ifelse(is.infinite(ci_high), plot_max, ci_high), plot_min), plot_max)
    ) %>%
    arrange(p_value)
  
  if (nrow(forest_df) == 0) return(NULL)
  
  tabletext <- rbind(
    c("Pathway group", "Class", "Case", "Control", "OR", "95% CI", "P value", "Adj.P value"),
    as.matrix(forest_df[, c("pathway_group", "class", "Case_text", "Control_text", "OR_text", "CI_text", "RawP_text", "Bonferroni_text")])
  )
  
  p <- forestplot(
    labeltext = tabletext,
    mean = c(NA, forest_df$odds_ratio_plot),
    lower = c(NA, forest_df$ci_low_plot),
    upper = c(NA, forest_df$ci_high_plot),
    zero = 1,
    xlog = TRUE,
    clip = c(plot_min, plot_max),
    xticks = c(0.05, 0.1, 0.5, 1, 2, 10, 50),
    boxsize = 0.18,
    lineheight = unit(0.55, "cm"),
    colgap = unit(4, "mm"),
    graph.pos = 5,
    xlab = "Odds Ratio",
    title = plot_title,
    is.summary = c(TRUE, rep(FALSE, nrow(forest_df))),
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.75),
      ticks = gpar(cex = 0.8),
      xlab = gpar(cex = 0.9),
      title = gpar(cex = 1.1, fontface = "bold")
    ),
    col = fpColors(box = "gray40", line = "gray30", summary = "gray20", zero = "black")
  )
  
  print(p)
  invisible(p)
}

p_child_gene <- make_forestplot(
  child_res,
  n_case_child,
  n_control_child,
  paste0("Case-control gene burden test: child | AF < ", af_threshold)
)

p_parent_gene <- make_forestplot(
  parent_res,
  n_case_parent,
  n_control_parent,
  paste0("Case-control gene burden test: parent | AF < ", af_threshold)
)

p_child_group <- make_group_forestplot(
  child_group_res,
  n_case_child,
  n_control_child,
  paste0("Case-control pathway burden test: child | AF < ", af_threshold)
)

p_parent_group <- make_group_forestplot(
  parent_group_res,
  n_case_parent,
  n_control_parent,
  paste0("Case-control pathway burden test: parent | AF < ", af_threshold)
)

case_target_variants <- all_variants %>%
  mutate(
    Gene.refGene = str_replace_all(Gene.refGene, "\\\\x3b", ";"),
    AAChange.refGene = str_replace_all(AAChange.refGene, "\\\\x3b", ";"),
    Func.refGene = str_replace_all(Func.refGene, "\\\\x3b", ";")
  ) %>%
  separate_rows(Gene.refGene, sep = ";|,") %>%
  mutate(Gene.refGene = str_trim(Gene.refGene)) %>%
  filter(
    group == "case",
    Gene.refGene %in% target_genes,
    variant_class %in% c("Missense", "LoF"),
    is.na(max_pop_af) | max_pop_af < af_threshold
  ) %>%
  rename(gene = Gene.refGene) %>%
  left_join(manual_groups, by = "gene") %>%
  relocate(pathway_group, .after = gene)

fwrite(
  case_target_variants,
  file.path(outdir, "case_manual_target_gene_rare_LoF_Missense_variants.tsv"),
  sep = "\t"
)

fwrite(
  child_res,
  file.path(outdir, "manual_gene_child_burden_results.tsv"),
  sep = "\t"
)

fwrite(
  parent_res,
  file.path(outdir, "manual_gene_parent_burden_results.tsv"),
  sep = "\t"
)

fwrite(
  child_group_res,
  file.path(outdir, "manual_group_child_burden_results.tsv"),
  sep = "\t"
)

fwrite(
  parent_group_res,
  file.path(outdir, "manual_group_parent_burden_results.tsv"),
  sep = "\t"
)

message("Saved data and burden test tables.")
message("Done.")

# ---- Save forest plots as SVG ----
# install once if needed:
save_forest_pdf <- function(plot_fn, args, file, width = 12, height = 6) {
  pdf(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  do.call(plot_fn, args)
}

save_forest_pdf(make_forestplot,
                list(child_res, n_case_child, n_control_child,
                     paste0("Case-control gene burden test: child | AF < ", af_threshold)),
                file.path(outdir, "forest_child_gene.pdf"), 12, 6)

save_forest_pdf(make_forestplot,
                list(parent_res, n_case_parent, n_control_parent,
                     paste0("Case-control gene burden test: parent | AF < ", af_threshold)),
                file.path(outdir, "forest_parent_gene.pdf"), 12, 6)

save_forest_pdf(make_group_forestplot,
                list(child_group_res, n_case_child, n_control_child,
                     paste0("Case-control pathway burden test: child | AF < ", af_threshold)),
                file.path(outdir, "forest_child_pathway.pdf"), 11, 5)

save_forest_pdf(make_group_forestplot,
                list(parent_group_res, n_case_parent, n_control_parent,
                     paste0("Case-control pathway burden test: parent | AF < ", af_threshold)),
                file.path(outdir, "forest_parent_pathway.pdf"), 11, 5)

message("Saved four PDF forest plots to: ", outdir)
