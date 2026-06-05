# ---- Packages ----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- Inputs ----
pop_betas <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/methlyation_analysis/metafora/METAFORA_new/old_12292025/METAFORA_1-02917/METAFORA_output_control/plot/Population_methylation.tissue_Blood.chrom_chr16.1.betas.mat.txt"

# 3 sample methylation beds to overlay
sample_beds <- c(
  "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/methlyation_analysis/metafora/METAFORA_new/old_12292025/METAFORA_1-14785/METAFORA_output_control/plot/PCGC_1-14785.tech_PacBio.cpg_methylation.combined.bed.txt",
  "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/methlyation_analysis/metafora/METAFORA_new/old_12292025/METAFORA_1-06216/METAFORA_output_control/plot/PCGC_1-06216.tech_PacBio.cpg_methylation.combined.bed.txt",
  "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/methlyation_analysis/metafora/METAFORA_new/old_12292025/METAFORA_1-02917/METAFORA_output_control/plot/PCGC_1-02917.tech_PacBio.cpg_methylation.combined.bed.txt"
)

# Explicit labels (for legend + ordering)
sample_labels <- c("1-14785", "1-06216", "1-02917")

# Color-blind safe warm tones (no yellow)
sample_colors <- c(
  "1-14785" = "#D55E00",  # vermillion (deep orange-red)
  "1-06216" = "#B2182B",  # dark crimson
  "1-02917" = "#7B3294"   # muted purple
)

chrom <- "chr16"

# Regions of interest (plot windows)
regions <- list(
  c(3442680L, 3444679L)
)

# Highlight regions (shaded)
highlight_regions <- list(
  c(3442827L, 3444461L)
)

# ---- Load population beta matrix ----
cat("Loading population matrix...\n")
pop <- fread(pop_betas)
cat("Population matrix loaded with", nrow(pop), "rows and", ncol(pop), "columns.\n")

stopifnot(all(c("chromosome", "start", "end") %in% names(pop)))

# Compute population mean & SD across sample columns (col 4+)
pop$pop_mean <- rowMeans(pop[, 4:ncol(pop), drop = FALSE], na.rm = TRUE)
pop$pop_sd   <- apply(pop[, 4:ncol(pop), drop = FALSE], 1, sd, na.rm = TRUE)

# ---- Helper: read & normalize ONE pb-cpg-tools combined.bed ----
read_one_sample <- function(bed_path, sample_label) {
  cat("Loading sample methylation:", sample_label, "\n")
  samp <- fread(bed_path)
  
  # Normalize pb-cpg-tools column names
  if ("#chrom" %in% names(samp)) setnames(samp, "#chrom", "chromosome")
  if ("chrom"  %in% names(samp)) setnames(samp, "chrom",  "chromosome")
  if ("begin"  %in% names(samp)) setnames(samp, "begin",  "start")
  if ("stop"   %in% names(samp)) setnames(samp, "stop",   "end")
  
  # Keep only overall methylation rows
  if ("type" %in% names(samp)) samp <- dplyr::filter(samp, type == "Total")
  
  # Compute beta if not present (mod_score is 0..100)
  if (!"beta" %in% names(samp) && "mod_score" %in% names(samp)) {
    samp <- dplyr::mutate(samp, beta = pmin(pmax(mod_score / 100, 0), 1))
  }
  
  stopifnot(all(c("chromosome", "start", "end", "beta") %in% names(samp)))
  
  samp %>%
    dplyr::select(chromosome, start, beta) %>%
    dplyr::distinct() %>%
    dplyr::mutate(sample_label = sample_label)
}

# ---- Load ALL sample beds and stack long ----
samp_all <- dplyr::bind_rows(
  Map(read_one_sample, sample_beds, sample_labels)
)

# enforce factor order for legend
samp_all$sample_label <- factor(samp_all$sample_label, levels = sample_labels)
cat("All samples combined:", length(unique(as.character(samp_all$sample_label))), "samples\n")

# ---- Region prep (join pop to ALL samples by chromosome+start) ----
prep_region_df_multi <- function(pop, samp_all, chr, start_pos, end_pos) {
  pop_sub <- pop %>%
    dplyr::filter(.data$chromosome == chr,
                  .data$start >= start_pos, .data$start <= end_pos) %>%
    dplyr::select(chromosome, start, end, pop_mean, pop_sd)
  
  samp_sub <- samp_all %>%
    dplyr::filter(.data$chromosome == chr,
                  .data$start >= start_pos, .data$start <= end_pos)
  
  pop_sub %>%
    dplyr::left_join(samp_sub, by = c("chromosome", "start")) %>%
    dplyr::mutate(
      position    = .data$start,
      sample_beta = .data$beta
    ) %>%
    dplyr::arrange(.data$position, .data$sample_label)
}

# ---- Plot function (title shows highlighted region coords) ----
plot_region_multi <- function(df, chr, start_pos, end_pos,
                              highlight_regions, palette, sample_labels) {
  
  df$sample_label <- factor(df$sample_label, levels = sample_labels)
  
  # Build highlight label string (supports multiple highlight intervals)
  highlight_label <- paste(
    sprintf(
      "%s:%s–%s",
      chr,
      format(sapply(highlight_regions, `[[`, 1), big.mark = ","),
      format(sapply(highlight_regions, `[[`, 2), big.mark = ",")
    ),
    collapse = ", "
  )
  
  ggplot(df, aes(x = position)) +
    
    # Highlight regions
    geom_rect(
      data = data.frame(
        xmin = sapply(highlight_regions, `[[`, 1),
        xmax = sapply(highlight_regions, `[[`, 2),
        ymin = 0, ymax = 1
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "yellow", alpha = 0.25
    ) +
    
    # Population envelope and mean
    geom_ribbon(
      aes(ymin = pmax(pop_mean - pop_sd, 0),
          ymax = pmin(pop_mean + pop_sd, 1)),
      fill = "grey70", alpha = 0.4
    ) +
    geom_line(aes(y = pop_mean), color = "black", linewidth = 1.1) +
    
    # Sample lines + points
    geom_line(aes(y = sample_beta, color = sample_label),
              linewidth = 0.9, alpha = 0.95, na.rm = TRUE) +
    geom_point(aes(y = sample_beta, color = sample_label),
               size = 0.6, alpha = 0.9, na.rm = TRUE) +
    
    scale_color_manual(
      name = "Sample",
      values = palette,
      breaks = sample_labels
    ) +
    
    scale_x_continuous(labels = scales::comma) +
    coord_cartesian(xlim = c(start_pos, end_pos), ylim = c(0, 1)) +
    
    labs(
      title = "Population vs 3 samples methylation",
      subtitle = paste(highlight_label),
      x = sprintf("%s Position", chr),
      y = "Methylation",
      caption = "Black: population mean • Gray: mean ± SD • Colored: samples • Yellow: Methlayation Outlier"
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
}

# ---- Run & Plot ----
plots <- vector("list", length(regions))

for (i in seq_along(regions)) {
  r <- regions[[i]]
  cat(sprintf("Plotting region %s:%d-%d\n", chrom, r[1], r[2]))
  
  df <- prep_region_df_multi(pop, samp_all, chrom, r[1], r[2])
  
  plots[[i]] <- plot_region_multi(
    df = df,
    chr = chrom,
    start_pos = r[1],
    end_pos = r[2],
    highlight_regions = highlight_regions,
    palette = sample_colors,
    sample_labels = sample_labels
  )
}

print(plots[[1]])

# ---- Save ----
ggsave(
  filename = "/Users/nanakong/Library/CloudStorage/Box-Box/Uniparental Disomy in Congenital Heart Disease/Figure/Fig5_source/Fig5B_chr16_case_control_methly.svg",
  plot = plots[[1]],
  device = "svg",
  width = 5,
  height = 6,
  units = "in"
)
