# ---- Packages ----
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# ---- Inputs ----
pop_betas <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/methlyation_analysis/metafora/METAFORA_new/old_12292025/METAFORA_1-03820/METAFORA_output_control/plot/Population_methylation.tissue_Blood.chrom_chr15.1.betas.mat.txt"
sample_bed <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/methlyation_analysis/metafora/METAFORA_new/old_12292025/METAFORA_1-03820/METAFORA_output_control/plot/PCGC_1-03820.tech_PacBio.cpg_methylation.combined.bed.txt"

sample_label <- "PCGC_1-03820"
chrom <- "chr15"

# Regions of interest
regions <- list(
#  c(3442680L, 3444679L),
#  c(29317409L, 29318922L),
  c(24953500L, 24958000L)
  
)

# ---- Load population beta matrix ----
cat("Loading population matrix...\n")
pop <- fread(pop_betas)
cat("Population matrix loaded with", nrow(pop), "rows and", ncol(pop), "columns.\n")

# Sanity check
stopifnot(all(c("chromosome", "start", "end") %in% names(pop)))

# Compute population mean & SD across sample columns (col 4+)
pop$pop_mean <- rowMeans(pop[, 4:ncol(pop), drop = FALSE], na.rm = TRUE)
pop$pop_sd   <- apply(pop[, 4:ncol(pop), drop = FALSE], 1, sd, na.rm = TRUE)

# ---- Load & normalize sample methylation (pb-cpg-tools combined.bed) ----
cat("Loading sample methylation file...\n")
samp <- fread(sample_bed)

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

cat("Sample methylation loaded with", nrow(samp), "rows and", ncol(samp), "columns.\n")

# Verify required columns
stopifnot(all(c("chromosome", "start", "end", "beta") %in% names(samp)))

# (Optional) If you prefer a 3-col join, uncomment the next line to normalize to single-base intervals:
# samp <- dplyr::mutate(samp, end = ifelse(end == start + 1L, start, end))

# ---- Region prep (join by chromosome+start to handle 1-bp intervals) ----
prep_region_df <- function(pop, samp, chr, start_pos, end_pos) {
  # population subset
  pop_sub <- pop %>%
    dplyr::filter(.data$chromosome == chr,
                  .data$start >= start_pos, .data$start <= end_pos) %>%
    dplyr::select(chromosome, start, end, pop_mean, pop_sd)
  
  # sample subset (join by start)
  samp_sub <- samp %>%
    dplyr::filter(.data$chromosome == chr,
                  .data$start >= start_pos, .data$start <= end_pos) %>%
    dplyr::select(chromosome, start, beta) %>%
    dplyr::distinct()
  
  pop_sub %>%
    dplyr::left_join(samp_sub, by = c("chromosome", "start")) %>%
    dplyr::mutate(
      position    = .data$start,
      sample_beta = .data$beta,
      delta       = .data$sample_beta - .data$pop_mean
    ) %>%
    dplyr::arrange(.data$position)
}

# Regions of interest
regions <- list(
#  c(3442680L, 3444679L),
#  c(29317409L, 29318922L)
  c(24953500L, 24958000L)
)

# Highlight regions (yellow)
highlight_regions <- list(
#  c(3442827L, 3444461L),
#  c(29317783L, 29318549L)
  c(24954396L, 24956827L)
)
# ---- Plot function ----
plot_region <- function(df, chr, start_pos, end_pos, sample_label, highlight_regions) {
  g <- ggplot(df, aes(x = position)) +
    # Highlight regions
    geom_rect(
      data = data.frame(
        xmin = sapply(highlight_regions, `[[`, 1),
        xmax = sapply(highlight_regions, `[[`, 2),
        ymin = 0, ymax = 1
      ),
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = "yellow", alpha = 0.3
    ) +
    geom_ribbon(aes(ymin = pmax(pop_mean - pop_sd, 0),
                    ymax = pmin(pop_mean + pop_sd, 1)),
                alpha = 0.2) +
    geom_line(aes(y = pop_mean), linewidth = 1.0, color = "black") +
    geom_line(aes(y = sample_beta), color = "red", linewidth = 0.7, na.rm = TRUE) +
    geom_point(aes(y = sample_beta), color = "red", size = 0.6, na.rm = TRUE) +
    scale_x_continuous(
      labels = label_comma(),
      expand = expansion(mult = 0, add = 0)
    ) +
    coord_cartesian(xlim = c(start_pos, end_pos), ylim = c(0, 1)) +
    coord_cartesian(xlim = c(start_pos, end_pos), ylim = c(0, 1)) +
    labs(
      title = sprintf("Population vs. %s methylation\n%s:%s-%s",
                      sample_label, chr,
                      format(start_pos, big.mark=","), format(end_pos, big.mark=",")),
      x = sprintf("%s Position", chr),
      y = "Methylation β",
      caption = "Gray: population mean ± SD   •   Red: observed sample β   •   Yellow: highlighted interval"
    ) +
    theme_bw(base_size = 12)
  g
}

# ---- Run & Plot ----
plots <- vector("list", length(regions))
for (i in seq_along(regions)) {
  r <- regions[[i]]
  cat(sprintf("Plotting region %s:%d-%d\n", chrom, r[1], r[2]))
  df <- prep_region_df(pop, samp, chrom, r[1], r[2])
  plots[[i]] <- plot_region(df, chrom, r[1], r[2], sample_label, highlight_regions)
}

# ---- Display ----
print(plots[[1]])
# ---- Save ----
ggsave(
  filename = "/Users/nanakong/Library/CloudStorage/Box-Box/Uniparental Disomy in Congenital Heart Disease/Figure/Fig5_source/Fig5D_chr15_case_control_methly.svg",
  plot = plots[[1]],
  device = "svg",
  width = 5,
  height = 6,
  units = "in",
  dpi = 300
)
