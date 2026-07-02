if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("karyoploteR", quietly = TRUE))
  BiocManager::install("karyoploteR")

if (!requireNamespace("matrixStats", quietly = TRUE))
  install.packages("matrixStats")

library(karyoploteR)
library(GenomicRanges)

iso_mat_color <- "#FFB6C1"
het_mat_color <- "#CC0000"
iso_pat_color <- "#ADD8E6"
het_pat_color <- "#00008B"
bg_color      <- "#E8E4D0"

wes_dir <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20240716_final_submission/triomix_results0525/new_breakpoint2/wes"
wgs_dir <- "/Volumes/jin810/Active/testing/Nahyun/UPD_PROJECT/20240716_final_submission/triomix_results0525/new_breakpoint/wgs_snp"
wgs_samples <- c("1-00180", "1-02917", "1-04975")

read_sample <- function(bed_file, sample_id) {
  df <- read.table(
    bed_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  
  if (nrow(df) == 0) return(NULL)
  
  df <- df[!is.na(df$UPD_subtype) &
             df$UPD_subtype %in% c("isodisomy", "heterodisomy"), ]
  
  if (nrow(df) == 0) return(NULL)
  
  df$type <- ifelse(df$UPD_subtype == "isodisomy", "iso", "het")
  df$origin <- ifelse(grepl("^paternal", df$UPD_final), "pat", "mat")
  
  make_gr <- function(subset_df) {
    if (nrow(subset_df) == 0) return(NULL)
    
    gr_df <- data.frame(
      chr = subset_df$chromosome,
      start = subset_df$start,
      end = subset_df$end
    )
    
    toGRanges(gr_df)
  }
  
  list(
    id      = sample_id,
    df      = df,
    chrs    = unique(df$chromosome),
    iso_mat = make_gr(df[df$type == "iso" & df$origin == "mat", ]),
    het_mat = make_gr(df[df$type == "het" & df$origin == "mat", ]),
    iso_pat = make_gr(df[df$type == "iso" & df$origin == "pat", ]),
    het_pat = make_gr(df[df$type == "het" & df$origin == "pat", ])
  )
}

# Load WES samples
samples <- list()

bed_files <- list.files(
  wes_dir,
  pattern = "\\.bed$",
  recursive = TRUE,
  full.names = TRUE
)

for (bf in bed_files) {
  sid <- basename(dirname(bf))
  s <- read_sample(bf, sid)
  
  if (!is.null(s)) {
    samples[[sid]] <- s
  }
}

# Override selected samples with WGS SNP results
for (sid in wgs_samples) {
  bf <- list.files(
    file.path(wgs_dir, sid),
    pattern = "\\.bed$",
    full.names = TRUE
  )
  
  if (length(bf) == 0) {
    warning(paste("No BED for", sid, "in wgs_snp"))
    next
  }
  
  s <- read_sample(bf[1], sid)
  
  if (!is.null(s)) {
    samples[[sid]] <- s
  }
}

message("Loaded samples:")
print(names(samples))

# Build chromosome-to-sample map
chr_samples <- list()

for (sid in names(samples)) {
  for (chr in samples[[sid]]$chrs) {
    chr_samples[[chr]] <- c(chr_samples[[chr]], sid)
  }
}

message("Chromosomes with UPD:")
print(names(chr_samples))

TRACK_H <- 0.28
TRACK_G <- 0.08

track_r <- function(case_idx) {
  r0 <- case_idx * (TRACK_H + TRACK_G)
  list(r0 = r0, r1 = r0 + TRACK_H)
}

custom_cytobands <- function() {
  cb <- getCytobands("hg19")
  orig <- cb$gieStain
  
  cb$gieStain <- "gneg"
  cb$gieStain[orig == "acen"] <- "acen"
  
  cb
}

custom_colors <- c(
  gneg    = "#CCCCCC",
  gpos25  = "#CCCCCC",
  gpos50  = "#CCCCCC",
  gpos75  = "#CCCCCC",
  gpos100 = "#CCCCCC",
  acen    = "#333333",
  gvar    = "#CCCCCC",
  stalk   = "#CCCCCC"
)

draw_case <- function(kp, chr, label,
                      iso_mat, het_mat,
                      iso_pat, het_pat,
                      r0, r1) {
  
  chr_len <- seqlengths(kp$genome)[chr]
  chr_gr <- toGRanges(data.frame(chr = chr, start = 1, end = chr_len))
  
  mid <- (r0 + r1) / 2
  half <- (r1 - r0) * 0.20
  
  top_r0 <- mid + half * 0.1
  top_r1 <- mid + half * 1.1
  bot_r0 <- mid - half * 1.1
  bot_r1 <- mid - half * 0.1
  
  kpPlotRegions(
    kp,
    data = chr_gr,
    r0 = r0,
    r1 = r1,
    col = bg_color,
    border = NA
  )
  
  if (!is.null(iso_mat)) {
    kpPlotRegions(kp, data = iso_mat, r0 = top_r0, r1 = top_r1,
                  col = iso_mat_color, border = NA)
    kpPlotRegions(kp, data = iso_mat, r0 = bot_r0, r1 = bot_r1,
                  col = iso_mat_color, border = NA)
  }
  
  if (!is.null(het_mat)) {
    kpPlotRegions(kp, data = het_mat, r0 = top_r0, r1 = top_r1,
                  col = iso_mat_color, border = NA)
    kpPlotRegions(kp, data = het_mat, r0 = bot_r0, r1 = bot_r1,
                  col = het_mat_color, border = NA)
  }
  
  if (!is.null(iso_pat)) {
    kpPlotRegions(kp, data = iso_pat, r0 = top_r0, r1 = top_r1,
                  col = iso_pat_color, border = NA)
    kpPlotRegions(kp, data = iso_pat, r0 = bot_r0, r1 = bot_r1,
                  col = iso_pat_color, border = NA)
  }
  
  if (!is.null(het_pat)) {
    kpPlotRegions(kp, data = het_pat, r0 = top_r0, r1 = top_r1,
                  col = iso_pat_color, border = NA)
    kpPlotRegions(kp, data = het_pat, r0 = bot_r0, r1 = bot_r1,
                  col = het_pat_color, border = NA)
  }
  
  kpText(
    kp,
    data = chr_gr,
    x = 0,
    y = r1 + 0.03,
    labels = label,
    cex = 0.40,
    pos = 4,
    col = "black"
  )
}

draw_upd_plot <- function() {
  
  if (length(samples) == 0) {
    stop("No samples loaded. Check BED files and UPD_subtype column.")
  }
  
  pp <- getDefaultPlotParams(plot.type = 2)
  pp$leftmargin     <- 0.06
  pp$rightmargin    <- 0.02
  pp$topmargin      <- 20
  pp$bottommargin   <- 20
  pp$ideogramheight <- 15
  pp$data1height    <- 500
  pp$data1inmargin  <- 10
  pp$data2height    <- 500
  pp$data2inmargin  <- 10
  pp$between.panel  <- 25
  
  kp <- plotKaryotype(
    plot.type      = 2,
    chromosomes    = "canonical",
    plot.params    = pp,
    cytobands      = custom_cytobands(),
    cytobandColors = custom_colors
  )
  
  for (chr in names(chr_samples)) {
    sids <- chr_samples[[chr]]
    
    for (i in seq_along(sids)) {
      sid <- sids[[i]]
      s <- samples[[sid]]
      
      filter_chr <- function(gr) {
        if (is.null(gr)) return(NULL)
        out <- gr[seqnames(gr) == chr]
        if (length(out) == 0) NULL else out
      }
      
      iso_mat <- filter_chr(s$iso_mat)
      het_mat <- filter_chr(s$het_mat)
      iso_pat <- filter_chr(s$iso_pat)
      het_pat <- filter_chr(s$het_pat)
      
      has_mat <- !is.null(iso_mat) || !is.null(het_mat)
      has_pat <- !is.null(iso_pat) || !is.null(het_pat)
      
      origin <- if (has_mat && has_pat) {
        "mat/pat"
      } else if (has_pat) {
        "pat"
      } else {
        "mat"
      }
      
      has_iso <- !is.null(iso_mat) || !is.null(iso_pat)
      has_het <- !is.null(het_mat) || !is.null(het_pat)
      
      type <- if (has_iso && has_het) {
        "mixed"
      } else if (has_iso) {
        "isodisomy"
      } else {
        "heterodisomy"
      }
      
      label <- paste0(
        sid,
        " (upd(",
        gsub("chr", "", chr),
        ")",
        origin,
        " ",
        type,
        ")"
      )
      
      t <- track_r(i - 1)
      
      draw_case(
        kp = kp,
        chr = chr,
        label = label,
        iso_mat = iso_mat,
        het_mat = het_mat,
        iso_pat = iso_pat,
        het_pat = het_pat,
        r0 = t$r0,
        r1 = t$r1
      )
    }
  }
  
  legend(
    "bottomright",
    legend = c(
      "Mat isodisomy",
      "Mat heterodisomy",
      "Pat isodisomy",
      "Pat heterodisomy"
    ),
    fill = c(
      iso_mat_color,
      het_mat_color,
      iso_pat_color,
      het_pat_color
    ),
    border = NA,
    bty = "n",
    cex = 0.75
  )
}

pdf("/Users/nanakong/Library/CloudStorage/Box-Box/Uniparental Disomy in Congenital Heart Disease/Figure/Fig1_source/UPD_karyotype.svg", width=18, height=14)

draw_upd_plot()
dev.off()
rstudioapi::viewer(normalizePath("/Users/nanakong/Library/CloudStorage/Box-Box/Uniparental Disomy in Congenital Heart Disease/Figure/Fig1_source/UPD_karyotype.svg"))
