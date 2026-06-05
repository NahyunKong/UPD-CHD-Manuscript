###############################################################################
# Lollipop plot — synaptonemal-complex variants (SYCP1, SYCP2, C14orf39)
# - fetches real domains from UniProt
# - all genes on one combined figure
# - coloured by consequence, sized by # samples, labelled with sample IDs
# - 12pt Helvetica, saved as SVG at 5 x 3 inches
###############################################################################

# ---- 0. packages ------------------------------------------------------------
pkgs <- c("httr", "jsonlite", "ggplot2", "dplyr", "stringr",
          "patchwork", "tidyr", "svglite")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1. your variants (one row per sample observation) ----------------------
variants_raw <- tribble(
  ~gene,      ~uniprot,  ~aa,  ~label,       ~type,        ~sample,      ~role,
  "SYCP2",    "Q9BX26",  454,  "I454V",      "Missense",   "1-00660",    "child",
  "SYCP2",    "Q9BX26",  454,  "I454V",      "Missense",   "1-00660-01", "mother",
  "SYCP2",    "Q9BX26",  1021, "T1021K",     "Missense",   "1-03820",    "child",
  "SYCP2",    "Q9BX26",  1019, "T1019N",     "Missense",   "1-03820",    "child",
  "SYCP1",    "Q15431",  629,  "G629Vfs*8",  "Truncating", "1-03573-02", "father",
  "C14orf39", "Q8IYX8",  138,  "Y138X",      "Truncating", "1-06216",    "child",
  "C14orf39", "Q8IYX8",  138,  "Y138X",      "Truncating", "1-06216-02", "father",
  "C14orf39", "Q8IYX8",  313,  "N313S",      "Missense",   "1-14523",    "child",
  "C14orf39", "Q8IYX8",  313,  "N313S",      "Missense",   "1-14523-01", "mother",
  "SYCP1",    "Q15431",  583,  "L583V",      "Missense",   "1-15162",    "child",
  "SYCP1",    "Q15431",  583,  "L583V",      "Missense",   "1-15162-02", "father"
)

variants <- variants_raw %>%
  mutate(samp_role = paste0(sample, " (", role, ")")) %>%
  group_by(gene, uniprot, aa, label, type) %>%
  summarise(n_samples = n(),
            samples   = paste(samp_role, collapse = "\n"),
            .groups   = "drop") %>%
  mutate(full_label = paste0(label, "\n", samples))

gene_uniprot <- c(SYCP1 = "Q15431", SYCP2 = "Q9BX26", C14orf39 = "Q8IYX8")
fallback_len <- c(SYCP1 = 976, SYCP2 = 1530, C14orf39 = 587)

# ---- 2. fetch protein length + domains from UniProt -------------------------
get_uniprot <- function(acc) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", acc, ".json")
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(res) || status_code(res) != 200) {
    message("UniProt fetch failed for ", acc, " — using fallback length, no domains.")
    return(list(length = NA, features = NULL))
  }
  dat  <- fromJSON(content(res, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  plen <- dat$sequence$length
  keep <- c("Domain", "Region", "Coiled coil", "Repeat",
            "Zinc finger", "DNA binding", "Motif")
  rows <- lapply(dat$features, function(f) {
    if (is.null(f$type) || !(f$type %in% keep)) return(NULL)
    s <- f$location$start$value; e <- f$location$end$value
    if (is.null(s) || is.null(e)) return(NULL)
    desc <- if (!is.null(f$description) && nzchar(f$description)) f$description else f$type
    data.frame(type = f$type, start = s, end = e, desc = desc,
               stringsAsFactors = FALSE)
  })
  rows <- do.call(rbind, rows)
  list(length = plen, features = rows)
}

info <- lapply(gene_uniprot, get_uniprot)

proteins <- tibble(
  gene = names(gene_uniprot),
  uniprot = unname(gene_uniprot),
  length = sapply(names(gene_uniprot), function(g) {
    L <- info[[g]]$length
    if (is.null(L) || is.na(L)) fallback_len[[g]] else L
  })
)

domains <- bind_rows(lapply(names(info), function(g) {
  f <- info[[g]]$features
  if (is.null(f)) return(NULL)
  f$gene <- g; f
}))

# ---- 3. plotting (12pt Helvetica, fixed lollipop size) ----------------------
type_cols <- c("Missense" = "#2c7fb8", "Truncating" = "#d7301f")

feat_levels <- if (nrow(domains)) sort(unique(domains$type)) else character(0)
feat_pal    <- setNames(grDevices::hcl.colors(max(length(feat_levels), 1), "Set2"),
                        feat_levels)

if (nrow(domains)) domains$type <- factor(domains$type, levels = feat_levels)
variants$type <- factor(variants$type, levels = names(type_cols))

make_panel <- function(g) {
  pr <- proteins %>% filter(gene == g)
  vr <- variants %>% filter(gene == g)
  dm <- if (nrow(domains)) domains %>% filter(gene == g) else domains[0, ]
  
  vr <- vr %>% arrange(aa) %>%
    mutate(stem_top = 1.0,
           lab_y    = 1.25 + 0.55 * (row_number() %% 2))
  
  ggplot() +
    geom_rect(data = pr,
              aes(xmin = 1, xmax = length, ymin = 0.42, ymax = 0.58),
              fill = "grey85", colour = "grey55", linewidth = 0.3) +
    { if (nrow(dm))
      geom_rect(data = dm,
                aes(xmin = start, xmax = end, ymin = 0.30, ymax = 0.70,
                    fill = type),
                colour = "grey30", linewidth = 0.2, alpha = 0.9) } +
    geom_segment(data = vr,
                 aes(x = aa, xend = aa, y = 0.58, yend = stem_top),
                 colour = "grey40", linewidth = 0.4) +
    # fixed size = 3 (no size aesthetic, so no # samples legend)
    geom_point(data = vr,
               aes(x = aa, y = stem_top, colour = type), size = 3) +
    geom_text(data = vr,
              aes(x = aa, y = lab_y, label = full_label),
              size = 12/.pt, vjust = 0, lineheight = 0.9, family = "Helvetica") +
    scale_colour_manual(values = type_cols, name = "Consequence",
                        limits = names(type_cols), drop = FALSE,
                        na.translate = FALSE) +
    scale_fill_manual(values = feat_pal, name = "Domain / region",
                      limits = feat_levels, drop = FALSE,
                      na.translate = FALSE) +
    scale_x_continuous(limits = c(1, pr$length),
                       expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(0.25, 2.1)) +
    labs(title = paste0(g, " (", pr$uniprot, ", ", pr$length, " aa)")) +
    xlab("Amino-acid position") +
    theme_minimal(base_size = 12, base_family = "Helvetica") +
    theme(axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid   = element_blank(),
          plot.title   = element_text(face = "bold", size = 12),
          axis.title.x = element_text(size = 8),
          axis.text.x  = element_text(size = 8),
          legend.title = element_text(size = 8),
          legend.text  = element_text(size = 8),
          legend.key.size = unit(0.4, "cm"))
}

panels <- lapply(proteins$gene, make_panel)

combined <- wrap_plots(panels, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Synaptonemal-complex gene variants by sample",
    theme = theme(plot.title = element_text(face = "bold", size = 12,
                                            family = "Helvetica"))
  ) &
  theme(legend.position = "right",
        text = element_text(family = "Helvetica"))

print(combined)

ggsave("/Users/nanakong/Library/CloudStorage/Box-Box/Uniparental Disomy in Congenital Heart Disease/Figure/Fig2_source/synaptonemal_lollipop_samples.svg", combined,
       device = svglite::svglite,
       width = 10, height = 6, units = "in")
