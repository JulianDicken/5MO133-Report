library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(patchwork)
library(tidyr)
library(ggrastr)
library(scales)
library(RColorBrewer)

# ── Constants ─────────────────────────────────────────────────────────────────

MULTIPLET_MIN_FRACTION   <- 0.05
RANDOM_ENTROPY_THRESHOLD <- 0.5
RANDOM_MIN_READS         <- 1e3

SPECIES_MAP <- c(
  "ATCC_10987"   = "Bacillus pacificus",
  "ATCC_12228"   = "Staphylococcus epidermidis",
  "ATCC_15703"   = "Bifidobacterium adolescentis",
  "ATCC_17029"   = "Cereibacter sphaeroides",
  "ATCC_33323"   = "Lactobacillus gasseri",
  "ATCC_35702"   = "Clostridium beijerinckii",
  "ATCC_47077"   = "Enterococcus faecalis",
  "ATCC_700610"  = "Streptococcus mutans",
  "ATCC_700926"  = "Escherichia coli",
  "ATCC_BAA-816" = "Deinococcus radiodurans",
  "NZ_CDNB"      = "Xanthomonas campestris",
  "phix"         = "PhiX"
)

SPECIES_EXCLUDE <- c("PhiX")

# ── Load & annotate data ──────────────────────────────────────────────────────

qc.composition.species <- bind_rows(
  read_table("data/mda/qc_composition.tsv") %>% mutate(group = "MDA"),
  read_table("data/pta/qc_composition.tsv") %>% mutate(group = "PTA")
) %>%
  mutate(
    species = case_when(
      reference == "phix"          ~ unname(SPECIES_MAP["phix"]),
      str_starts(reference, "NZ_") ~ unname(SPECIES_MAP["NZ_CDNB"]),
      .default = unname(SPECIES_MAP[str_remove(reference, "_contig_\\d+$")])
    )
  ) %>%
  group_by(group, id, species) %>%
  summarise(countof_reads = sum(countof_reads), .groups = "drop")

# ── QC flags ─────────────────────────────────────────────────────────────────

qc.multiplet <- qc.composition.species %>%
  group_by(group, id) %>%
  mutate(species.fraction = countof_reads / sum(countof_reads)) %>%
  summarise(
    n_contributors = sum(species.fraction > MULTIPLET_MIN_FRACTION),
    .groups = "drop"
  ) %>%
  mutate(is.multiplet = n_contributors >= 2)

qc.entropy <- qc.composition.species %>%
  group_by(group, id) %>%
  mutate(species.fraction = countof_reads / sum(countof_reads)) %>%
  summarise(
    shannon.entropy = -sum(ifelse(species.fraction > 0, species.fraction * log2(species.fraction), 0)) / log2(n()),
    total           = sum(countof_reads),
    .groups = "drop"
  ) %>%
  mutate(is.random = shannon.entropy > RANDOM_ENTROPY_THRESHOLD | total < RANDOM_MIN_READS)

# namespacing required here due to conflict with another package????
qc.dominant <- qc.composition.species %>%
  dplyr::group_by(group, id) %>%
  dplyr::mutate(species.fraction = countof_reads / sum(countof_reads)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, id, dplyr::desc(species.fraction)) %>%
  dplyr::distinct(group, id, .keep_all = TRUE) %>%
  dplyr::rename(dominant.species = species) %>%
  dplyr::select(group, id, dominant.species)

# ── Barnyard data ─────────────────────────────────────────────────────────────

qc.barnyard <- qc.composition.species %>%
  group_by(group, id) %>%
  mutate(fraction = countof_reads / sum(countof_reads)) %>%
  arrange(desc(fraction), .by_group = TRUE) %>%
  summarise(
    dominant.fraction    = max(fraction),
    subdominant.fraction = sort(fraction, decreasing = TRUE)[2],
    remaining.fraction   = 1 - max(fraction),
    .groups = "drop"
  ) %>%
  left_join(qc.multiplet, by = c("group", "id")) %>%
  left_join(qc.entropy,   by = c("group", "id")) %>%
  left_join(qc.dominant,  by = c("group", "id")) %>%
  mutate(
    cell.classification = case_when(
      is.random    ~ "background",
      is.multiplet ~ "multiplet",
      TRUE         ~ "singlet"
    ),
    x.reads = dominant.fraction    * total,
    y.sub   = subdominant.fraction * total,
    y.rem   = remaining.fraction   * total
  )

# ── Contamination heatmap  ────────────────────────────────────────────────

qc.contamination <- qc.composition.species %>%
  left_join(
    qc.barnyard %>% select(group, id, dominant.species, cell.classification),
    by = c("group", "id")
  ) %>%
  filter(
    cell.classification == "singlet",
    species != dominant.species,
    !species %in% SPECIES_EXCLUDE,
    !dominant.species %in% SPECIES_EXCLUDE
  ) %>%
  rename(contaminant.species = species) %>%
  group_by(group, dominant.species, contaminant.species) %>%
  summarise(contamination_reads = sum(countof_reads), .groups = "drop") %>%
  group_by(group, dominant.species) %>%
  mutate(contamination_fraction = contamination_reads / sum(contamination_reads)) %>%
  ungroup()

# ── Pre-log transform for barnyard ───────────────────────────────────────────

qc.barnyard.log <- qc.barnyard %>%
  mutate(
    x.reads = log10(x.reads + 1),
    y.sub   = log10(y.sub   + 1),
    y.rem   = log10(y.rem   + 1)
  )

barnyard.limits <- range(c(
  qc.barnyard.log$x.reads,
  qc.barnyard.log$y.sub,
  qc.barnyard.log$y.rem
), na.rm = TRUE)

barnyard.log.breaks <- seq(floor(barnyard.limits[1]), ceiling(barnyard.limits[2]))
barnyard.scale.x    <- scale_x_continuous(limits = barnyard.limits, breaks = barnyard.log.breaks)
barnyard.scale.y    <- scale_y_continuous(limits = barnyard.limits, breaks = barnyard.log.breaks)

barnyard.species.counts <- qc.barnyard %>%
  count(dominant.species, name = "species.count")

barnyard.centroids <- bind_rows(
  qc.barnyard.log %>%
    group_by(group) %>%
    summarise(x = median(x.reads), y.sub = median(y.sub), y.rem = median(y.rem), .groups = "drop") %>%
    mutate(type = "all"),
  qc.barnyard.log %>%
    filter(cell.classification == "singlet") %>%
    group_by(group) %>%
    summarise(x = median(x.reads), y.sub = median(y.sub), y.rem = median(y.rem), .groups = "drop") %>%
    mutate(type = "singlet")
)

# ── Barnyard plot ─────────────────────────────────────────────────────────────
plot.barnyard <- function(data, y.col, title, x.lab, y.lab) {
  data <- data %>%
    left_join(barnyard.species.counts, by = "dominant.species") %>%
    arrange(desc(species.count))
  
  error.lines <- bind_rows(
    data.frame(
      error.rate = c(0.01, 0.05, 0.10),
      label      = factor(c("1%", "5%", "10%"), levels = c("1%", "5%", "10%")),
      type       = "observed"
    ),
    # worst case sub: 50% of all reads are from other species, worst case rem 90% of reads are from other species (even dist)
    data.frame(
      error.rate = c(0.50, 0.90),
      label      = factor(c("50%", "90%"), levels = c("50%", "90%")),
      type       = "worst-case"
    )
  ) %>%
    mutate(
      # y = x represents 50% error not 100%
      intercept = log10(error.rate / (1 - error.rate)),
      x1 = pmax(barnyard.limits[1], barnyard.limits[1] - intercept),
      x2 = pmin(barnyard.limits[2], barnyard.limits[2] - intercept),
      y1 = x1 + intercept,
      y2 = x2 + intercept
    ) %>%
    mutate(
      x1 = pmax(x1, barnyard.limits[1]),
      x2 = pmin(x2, barnyard.limits[2]),
      y1 = pmax(y1, barnyard.limits[1]),
      y2 = pmin(y2, barnyard.limits[2])
    )
  
  ggplot(data, aes(x = x.reads, y = .data[[y.col]])) +
    rasterize(
      geom_point(data = data %>% filter(cell.classification != "singlet"), color = "gray32", size = 0.8),
      dpi = 300, dev = "ragg_png"
    ) +
    rasterize(
      geom_point(data = data %>% filter(cell.classification == "singlet"), color = "white", size = 1.25),
      dpi = 300, dev = "ragg_png"
    ) +
    rasterize(
      geom_point(data = data %>% filter(cell.classification == "singlet"), color = "black", size = 0.8),
      dpi = 300, dev = "ragg_png"
    ) +
    geom_segment(
      data = error.lines %>% filter(type == "observed"),
      aes(x = x1, xend = x2, y = y1, yend = y2, linetype = label),
      color = "#2078B4", linewidth = 0.4, inherit.aes = FALSE
    ) +
    geom_segment(
      data = error.lines %>% filter(type == "worst-case"),
      aes(x = x1, xend = x2, y = y1, yend = y2, linetype = label),
      color = "#E2201D", linewidth = 0.4, inherit.aes = FALSE
    ) +
    scale_linetype_manual(
      name   = "Contamination",
      values = c("1%" = "dotted", "5%" = "dashed", "10%" = "longdash", "50%" = "solid", "90%" = "solid")
    ) +
    annotation_logticks(
      sides = "bl", outside = TRUE,
      short = unit(4, "pt"), mid = unit(6, "pt"), long = unit(4, "pt"),
      linewidth = 0.2
    ) +
    barnyard.scale.x +
    barnyard.scale.y +
    scale_color_brewer(
      palette = "Paired",
      guide = "none"
    ) +
    coord_fixed(ratio = 1, clip = "off") +
    labs(title = title, x = x.lab, y = y.lab, color = NULL) +
    theme_thesis() +
    theme(
      legend.position.inside = c(0.02, 1),
      legend.justification   = c(0, 1),
      legend.position        = "inside",
      legend.background      = element_blank(),
      title = element_text(hjust = 1)
    )
}

# ── Error fraction plot ───────────────────────────────────────────────────────
# additive component?

plot.error.fraction <- function(data, title) {
  singlets <- data %>% filter(cell.classification == "singlet")
  
  ggplot(data, aes(x = total, y = remaining.fraction)) +
    rasterize(
      geom_point(data = data %>% filter(cell.classification != "singlet"), color = "gray48", size = 0.8),
      dpi = 300, dev = "ragg_png"
    ) +
    rasterize(
      geom_point(data = singlets, color = "black", size = 1.25),
      dpi = 300, dev = "ragg_png"
    ) +
    rasterize(
      geom_point(data = singlets, aes(color = dominant.species), size = 0.8),
      dpi = 300, dev = "ragg_png"
    ) +
    scale_x_log10() +
    scale_y_continuous(labels = percent_format(accuracy = 0.1), limits = c(0, NA)) +
    scale_color_brewer(palette = "Paired") +
    labs(title = title, x = "Total reads", y = "Error fraction") +
    theme_thesis(base_size = 48) +
    theme(legend.position = "none")
}

# ── Contamination heatmap plot ────────────────────────────────────────────────

plot.contamination.heatmap <- function(g) {
  abbrev <- function(x) gsub("^(\\w)\\w+ ", "\\1. ", as.character(x))
  
  ggplot(
    filter(qc.contamination, group == g),
    aes(x = contaminant.species, y = dominant.species, fill = contamination_fraction)
  ) +
    geom_tile() +
    geom_text(aes(label = percent(contamination_fraction, accuracy = 0.1)), size = 2.5) +
    scale_fill_distiller(
      palette   = "Spectral",
      direction = -1,
      limits    = c(0, 0.5),
      oob       = scales::squish,
      labels    = percent_format(),
      guide     = guide_colorbar(
        title          = "Fraction of non-dominant reads",
        direction      = "horizontal",
        title.position = "top",
        barwidth       = 15,
        barheight      = 0.5,
      )
    ) +
    scale_x_discrete(position = "top", labels = abbrev) +
    scale_y_discrete(labels = abbrev) +
    theme_thesis(base_size = 48) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 0, face = "italic"),
      axis.text.y     = element_text(face = "italic"),
      legend.position = "none"
    ) +
    labs(x = "Contaminating species", y = "Dominant species")
}

# ── Save figures ──────────────────────────────────────────────────────────────

p.barnyard.mda.sub <- plot.barnyard(filter(qc.barnyard.log, group == "MDA"), "y.sub", "MDA", "Dominant species reads (log10)", "Subdominant species\nreads (log10)")
p.barnyard.pta.sub <- plot.barnyard(filter(qc.barnyard.log, group == "PTA"), "y.sub", "PTA", "", "")

p.barnyard.mda.rem <- plot.barnyard(filter(qc.barnyard.log, group == "MDA"), "y.rem", "MDA", "Dominant species reads (log10)", "Non-Dominant species\nreads (log10)")
p.barnyard.pta.rem <- plot.barnyard(filter(qc.barnyard.log, group == "PTA"), "y.rem", "PTA", "", "")

p.error.mda <- plot.error.fraction(filter(qc.barnyard, group == "MDA"), "MDA")
p.error.pta <- plot.error.fraction(filter(qc.barnyard, group == "PTA"), "PTA")

p.contam.mda <- plot.contamination.heatmap("MDA") 
p.contam.pta <- plot.contamination.heatmap("PTA") & theme(
  axis.text.y     = element_blank()
)

ggsave("contamination_barnyard.pdf",
       p.barnyard.mda.sub | p.barnyard.pta.sub | p.barnyard.mda.rem | p.barnyard.pta.rem,
       width = 14, height = 8
)
ggsave("contamination_errfraction.pdf",     p.error.mda    / p.error.pta,    width = 8,  height = 8)
ggsave("contamination_heatmap.pdf",         (p.contam.mda | p.contam.pta) + plot_layout(axes = "collect"),   width = 14, height = 6)
