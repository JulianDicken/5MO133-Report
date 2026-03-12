library(ggplot2)
library(dplyr)
library(readr)
library(ggrepel)
library(stringr)
library(ggpointdensity)
library(patchwork)
library(tidyr)
library(forcats)
library(cowplot)
library(duckdb)
library(DBI)
library(ggpubr)
library(ggh4x)
library(legendry)
library(ggsignif)
library(gghalves)
library(ggbeeswarm)

# ── Constants ─────────────────────────────────────────────────────────────────

MULTIPLET_MIN_FRACTION  <- 0.05
RANDOM_ENTROPY_THRESHOLD <- 0.5
RANDOM_MIN_READS        <- 1e3
DOMINANT_GOOD_CUTOFF    <- 0.99
HISTOGRAM_BINS          <- 20
PERCENTILE_THRESHOLD    <- 0.1

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

# ── Load & annotate data ──────────────────────────────────────────────────────

qc.composition <- bind_rows(
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
  left_join(
    qc.coverage %>% ungroup() %>% select(group, reference, sizeof_reference) %>% distinct(),
    join_by(group, reference)
  )
  # filter(sizeof_reference > 80e3 | species == "Xanthomonas campestris")

qc.composition.species <- qc.composition %>%
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
    entropy          = -sum(ifelse(species.fraction > 0, species.fraction * log2(species.fraction), 0)),
    shannon.entropy = entropy / log2(n()),
    total            = sum(countof_reads),
    .groups = "drop"
  ) %>%
  mutate(is.random = shannon.entropy > RANDOM_ENTROPY_THRESHOLD | total < RANDOM_MIN_READS)

# Dominant species per cell
# namespacing here required for some reason??
qc.dominant <- qc.composition.species %>%
  dplyr::group_by(group, id) %>%
  dplyr::mutate(
    countof_reads_total = sum(countof_reads),
    species.fraction    = countof_reads / countof_reads_total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, id, dplyr::desc(species.fraction)) %>%
  dplyr::distinct(group, id, .keep_all = TRUE) %>%
  dplyr::rename(dominant.species          = species,
                dominant.species.fraction = species.fraction) %>%
  dplyr::select(-countof_reads)

species.order <- qc.dominant %>%
  distinct(dominant.species) %>%
  mutate(not.atcc = dominant.species %in% c("NZ_CDNB", "PhiX")) %>%
  arrange(not.atcc, dominant.species) %>%
  pull(dominant.species)

qc.dominant <- qc.dominant %>%
  mutate(dominant.species = factor(dominant.species, levels = species.order))

qc.high_quality <- qc.multiplet %>%
  left_join(qc.entropy,  by = c("group", "id")) %>%
  left_join(qc.dominant, by = c("group", "id")) %>%
  filter(!is.multiplet, !is.random) %>%
  select(group, id)

# ── Histogram plot ────────────────────────────────────────────────────────────

hist.good.pos <- DOMINANT_GOOD_CUTOFF + 1 * (1/HISTOGRAM_BINS)
hist.shared.breaks <- c(seq(0, DOMINANT_GOOD_CUTOFF, 0.1), hist.good.pos)
hist.shared.limits <- c(0, hist.good.pos + 1 / HISTOGRAM_BINS)

qc.hist.bins_all <- c(
  round(seq(0, DOMINANT_GOOD_CUTOFF - 1/HISTOGRAM_BINS, 1 / HISTOGRAM_BINS), 4),
  hist.good.pos
)

qc.hist.data <- qc.dominant %>%
  # semi_join(qc.high_quality, by = c("group", "id")) %>%
  mutate(
    bin = if_else(
      dominant.species.fraction > DOMINANT_GOOD_CUTOFF,
      hist.good.pos,
      round(floor(dominant.species.fraction * HISTOGRAM_BINS) / HISTOGRAM_BINS, 4)
    )
  ) %>%
  group_by(group, bin, dominant.species) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(
    group            = unique(qc.dominant$group),
    bin              = qc.hist.bins_all,
    dominant.species = levels(qc.dominant$dominant.species),
    fill = list(n = 0)
  ) %>%
  mutate(dominant.species = factor(dominant.species, levels = levels(qc.dominant$dominant.species)))

hist.y.scale <- qc.hist.data %>%
  group_by(group, bin) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  pull(n) %>%
  max()

qc.hist.cumulative <- qc.hist.data %>%
  group_by(group, bin) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  arrange(group, bin) %>%
  group_by(group) %>%
  mutate(cumulative = cumsum(n) / sum(n)) %>%
  bind_rows(
    distinct(., group) %>% mutate(bin = 0, n = 0, cumulative = 0)
  ) %>%
  arrange(group, bin)

qc.hist.percentile <- qc.dominant %>%
  arrange(group, dominant.species.fraction) %>%
  group_by(group) %>%
  mutate(n_cumulative = row_number()) %>%
  summarise(
    inflection = dominant.species.fraction[which.min(abs(n_cumulative / n() - PERCENTILE_THRESHOLD))]
  )

qc.hist.plot <- function(g) {
  inflection <- qc.hist.percentile %>% filter(group == g) %>% pull(inflection)
  
  proximity.threshold <- 0.0
  nearby.breaks <- abs(hist.shared.breaks - inflection) < proximity.threshold
  base.breaks <- hist.shared.breaks[!nearby.breaks]
  
  axis.breaks <- sort(c(base.breaks, inflection))
  axis.labels <- ifelse(
    axis.breaks == hist.good.pos,   paste0(">", DOMINANT_GOOD_CUTOFF),
    ifelse(axis.breaks == inflection, round(inflection, 2), as.character(axis.breaks))
  )
  
  ggplot(filter(qc.hist.data, group == g), aes(x = bin, y = n, fill = dominant.species)) +
    geom_bar(stat = "identity", width = 1 / HISTOGRAM_BINS) +
    geom_line(
      data = filter(qc.hist.cumulative, group == g),
      aes(
        x = bin, y = cumulative * hist.y.scale, 
        group = group
      ),
      inherit.aes = FALSE, color = "black", linewidth = 0.8
    ) +
    geom_vline(xintercept = inflection, linetype = "dotted") +
    annotate(
      "text", 
      x = inflection, y = Inf,
      label = paste0("P", ((1 - PERCENTILE_THRESHOLD) * 100)),
      hjust = -0.2, vjust = 1.5
    ) +
    scale_x_continuous(
      breaks = axis.breaks,
      minor_breaks = seq(0.05, DOMINANT_GOOD_CUTOFF, 0.1),
      labels = c()
    ) +
    scale_y_continuous(
      limits = c(0, hist.y.scale),
      name = "Number of cells",
      sec.axis = sec_axis(~ . / hist.y.scale, name = "Cumulative proportion")
    ) +
    scale_fill_brewer(
      palette = "Paired",
      guide   = if (g == "MDA") guide_legend() else "none"
    ) +
    coord_cartesian(xlim = hist.shared.limits) +
    theme_thesis() +
    theme(
      # axis.text.x = element_text(angle = 45, hjust = 1),
      # doesnt work, need to change axis title in svg designer
      # axis.title.y.left = element_text(margin = margin(r = -10000)),\
      axis.title.y.left = element_text(hjust = 0.5, vjust = -20),
      legend.position  = c(0.02, 0.995),
      legend.justification = c(0, 1),
      legend.text = element_text(face = "italic"),
      legend.background = element_blank()
    ) +
    labs(
      title = g, 
      x = "Dominant species fraction",
      fill = NULL
    )
}

qc.hist.heatmap <- qc.hist.data %>%
  group_by(group, bin) %>%
  mutate(
    total          = sum(n),
    fraction       = if_else(total > 0, n / total, 0),
    fraction_log10 = if_else(fraction > 0, fraction, 0)
  ) %>%
  select(-total) %>%
  ungroup()

qc.hist.heatmap.plot <- function(g, axis.breaks, axis.labels) {
  inflection <- qc.hist.percentile %>% filter(group == g) %>% pull(inflection)
  
  proximity.threshold <- 0.0
  nearby.breaks <- abs(hist.shared.breaks - inflection) < proximity.threshold
  axis.breaks <- sort(c(hist.shared.breaks[!nearby.breaks], inflection))
  axis.labels <- ifelse(
    axis.breaks == hist.good.pos,     paste0(">", DOMINANT_GOOD_CUTOFF),
    ifelse(axis.breaks == inflection, round(inflection, 2), as.character(axis.breaks))
  )
  
  ggplot(qc.hist.heatmap %>% filter(group == g),
         aes(
           xmin = bin - 1 / (2 * HISTOGRAM_BINS),
           xmax = bin + 1 / (2 * HISTOGRAM_BINS),
           ymin = as.numeric(dominant.species) - 0.5,
           ymax = as.numeric(dominant.species) + 0.5,
           fill = fraction_log10
         )
  ) +
    geom_rect(interpolate = FALSE) +
    scale_x_continuous(
      breaks = axis.breaks,
      minor_breaks = seq(0.05, DOMINANT_GOOD_CUTOFF, 0.1),
      labels = axis.labels,
    ) +
    scale_y_reverse(
      breaks = seq_along(levels(qc.hist.heatmap$dominant.species)),
      labels = levels(qc.hist.heatmap$dominant.species),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = hist.shared.limits) +
    scale_fill_distiller(
      palette   = "Spectral",
      direction = -1,
      limits    = c(0, 1),
      trans     = "sqrt",
      guide     = guide_colorbar(
        title          = "Relative Abundance in Bin",
        direction      = "horizontal",
        title.position = "top",
        label.position = "bottom",
        barwidth       = 15,
        barheight      = 0.5,
      )
    ) +
    theme_thesis() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(face = "italic")
    ) +
    labs(x = "Dominant species fraction", y = NULL)
}

legend.carrier <- ggplot(data.frame(x = 0, y = 0, fill = 0.5),
                         aes(x = x, y = y, fill = fill)) +
  geom_point(alpha = 0) +
  scale_fill_distiller(
    palette        = "Spectral",
    direction      = -1,
    limits         = c(0, 1),
    trans          = "sqrt",
    guide          = guide_colorbar(
      title          = "Relative Abundance in Bin",
      direction      = "horizontal",
      title.position = "top",
      label.position = "bottom",
      barwidth       = 10,
      barheight      = 0.5
    )
  ) +
  theme_void() +
  theme(legend.position = "bottom")

qc.hist.species.global <- qc.dominant %>%
  group_by(group, dominant.species) %>%
  summarise(n = n(), .groups = "drop")

hist.species.global.plot <- function(g) {
  ggplot(filter(qc.hist.species.global, group == g), aes(x = n, y = dominant.species, fill = dominant.species)) +
    geom_bar(stat = "identity") +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    scale_y_discrete(limits = rev(levels(qc.dominant$dominant.species))) +
    scale_fill_brewer(palette = "Paired", guide = "none") +
    theme_thesis() +
    theme(axis.text.y = element_blank()) +
    labs(x = "Cell count", y = NULL, title = g)
}

hist.groups  <- unique(qc.hist.percentile$group)

hist.plots   <- wrap_plots(lapply(hist.groups, qc.hist.plot), axes = "collect") &
  theme(axis.title.x = element_blank())

heat.plots   <- wrap_plots(lapply(hist.groups, qc.hist.heatmap.plot), axes = "collect") &
  theme(legend.position = "none")

global.plots <- wrap_plots(
  c(lapply(hist.groups, hist.species.global.plot), list(wrap_elements(get_legend(legend.carrier)))),
  ncol    = 1,
  heights = c(1, 1, 0.15)
)

((hist.plots / heat.plots) + plot_layout(heights = c(1.5, 1)) |
    global.plots) +
  plot_layout(widths = c(1, 0.3))


# ── Violin data ───────────────────────────────────────────────────────────────

violin.colors <- c(
  setNames(brewer.pal(length(species.order), "Paired"), species.order),
  "Total"  = "black",
  "spacer" = NA
)

species.levels.violin <- c(
  "Total",
  species.order %>%
    as.character() %>%
    Filter(function(x) !(x %in% c("Xanthomonas campestris", "PhiX")), .)
)

qc.violin <- qc.barnyard %>%
  filter(
    # cell.classification == "singlet",
    !(dominant.species %in% c("Xanthomonas campestris", "PhiX"))
  ) %>%
  mutate(separation.ratio = log10((x.reads + 1) / (y.rem + 1))) %>%
  bind_rows(mutate(., dominant.species = "Total")) %>%
  mutate(
    dominant.species = factor(as.character(dominant.species),
                              levels = species.levels.violin)
  )

violin.stats <- qc.violin %>%
  group_by(dominant.species, group) %>%
  summarise(
    mean   = mean(separation.ratio,   na.rm = TRUE),
    median = median(separation.ratio, na.rm = TRUE),
    ymin   = median(separation.ratio, na.rm = TRUE) - 1.1 * IQR(separation.ratio, na.rm = TRUE),
    ymax   = median(separation.ratio, na.rm = TRUE) + 1.1 * IQR(separation.ratio, na.rm = TRUE),
    .groups = "drop"
  )

violin.annotate <- qc.violin %>%
  group_by(dominant.species) %>%
  summarise(
    p.value = wilcox.test(
      separation.ratio[group == "MDA"],
      separation.ratio[group == "PTA"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p.adj    = p.adjust(p.value, method = "bonferroni"),
    p.signif = case_when(
      p.adj < 0.00001 ~ "*****",
      p.adj < 0.0001  ~ "****",
      p.adj < 0.001   ~ "***",
      p.adj < 0.01    ~ "**",
      p.adj < 0.05    ~ "*",
      TRUE            ~ "ns"
    ),
    dominant.species = factor(as.character(dominant.species),
                              levels = species.levels.violin)
  ) %>%
  arrange(dominant.species) %>%
  left_join(
    violin.stats %>%
      select(dominant.species, group, mean) %>%
      pivot_wider(names_from = group, values_from = mean) %>%
      mutate(fold.change = 10^(PTA - MDA)),
    by = "dominant.species"
  ) %>%
  mutate(fold.label = sprintf("%.2fx", fold.change))

# ── Sort species ───────────────────────────────────────────────

species.order.violin <- violin.annotate %>%
  filter(dominant.species != "Total") %>%
  arrange(p.adj) %>%
  pull(dominant.species) %>%
  as.character()

species.levels.sorted <- c("Total", species.order.violin)

qc.violin <- qc.violin %>%
  mutate(dominant.species = factor(as.character(dominant.species),
                                   levels = species.levels.sorted))

violin.stats <- violin.stats %>%
  mutate(dominant.species = factor(as.character(dominant.species),
                                   levels = species.levels.sorted))

violin.annotate <- violin.annotate %>%
  mutate(dominant.species = factor(as.character(dominant.species),
                                   levels = species.levels.sorted))

n.half <- ceiling(length(species.order.violin) / 2)
species.first.half  <- species.order.violin[1:n.half]
species.second.half <- species.order.violin[(n.half + 1):length(species.order.violin)]

plot.violin.group <- function(species.subset, include.total = TRUE) {
  levels.subset <- c(
    if (include.total) "Total",
    species.subset
  )
  
  data  <- qc.violin %>%
    filter(dominant.species %in% levels.subset) %>%
    mutate(dominant.species = factor(as.character(dominant.species), levels = levels.subset))
  
  stats <- violin.stats %>%
    filter(dominant.species %in% levels.subset) %>%
    mutate(dominant.species = factor(as.character(dominant.species), levels = levels.subset)) %>%
    mutate(color = ifelse(dominant.species == "Total", "outline_white", "outline_black"))
  
  annot <- violin.annotate %>%
    filter(dominant.species %in% levels.subset) %>%
    mutate(dominant.species = factor(as.character(dominant.species), levels = levels.subset))
  
  outline <- data %>% 
    filter(group == "PTA") %>%
    mutate(outline = ifelse(dominant.species == "Total", "outline_white", "outline_black"))
  
  ggplot(data, aes(x = dominant.species, y = separation.ratio, fill = dominant.species)) +
    rasterise(geom_beeswarm(
      data = data %>% filter(group == "MDA"),
      aes(color = dominant.species),
      corral = "wrap", corral.width = 0.4,
      side = -1, shape = 16, size = 0.8
    ), dev = "ragg_png", dpi = 300) +
    rasterise(geom_beeswarm(
      data = outline %>% filter(group == "PTA"),
      aes(color = outline),
      corral = "wrap", corral.width = 0.4,
      side = 1, shape = 16, size = 1.25
    ), dev = "ragg_png", dpi = 300) +
    rasterise(geom_beeswarm(
      data = data %>% filter(group == "PTA"),
      aes(color = dominant.species),
      corral = "wrap", corral.width = 0.4,
      side = 1, shape = 16, size = 0.8
    ), dev = "ragg_png", dpi = 300) +
    geom_errorbar(
      data = stats,
      aes(
        x     = dominant.species,
        ymin  = ymin,
        ymax  = ymax,
        group = interaction(dominant.species, group),
        color = color
      ),
      width       = 0.2,
      linewidth   = 0.4,
      position    =  position_dodge(0.8),
      inherit.aes = FALSE
    ) +
    geom_errorbar(
      data = stats,
      aes(
        x     = dominant.species,
        ymin  = median,
        ymax  = median,
        group = interaction(dominant.species, group),
        color = color
      ),
      width       = 0.4,
      linewidth   = 0.4,
      position    =  position_dodge(0.8),
      inherit.aes = FALSE
    ) +
    geom_point(
      data = stats,
      aes(
        x     = dominant.species,
        y     = mean,
        group = interaction(dominant.species, group)
      ),
      shape       = 22,
      size        = 2,
      fill        = "white",
      color       = "black",
      stroke      = 0.4,
      position    = position_dodge(0.8),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = annot,
      aes(x = as.numeric(dominant.species), y = 4, label = p.signif),
      size        = 5,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = annot,
      aes(x = as.numeric(dominant.species), y = 3.8, label = fold.label),
      size        = 3,
      color       = "grey20",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = data %>%
        distinct(dominant.species) %>%
        mutate(
          x        = as.numeric(dominant.species),
          label    = ifelse(dominant.species == "Total", "Total",
                            gsub("^(\\w)\\w+ ", "\\1. ", as.character(dominant.species))),
          fontface = ifelse(dominant.species == "Total", "plain", "italic")
        ),
      aes(x = x, y = -2, label = label, fontface = fontface),
      size        = 3.5,
      inherit.aes = FALSE
    ) +
    scale_y_continuous(breaks = seq(-2, 4, 1)) +
    scale_color_manual(values = c(violin.colors, "outline_white" = "white", "outline_black" = "black")) +
    scale_fill_manual(values = violin.colors) +
    coord_cartesian(ylim = c(-2, 4), clip = "off") +
    annotation_logticks(
      sides     = "l",
      outside   = TRUE,
      short     = unit(8, "pt"),
      mid       = unit(12, "pt"),
      long      = unit(18, "pt"),
      linewidth = 0.5
    ) +
    theme_thesis() +
    theme(
      legend.position = "none",
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.text.y     = element_text(margin = margin(r = 20)),
      panel.border    = element_blank()
    ) +
    labs(y = "Dominant / Non-Dominant\nReads (log10)", x = NULL)
}

ggsave("violin_separation_upper.pdf", plot.violin.group(species.first.half, include.total = TRUE), width = 8, height = 6)
ggsave("violin_separation_lower.pdf", plot.violin.group(species.second.half, include.total = TRUE), width = 8, height = 6)