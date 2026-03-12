library(Seurat)
library(Signac)
library(stringr)
library(ggplot2)
library(Matrix)
library(GenomicRanges)
library(patchwork)
library(Zorn)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)
library(ggrepel)
library(ggpointdensity)
library(cowplot)
library(ggrastr)

source("./style.R")

plan("multicore")
options(future.globals.maxSize = 1024 * 1024 * 1024 * 64)  # 64 GB
options(future.seed = TRUE)

# ── Constants ─────────────────────────────────────────────────────────────────

MULTIPLET_MIN_FRACTION   <- 0.05
RANDOM_ENTROPY_THRESHOLD <- 0.5
RANDOM_MIN_READS         <- 1e3
RDATA <- "analysis_cached.RData"

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

# ── Species colors ────────────────────────────────────────────────────────────

species.order <- tibble(species = unname(SPECIES_MAP)) %>%
  filter(species != "PhiX") %>%
  mutate(not.atcc = species %in% c("Xanthomonas campestris", "PhiX")) %>%
  arrange(not.atcc, species) %>%
  pull(species)

species.colors <- brewer.pal(length(species.order), "Paired")
names(species.colors) <- species.order

# ── Load & annotate composition ───────────────────────────────────────────────

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
  )

qc.composition.species <- qc.composition %>%
  group_by(group, id, species) %>%
  summarise(countof_reads = sum(countof_reads), .groups = "drop")

# ── QC: multiplet detection ──────────────────────────────────────────────────

qc.multiplet <- qc.composition.species %>%
  group_by(group, id) %>%
  mutate(species.fraction = countof_reads / sum(countof_reads)) %>%
  summarise(
    n.contributors = sum(species.fraction > MULTIPLET_MIN_FRACTION),
    .groups = "drop"
  ) %>%
  mutate(is.multiplet = n.contributors >= 2)

# ── QC: entropy / background detection ───────────────────────────────────────

qc.entropy <- qc.composition.species %>%
  group_by(group, id) %>%
  mutate(species.fraction = countof_reads / sum(countof_reads)) %>%
  summarise(
    entropy         = -sum(ifelse(species.fraction > 0, species.fraction * log2(species.fraction), 0)),
    shannon.entropy = entropy / log2(n()),
    total           = sum(countof_reads),
    .groups = "drop"
  ) %>%
  mutate(is.random = shannon.entropy > RANDOM_ENTROPY_THRESHOLD | total < RANDOM_MIN_READS)

# ── QC: dominant species per cell ─────────────────────────────────────────────

qc.dominant <- qc.composition.species %>%
  dplyr::group_by(group, id) %>%
  dplyr::mutate(
    countof.reads.total = sum(countof_reads),
    species.fraction    = countof_reads / countof.reads.total
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, id, dplyr::desc(species.fraction)) %>%
  dplyr::distinct(group, id, .keep_all = TRUE) %>%
  dplyr::rename(
    dominant.species          = species,
    dominant.species.fraction = species.fraction
  ) %>%
  dplyr::select(-countof_reads)

# ── QC: high-quality singlet filter ──────────────────────────────────────────

qc.high.quality <- qc.multiplet %>%
  left_join(qc.entropy,  by = c("group", "id")) %>%
  left_join(qc.dominant, by = c("group", "id")) %>%
  filter(!is.multiplet, !is.random) %>%
  select(group, id)

# ── Species fraction matrix per cell ─────────────────────────────────────────

qc.fractions <- qc.composition.species %>%
  group_by(group, id) %>%
  mutate(species.fraction = countof_reads / sum(countof_reads)) %>%
  ungroup() %>%
  select(group, id, species, species.fraction) %>%
  pivot_wider(names_from = species, values_from = species.fraction, values_fill = 0)

# ══════════════════════════════════════════════════════════════════════════════
# UMAP COMPUTATION
# ══════════════════════════════════════════════════════════════════════════════

# ── Load + annotate ground truth ─────────────────────────────────────────────

load.adata <- function(g) {
  bascet.root <- paste0("./data/", tolower(g))
  adata <- BascetLoadCountSketchMatrix(bascet.root)
  adata$log10.celldepth <- log10(adata$celldepth)
  
  dominant.map <- qc.dominant %>%
    filter(group == g) %>%
    mutate(id = str_remove(id, "BASCET_"))
  
  adata$ground.truth <- dominant.map$dominant.species[
    match(colnames(adata), dominant.map$id)
  ]
  
  adata
}

# ── Filter to high-quality singlets ──────────────────────────────────────────

filter.adata <- function(adata, g) {
  hq.ids <- qc.high.quality %>%
    filter(group == g) %>%
    mutate(id = str_remove(id, "BASCET_")) %>%
    pull(id)
  hq.ids <- intersect(hq.ids, colnames(adata))
  adata[, hq.ids]
}

# ── Unfiltered kmer UMAP ────────────────────────────────────────────────────

compute.umap.unfiltered <- function(adata) {
  n.dims <- ncol(adata@reductions[["kmersketch"]]@cell.embeddings)
  adata <- RunUMAP(adata, dims = 1:n.dims, reduction = "kmersketch",
                   metric = "cosine", reduction.name = "umap.kmer",
                   reduction.key = "UMAPkmer_")
  adata <- FindNeighbors(adata, reduction = "kmersketch", dims = 1:n.dims)
  adata <- FindClusters(adata, resolution = 0.8)
  adata
}

# ── Filtered kmer UMAP + clusters ────────────────────────────────────────────

compute.umap.filtered <- function(adata) {
  n.dims <- ncol(adata@reductions[["kmersketch"]]@cell.embeddings)
  adata <- RunUMAP(adata, dims = 1:n.dims, reduction = "kmersketch",
                   metric = "cosine", reduction.name = "umap.kmer",
                   reduction.key = "UMAPkmer_")
  adata <- FindNeighbors(adata, reduction = "kmersketch", dims = 1:n.dims)
  adata <- FindClusters(adata, resolution = 0.8)
  adata
}

# ── Species fraction UMAP ───────────────────────────────────────────────────

compute.umap.species.fraction <- function(adata, g) {
  fractions <- qc.fractions %>%
    filter(group == g) %>%
    mutate(id = str_remove(id, "BASCET_")) %>%
    filter(id %in% colnames(adata))
  
  frac.mat <- fractions %>%
    select(-group, -id) %>%
    as.matrix()
  rownames(frac.mat) <- fractions$id
  frac.mat <- frac.mat[colnames(adata), , drop = FALSE]
  colnames(frac.mat) <- paste0("SF_", seq_len(ncol(frac.mat)))
  
  adata[["specfrac"]] <- CreateDimReducObject(
    embeddings = frac.mat,
    key        = "SF_",
    assay      = DefaultAssay(adata)
  )
  
  RunUMAP(adata, reduction = "specfrac",
          dims = 1:ncol(frac.mat), metric = "cosine",
          reduction.name = "umap.species.fraction",
          reduction.key = "UMAPfrac_")
}

# ── Run all ──────────────────────────────────────────────────────────────────

run.analysis <- function(g) {
  adata.all  <- load.adata(g)
  adata.all  <- compute.umap.unfiltered(adata.all)
  adata.filt <- filter.adata(adata.all, g)
  adata.filt <- compute.umap.filtered(adata.filt)
  adata.filt <- compute.umap.species.fraction(adata.filt, g)
  list(all = adata.all, filt = adata.filt)
}

if (!interactive()) {
  adata.mda <- run.analysis("MDA")
  adata.pta <- run.analysis("PTA")
  save.image(RDATA)
} else {
  if (file.exists(RDATA)) {
    message("Loading cached environment ", RDATA)
    load(RDATA)
  } else {
    message("No cache found, running analysis...")
    adata.mda <- run.analysis("MDA")
    adata.pta <- run.analysis("PTA")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# PLOT STYLING
# ══════════════════════════════════════════════════════════════════════════════

centroid.labels <- function(df, group.col) {
  df %>%
    group_by(.data[[group.col]]) %>%
    summarise(umap.1 = median(umap.1), umap.2 = median(umap.2), .groups = "drop")
}

style.umap.base <- function(p) {
  p +
    theme_thesis() +
    theme(
      axis.ticks      = element_blank(),
      axis.text       = element_blank(),
      axis.title      = element_blank(),
      axis.line       = element_blank(),
      panel.border    = element_blank(),
      plot.title      = element_text(size = 16, hjust = 1),
      plot.subtitle   = element_text(size = 14, hjust = 1),
    )
}

style.umap.species <- function(p) {
  style.umap.base(p) +
    scale_color_manual(values = species.colors) +
    theme(legend.text = element_text(face = "italic"))
}

style.umap.clusters <- function(p) {
  n <- nlevels(p$data$cluster)
  style.umap.base(p) +
    scale_color_manual(values = colorRampPalette(brewer.pal(8, "Set2"))(n))
}

# ── Individual UMAP plots ─────────────────────────────────────────────────────

plot.umap.unfiltered <- function(adata, g) {
  coords <- adata@reductions[["umap.kmer"]]@cell.embeddings
  n <- ncol(adata)
  df <- data.frame(umap.1  = coords[, 1], umap.2 = coords[, 2],
                   cluster = adata$seurat_clusters)
  labels <- centroid.labels(df, "cluster")
  
  p <- ggplot(df, aes(x = umap.1, y = umap.2, color = cluster)) +
    rasterize(list(
      geom_point(size = 1.25, color = "black"),
      geom_point(size = 0.8)
    ), dpi = 300, dev = "ragg_png") +
    geom_text_repel(
      data = labels, aes(label = cluster),
      color = "black", size = 2.5,
      bg.color = "white", bg.r = 0.15, show.legend = FALSE
    ) +
    labs(title = "Unfiltered", subtitle = paste0("n = ", n),
         x = "UMAP 1", y = "UMAP 2", color = NULL)
  
  style.umap.clusters(p)
}

plot.umap.species.fraction <- function(adata, g) {
  coords <- adata@reductions[["umap.species.fraction"]]@cell.embeddings
  n <- ncol(adata)
  df <- data.frame(umap.1  = coords[, 1], umap.2 = coords[, 2],
                   species = adata$ground.truth)
  labels <- centroid.labels(df, "species")
  
  p <- ggplot(df, aes(x = umap.1, y = umap.2, color = species)) +
    rasterize(list(
      geom_point(size = 1.25, color = "black"),
      geom_point(size = 0.8)
    ), dpi = 300, dev = "ragg_png") +
    geom_text_repel(
      data = labels, aes(label = species),
      color = "black", size = 2.5, fontface = "italic",
      bg.color = "white", bg.r = 0.15, show.legend = FALSE
    ) +
    labs(title = "Alignment", subtitle = paste0("n = ", n),
         x = "UMAP 1", y = "UMAP 2", color = NULL)
  
  style.umap.species(p)
}

plot.umap.clusters <- function(adata, g) {
  coords <- adata@reductions[["umap.kmer"]]@cell.embeddings
  n <- ncol(adata)
  df <- data.frame(umap.1  = coords[, 1], umap.2 = coords[, 2],
                   cluster = adata$seurat_clusters)
  labels <- centroid.labels(df, "cluster")
  
  p <- ggplot(df, aes(x = umap.1, y = umap.2, color = cluster)) +
    rasterize(list(
      geom_point(size = 1.25, color = "black"),
      geom_point(size = 0.8)
    ), dpi = 300, dev = "ragg_png") +
    geom_text_repel(
      data = labels, aes(label = cluster),
      color = "black", size = 2.5,
      bg.color = "white", bg.r = 0.15, show.legend = FALSE
    ) +
    labs(title = "Filtered", subtitle = paste0("n = ", n),
         x = "UMAP 1", y = "UMAP 2", color = NULL)
  
  style.umap.clusters(p)
}

plot.umap.ground.truth <- function(adata, g) {
  coords <- adata@reductions[["umap.kmer"]]@cell.embeddings
  n <- ncol(adata)
  df <- data.frame(umap.1  = coords[, 1], umap.2 = coords[, 2],
                   species = adata$ground.truth)
  labels <- centroid.labels(df, "species")
  
  p <- ggplot(df, aes(x = umap.1, y = umap.2, color = species)) +
    rasterize(list(
      geom_point(size = 1.25, color = "black"),
      geom_point(size = 0.8)
    ), dpi = 300, dev = "ragg_png") +
    geom_text_repel(
      data = labels, aes(label = species),
      color = "black", size = 2.5, fontface = "italic",
      bg.color = "white", bg.r = 0.15, show.legend = FALSE
    ) +
    labs(title = "Ground Truth", subtitle = paste0("n = ", n),
         x = "UMAP 1", y = "UMAP 2", color = NULL)
  
  style.umap.species(p)
}

plot.umap.density <- function(adata, g) {
  coords <- adata@reductions[["umap.kmer"]]@cell.embeddings
  n <- ncol(adata)
  df <- data.frame(umap.1 = coords[, 1], umap.2 = coords[, 2])
  
  p <- ggplot(df, aes(x = umap.1, y = umap.2)) +
    rasterize(list(
      geom_point(size = 1.25, color = "black"),
      geom_pointdensity(size = 0.8)
    ), dpi = 300, dev = "ragg_png") +
    scale_color_distiller(
      palette = "Spectral",
      guide   = guide_colorbar(
        title          = "Point Density",
        title.position = "right",
        label.position = "left",
        barwidth       = 0.5,
        barheight      = 10
      )
    ) +
    labs(title = "Point Density", subtitle = paste0("n = ", n),
         x = "UMAP 1", y = "UMAP 2")
  
  style.umap.base(p) + theme(legend.title = element_text(angle = 270, hjust = 0, vjust = 1))
}

plot.umap.depth <- function(adata, g) {
  coords <- adata@reductions[["umap.kmer"]]@cell.embeddings
  n <- ncol(adata)
  df <- data.frame(
    umap.1 = coords[, 1], umap.2 = coords[, 2], 
    depth = adata$log10.celldepth
  ) %>% arrange(depth)
  
  p <- ggplot(df, aes(x = umap.1, y = umap.2, color = depth)) +
    rasterize(list(
      geom_point(size = 1.25, color = "black"),
      geom_point(size = 0.8)
    ), dpi = 300, dev = "ragg_png") +
    scale_color_distiller(
      palette = "Spectral",
      guide   = guide_colorbar(
        title          = "Read Depth (log10)",
        title.position = "right",
        label.position = "left",
        barwidth       = 0.5,
        barheight      = 10
      )
    ) +
    labs(title = "Read Depth", subtitle = paste0("n = ", n),
         x = "UMAP 1", y = "UMAP 2")
  
  style.umap.base(p) + theme(legend.title = element_text(angle = 270, hjust = 1, vjust = -1))
}

make.group.figure <- function(res, g) {
  p.a <- plot.umap.unfiltered(res$all, g)        + theme(legend.position = "none")
  p.b <- plot.umap.species.fraction(res$filt, g) + theme(legend.position = "none")
  p.c <- plot.umap.clusters(res$filt, g)         + theme(legend.position = "none")
  p.d <- plot.umap.ground.truth(res$filt, g)     + theme(legend.position = "none")
  p.e <- plot.umap.depth(res$filt, g)
  
  p.legend <- wrap_elements(get_legend(p.e))
  p.e      <- p.e + theme(legend.position = "none")
  
  (p.a | p.c | p.d) / (p.b | p.e | p.legend) +
    plot_annotation(
      title = g,
      theme = theme(plot.title = element_text(size = 24, hjust = 1))
    )
}

p.mda <- make.group.figure(adata.mda, "MDA")
p.pta <- make.group.figure(adata.pta, "PTA")

ggsave("umap_mda.pdf", p.mda, width = 12, height = 8)
ggsave("umap_pta.pdf", p.pta, width = 12, height = 8)