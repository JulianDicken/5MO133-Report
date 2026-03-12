library(duckdb)
library(DBI)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)
library(patchwork)
library(ggh4x)
library(RColorBrewer)
library(farver)
library(ggtext)
library(ggrastr)
library(tibble)
library(ggpattern)

# ── Constants ─────────────────────────────────────────────────────────────────

# BIN_SIZE <- 1000L
BIN_SIZE <- 10000L

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

CSV_OPTS <- paste(
  "delim='\\t'",
  "header=true",
  "max_line_size=100000000",
  "strict_mode=false",
  "quote=''",
  "auto_detect=true",
  sep = ", "
)

# ── DuckDB setup ─────────────────────────────────────────────────────────────

con <- dbConnect(duckdb())

dbExecute(con, sprintf("
  CREATE VIEW qc_raw AS
  SELECT *,
    CASE
      WHEN reference = 'phix' THEN 'phix'
      WHEN reference LIKE 'NZ_%%' THEN 'NZ_CDNB'
      ELSE regexp_replace(reference, '_contig_\\d+$', '')
    END AS organism_id
  FROM (
    SELECT *, 'MDA' AS \"group\" FROM read_csv('./data/mda/qc_coverage.tsv', %s)
    UNION ALL
    SELECT *, 'PTA' AS \"group\" FROM read_csv('./data/pta/qc_coverage.tsv', %s)
  )
", CSV_OPTS, CSV_OPTS))

# ── Contig offsets ────────────────────────────────────────────────────────────

dbExecute(con, "
  CREATE TABLE contig_offsets AS
  WITH distinct_contigs AS (
    SELECT DISTINCT organism_id, reference, sizeof_reference
    FROM qc_raw
  ),
  ordered AS (
    SELECT *,
      ROW_NUMBER() OVER (PARTITION BY organism_id ORDER BY reference) AS contig_rank
    FROM distinct_contigs
  )
  SELECT *,
    COALESCE(
      SUM(sizeof_reference) OVER (
        PARTITION BY organism_id
        ORDER BY contig_rank
        ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
      ),
      0
    ) AS contig_offset
  FROM ordered
")

# ── RLE decode -> genome bins ───────────────────────────────────────

dbExecute(con, sprintf("
  CREATE TABLE pileup_binned AS
  WITH genome_sizes AS (
    SELECT organism_id, SUM(sizeof_reference) AS genome_size
    FROM contig_offsets
    GROUP BY organism_id
  ),
  contig_bins AS (
    SELECT c.organism_id, c.reference, c.sizeof_reference, c.contig_offset,
      g.genome_size,
      GREATEST(CEIL(c.sizeof_reference::DOUBLE / %d), 1)::BIGINT AS n_contig_bins
    FROM contig_offsets c
    JOIN genome_sizes g USING (organism_id)
  ),
  runs AS (
    SELECT
      q.\"group\", q.id, q.organism_id,
      cb.contig_offset, cb.sizeof_reference,
      cb.genome_size, cb.n_contig_bins,
      UNNEST(string_split(q.pileup, ',')) AS rle_token
    FROM qc_raw q
    JOIN contig_bins cb USING (organism_id, reference)
    WHERE q.pileup IS NOT NULL
      AND q.pileup != ''
  ),
  parsed AS (
    SELECT
      \"group\", id, organism_id, genome_size,
      contig_offset, sizeof_reference, n_contig_bins,
      CAST(string_split(rle_token, ':')[1] AS BIGINT) AS local_start,
      CAST(string_split(rle_token, ':')[2] AS DOUBLE)  AS run_depth,
      CAST(string_split(rle_token, ':')[3] AS BIGINT)  AS run_length
    FROM runs
    WHERE rle_token != ''
  ),
  with_ends AS (
    SELECT *,
      local_start + run_length AS local_end,
      LEAST(CAST(FLOOR(local_start::DOUBLE / sizeof_reference * n_contig_bins) AS BIGINT), n_contig_bins - 1) AS first_bin,
      LEAST(CAST(FLOOR((local_start + run_length - 1)::DOUBLE / sizeof_reference * n_contig_bins) AS BIGINT), n_contig_bins - 1) AS last_bin
    FROM parsed
    WHERE run_length > 0
  ),
  expanded AS (
    SELECT
      \"group\", id, organism_id, genome_size,
      contig_offset, sizeof_reference, n_contig_bins, run_depth,
      local_start, local_end,
      UNNEST(generate_series(first_bin, last_bin)) AS local_bin
    FROM with_ends
  )
  SELECT
    \"group\", id, organism_id,
    contig_offset + (local_bin + 0.5) * sizeof_reference / n_contig_bins AS position,
    SUM(
      run_depth * (
        LEAST(local_end, (local_bin + 1)::DOUBLE * sizeof_reference / n_contig_bins)
        - GREATEST(local_start, local_bin::DOUBLE * sizeof_reference / n_contig_bins)
      ) / (sizeof_reference::DOUBLE / n_contig_bins)
    ) AS weighted_depth
  FROM expanded
  GROUP BY \"group\", id, organism_id, contig_offset, sizeof_reference, n_contig_bins, local_bin, genome_size
", BIN_SIZE))

# ── Retrieve pileups ──────────────────────────────────────────────────────────

pileup.raw <- dbGetQuery(con, "SELECT * FROM pileup_binned")

contigs <- dbGetQuery(con, "
  SELECT organism_id, reference, contig_offset, sizeof_reference
  FROM contig_offsets
  ORDER BY organism_id, contig_offset
")

genome.sizes <- dbGetQuery(con, "
  SELECT organism_id, SUM(sizeof_reference) AS genome_size
  FROM contig_offsets
  GROUP BY organism_id
")

dbDisconnect(con)

# ── Annotate ──────────────────────────────────────────────────────────

pileup <- pileup.raw %>%
  mutate(species = SPECIES_MAP[organism_id])

contigs <- contigs %>%
  left_join(genome.sizes, by = "organism_id") %>%
  mutate(
    species  = SPECIES_MAP[organism_id],
    bp.start = contig_offset,
    bp.end   = contig_offset + sizeof_reference
  )

# ── Cell ordering ─────────────────────────────────────────────

pileup <- pileup %>%
  group_by(group, id, species) %>%
  mutate(
    mean.depth  = mean(weighted_depth),
    total.depth = sum(weighted_depth)
  ) %>%
  ungroup()

pileup.ranks <- pileup %>%
  distinct(group, id, species, mean.depth, total.depth) %>%
  arrange(group, species, total.depth) %>%
  group_by(group, species) %>%
  mutate(rank = row_number()) %>%
  select(group, id, species, rank) %>%
  ungroup()

pileup <- pileup %>%
  left_join(pileup.ranks, by = join_by(group, id, species))

# ── rescale to per-group per-species  ────────────────────────────────────────────

pileup <- pileup %>%
  group_by(group, species) %>%
  mutate(
    # normalised value used in most places so scaling +1e-9 doesnt result in negative values
    log.depth  = log10(weighted_depth + 1e-9),
    norm.depth = (log.depth - min(log.depth)) / (max(log.depth) - min(log.depth))
  ) %>%
  ungroup()

# ── color fill to match species in earlier plots ─────────────────────────────────────────────────────────

species.order <- tibble(species = unname(SPECIES_MAP)) %>%
  filter(species != "PhiX") %>%
  mutate(not.atcc = species %in% c("PhiX")) %>%
  arrange(not.atcc, species) %>%
  pull(species)

species.colors <- brewer.pal(length(species.order), "Paired")
names(species.colors) <- species.order

species.oklab <- convert_colour(t(col2rgb(species.colors)), from = "rgb", to = "oklab")
rownames(species.oklab) <- names(species.colors)

pileup <- pileup %>%
  filter(species != "PhiX") %>%
  mutate(
    fill = {
      idx <- match(species, rownames(species.oklab))
      ok.l <- species.oklab[idx, "l"]
      ok.a <- species.oklab[idx, "a"]
      ok.b <- species.oklab[idx, "b"]
      
      mid <- 0.5
      t.low  <- pmin(norm.depth / mid, 1)
      t.high <- pmax((norm.depth - mid) / (1 - mid), 0)
      
      new.l <- ifelse(norm.depth <= mid,
                      1.0 - t.low * (1.0 - ok.l),
                      ok.l - t.high * (ok.l - 0.15))
      new.a <- ifelse(norm.depth <= mid,
                      t.low * ok.a,
                      ok.a * (1 + t.high * 0.8))
      new.b <- ifelse(norm.depth <= mid,
                      t.low * ok.b,
                      ok.b * (1 + t.high * 0.8))
      
      new.rgb <- convert_colour(cbind(l = new.l, a = new.a, b = new.b),
                                from = "oklab", to = "rgb")
      new.rgb <- pmin(pmax(new.rgb, 0), 255)
      rgb(new.rgb[, 1], new.rgb[, 2], new.rgb[, 3], maxColorValue = 255)
    }
  )

bottom.limits <- pileup %>%
  group_by(group, species, position) %>%
  summarise(total.depth = sum(weighted_depth), .groups = "drop") %>%
  mutate(total.depth.t = log10(total.depth + 1)) %>%
  group_by(species) %>%
  summarise(ymax = max(total.depth.t)) %>%
  deframe()

side.limits <- pileup %>%
  distinct(group, id, species, total.depth) %>%
  mutate(total.depth.t = log10(total.depth + 1)) %>%
  group_by(species) %>%
  summarise(xmax = max(total.depth.t)) %>%
  deframe()

# ── plot functions ───────────────────────────────────────────────

pileup.species.heatmap <- function(g, sp) {
  df.heat <- pileup %>%
    filter(species == sp & group == g)
  
  ggplot(df.heat, aes(x = position, y = rank, fill = fill)) +
    rasterise(geom_tile(width = BIN_SIZE), dpi = 300, dev = "ragg_png") +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = scales::breaks_width(500000),
      labels = function(x) paste0(x / 1e6, " Mbp"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(breaks = c(), expand = c(0, 0)) +
    labs(x = "", y = "") +
    theme_thesis() +
    theme(
      axis.line = element_blank(),
      axis.line.x.bottom = element_line(colour = "#000000", linewidth = 0.25),
      plot.margin = margin(l = 10, r = 10, t = 10, b = 10)
    )
}

pileup.species.bottom.hist <- function(g, sp) {
  n.cells <- pileup %>%
    filter(species == sp & group == g) %>%
    distinct(id) %>%
    nrow()
  
  df.bhist <- pileup %>%
    filter(species == sp & group == g) %>%
    group_by(position) %>%
    summarise(total.depth = sum(weighted_depth), .groups = "drop") %>%
    mutate(total.depth.t = log10(total.depth + 1))
  
  ggplot(df.bhist, aes(x = position, y = total.depth.t)) +
    rasterise(geom_col(width = BIN_SIZE, fill = species.colors[sp], position = "identity"), dpi = 300, dev = "ragg_png") +
    scale_x_continuous(
      breaks = scales::breaks_width(500000),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      position = "right",
      expand = c(0, 0),
      limits = c(0, bottom.limits[sp])
    ) +
    annotation_logticks(
      sides = "r", outside = TRUE,
      short = unit(1, "pt"), mid = unit(2, "pt"), long = unit(3, "pt"),
      linewidth = 0.2
    ) +
    coord_cartesian(clip = "off") +
    labs(title = paste0(g, " (n = ", n.cells, ")"), x = "", y = "") +
    theme_thesis() +
    theme(
      plot.title = element_text(size = 24),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
}

pileup.species.side.hist <- function(g, sp) {
  df.shist <- pileup %>%
    filter(species == sp & group == g) %>%
    distinct(id, rank, total.depth) %>%
    mutate(total.depth.t = log10(total.depth + 1))
  
  ggplot(df.shist, aes(x = total.depth.t, y = rank)) +
    rasterise(geom_segment(aes(x = 0, xend = total.depth.t, yend = rank), color = species.colors[sp], linewidth = 0.3), dpi = 300, dev = "ragg_png") +
    scale_x_continuous(
      expand = c(0, 0),
      position = "top",
      limits = c(0, side.limits[sp]),
    ) +
    scale_y_continuous(breaks = c(), expand = c(0, 0)) +
    annotation_logticks(
      sides = "t", outside = TRUE,
      short = unit(1, "pt"), mid = unit(2, "pt"), long = unit(3, "pt"),
      linewidth = 0.2
    ) +
    coord_cartesian(clip = "off") +
    labs(x = expression(log[10](depth)), y = "") +
    theme_thesis() +
    theme(
      axis.ticks.x = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(0, 0, 0, 0)
    )
}

# ── Render ───────────────────────────────────────────────────────────────────

pileup.species.block <- function(g, sp) {
  p.bhist <- pileup.species.bottom.hist(g, sp)
  p.heat  <- pileup.species.heatmap(g, sp)
  p.shist <- pileup.species.side.hist(g, sp)
  
  design <- "
    AB
    CD
  "
  
  wrap_plots(
    A = p.bhist, 
    B = plot_spacer(), 
    C = p.heat, 
    D = p.shist,
    design = design, 
    widths = c(4, 1), heights = c(1, 3)
  )
}

species.plot <- setdiff(unname(SPECIES_MAP), c("PhiX"))

for (sp in species.plot) {
  p <- (pileup.species.block("MDA", sp) / pileup.species.block("PTA", sp)) +
    plot_annotation(
      title = paste0("*", sp, "*"),
      theme = theme(plot.title = element_markdown(size = 48))
    ) &
    theme(panel.border = element_blank())
  
  sp.filename <- gsub(" ", "_", tolower(sp))
  ggsave(paste0("pileup_", sp.filename, ".pdf"), p, width = 12, height = 14)
}

# ──  Global dist ────────────────────────── 

df.dist <- pileup %>%
  filter(species %in% "PhiX") %>%
  distinct(group, id, species, total.depth) %>%
  group_by(group, species) %>%
  summarise(total.depth = sum(total.depth), .groups = "drop") %>%
  mutate(
    total.depth.t = log10(total.depth),
    species       = factor(species, levels = rev(species.order)),
    group         = factor(group, levels = c("PTA", "MDA"))
  )

p.dist <- ggplot(df.dist, aes(x = total.depth.t, y = species, fill = species, pattern = group)) +
  rasterise(geom_col_pattern(
    aes(group = interaction(species, group)),
    position          = position_dodge(width = 0.8),
    width             = 0.7,
    pattern_fill      = NA,
    pattern_colour    = "black",
    pattern_density   = 0.1,
    pattern_spacing   = 0.01,
    pattern_size      = 0.5,
    pattern_angle     = 45
  ), dev = "ragg_png", dpi = 300) +
  scale_fill_manual(values = species.colors, guide = "none") +
  scale_pattern_manual(values = c("MDA" = "crosshatch", "PTA" = "none"), name = NULL, guide = "none") +
  scale_x_continuous(name = "Total Reads (log10)", expand = c(0, NA)) +
  scale_y_discrete(name = NULL) +
  annotation_logticks(
    sides = "b", outside = TRUE,
    short = unit(6, "pt"), mid = unit(16, "pt"), long = unit(20, "pt"),
    linewidth = 1.8
  ) +
  coord_cartesian(clip = "off") +
  theme_thesis(base_size = 96, font_scale = 4) +
  theme(
    axis.text.y     = element_blank(),
    panel.border    = element_blank(),
  )

ggsave("pileup_global_read_distribution.pdf", p.dist, width = 12, height = 14)