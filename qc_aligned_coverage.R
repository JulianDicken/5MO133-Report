library(duckdb)
library(DBI)
library(dplyr)
library(ggplot2)

# ── Constants ─────────────────────────────────────────────────────────────────

CSV_OPTS <- "delim='\\t', header=true, max_line_size=100000000, strict_mode=false, quote='', auto_detect=true"
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

con <- dbConnect(duckdb())

dbExecute(con, sprintf("
  CREATE VIEW qc AS
  SELECT *,
    CASE
      WHEN reference = 'phix'     THEN 'phix'
      WHEN reference LIKE 'NZ_%%' THEN 'NZ_CDNB'
      ELSE regexp_replace(reference, '_contig_\\d+$', '')
    END AS organism_id
  FROM (
    SELECT reference, sizeof_reference, union_aligned_bases, sumof_aligned_bases, 'MDA' AS group
    FROM read_csv('./data/mda/qc_coverage.tsv', %s)
    UNION ALL
    SELECT reference, sizeof_reference, union_aligned_bases, sumof_aligned_bases, 'PTA' AS group
    FROM read_csv('./data/pta/qc_coverage.tsv', %s)
  )
", CSV_OPTS, CSV_OPTS))

qc.coverage <- dbGetQuery(con, "
  SELECT reference, sizeof_reference, union_aligned_bases, sumof_aligned_bases, organism_id, \"group\"
  FROM qc
") %>%
  mutate(
    species = case_when(
      reference == "phix"          ~ unname(SPECIES_MAP["phix"]),
      str_starts(reference, "NZ_") ~ unname(SPECIES_MAP["NZ_CDNB"]),
      .default = unname(SPECIES_MAP[str_remove(reference, "_contig_\\d+$")])
    )
  )

dbDisconnect(con)

# ── Data processing ───────────────────────────────────────────────────────────

qc.coverage <- qc.coverage %>%
  filter(!(species %in% c("PhiX", "Xanthomonas campestris"))) %>%
  group_by(group, species) %>%
  mutate(
    is.short = sizeof_reference < 80e3,
    coverage = (union_aligned_bases / sizeof_reference) * 100,
    depth    = sumof_aligned_bases / sizeof_reference
  )

qc.coverage.gen <- qc.coverage %>%
  group_by(group, species) %>%
  summarise(gen.size = sum(sizeof_reference), gen.seq = sum(sumof_aligned_bases)) %>%
  select(group, species, gen.size, gen.seq)
  
qc.depth <- qc.coverage %>%
  left_join(qc.coverage.gen, join_by(group, species)) %>%
  mutate(gen.depth = gen.seq / gen.size)

# ── Plot ──────────────────────────────────────────────────────────────────────

coverage.plot <- ggplot(qc.depth) +
  rasterize(
    geom_point(data = qc.depth %>% filter(is.short), aes(x = coverage, y = depth, color = species), size = 0.25, adjust = 4),
    dpi = 300, dev = "ragg_png"
  ) +
  rasterize(
    geom_point(data = qc.depth %>% filter(!is.short), aes(x = coverage, y = depth), size = 0.8, color = "black"),
    dpi = 300, dev = "ragg_png"
  ) +
  rasterize(
    geom_point(data = qc.depth %>% filter(!is.short), aes(x = coverage, y = depth, color = species), size = 0.25, adjust = 4),
    dpi = 300, dev = "ragg_png"
  ) +
  # scale_x_log10() +
  scale_y_log10() +
  scale_color_brewer(
    palette = "Paired"
    # trans = "sqrt"
  ) +
  facet_grid(. ~ group) +
  theme_thesis() +
  theme(
    plot.title    = element_text(hjust = 1),
    plot.subtitle = element_text(hjust = 1)
  )

ggsave("coverage.pdf", coverage.plot, width = 12, height = 4)