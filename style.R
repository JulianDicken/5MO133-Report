library(showtext)

font_add(
  "Verdana",
   regular    = "/usr/share/fonts/truetype/msttcorefonts/Verdana.ttf",
   bold       = "/usr/share/fonts/truetype/msttcorefonts/Verdana_Bold.ttf",
   italic     = "/usr/share/fonts/truetype/msttcorefonts/Verdana_Italic.ttf",
   bolditalic = "/usr/share/fonts/truetype/msttcorefonts/Verdana_Bold_Italic.ttf"
)
showtext_auto()

theme_thesis <- function(base_size = 24, font_scale = 1) {
  title_size <- 24 * font_scale
  text_size <- 20 * font_scale
  axis_text_size <- 14 * font_scale
  axis_title_size <- 18 * font_scale
  legend_title_size <- 14 * font_scale
  legend_text_size <- 14 * font_scale
  
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      panel.background    = element_rect(fill = "#FFFFFF", colour = NA),
      plot.background     = element_rect(fill = "#FFFFFF", colour = NA),
      strip.background    = element_blank(),
      legend.background   = element_rect(fill = "#FFFFFF", colour = NA),
      legend.key          = element_rect(fill = "#FFFFFF", colour = NA),
      legend.title        = element_text(
        size = legend_title_size, colour = "#000000",
        margin = margin(t = 2.5, b = 2.5), hjust = 0),
      legend.text         = element_text(size = legend_text_size, colour = "#000000"),
      legend.spacing.y    = unit(0, "cm"),
      legend.box.margin   = margin(),
      legend.margin       = margin(),
      text                = element_text(colour = "#000000", size = text_size),
      axis.text           = element_text(colour = "#000000", size = axis_text_size),
      axis.title          = element_text(
        colour = "#000000", size = axis_title_size,
        margin = margin(t = 1, r = 1, b = 1, l = 1)),
      plot.title          = element_text(size = title_size, margin = margin(b = 5)),
      plot.title.position = "plot",
      axis.line  = element_line(colour = "#000000", linewidth = base_size / 48),
      axis.ticks = element_line(colour = "#000000", linewidth = base_size / 48),
      panel.border = element_blank()
    )
}