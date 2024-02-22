library(ggplot2)
library(ggrepel)
library(showtext)
library(ggdist)
library(ggh4x)

font_add(family='Barlow Semi Condensed', regular = '~/Downloads/Barlow_Semi_Condensed/BarlowSemiCondensed-Regular.ttf')
showtext_auto()
theme_clean <- function(axis_title=T) {
  t <- theme_minimal(base_family = "Barlow Semi Condensed", base_size=16) +
    theme(panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", hjust = 0),
          strip.background = element_rect(fill = "grey80", color = NA),
          legend.title = element_blank(),#element_text(face = "bold"),
          legend.position = 'top')
  if (!axis_title) t <- t + theme(axis.title.x = element_blank())
  t
}

# Make labels use Barlow by default
update_geom_defaults("label_repel",
                     list(family = "Barlow Semi Condensed",
                          fontface = "bold"))
update_geom_defaults("label",
                     list(family = "Barlow Semi Condensed",
                          fontface = "bold"))

nested_settings <- strip_nested(
  text_x = list(element_text(family = "Barlow Semi Condensed",
                             face = "plain"), NULL),
  background_x = list(element_rect(fill = "grey93"), NULL),
  by_layer_x = TRUE)
