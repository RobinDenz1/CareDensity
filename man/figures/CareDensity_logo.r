
library(ggplot2)
library(hexSticker)
library(sysfonts)

# 4 patients
patients <- data.frame(
  x = 0,
  y = c(3.6, 2.4, 1.2, 0)
)

# 3 providers
providers <- data.frame(
  x = 3,
  y = c(3.0, 1.8, 0.6)
)

# edges between them
edges <- data.frame(
  x = c(0, 0, 0, 0, 0, 0),
  y = c(3.6, 3.6, 2.4, 1.2, 1.2, 0),
  xend = c(3, 3, 3, 3, 3, 3),
  yend = c(3.0, 1.8, 3.0, 1.8, 0.6, 0.6)
)

# define the colors
color_patient <- "#0F766E"
color_provider <- "#F59E0B"
color_connect <- "#94A3B8"
color_border <- color_text <- "#134E4A"
color_background <- "#E7F5F0"

# the main plot
p <- ggplot() +
  geom_segment(data=edges, aes(x=x, y=y, xend=xend, yend=yend),
               linewidth=1, color=color_connect, alpha=0.8, lineend="round") +
  geom_point(data=patients, aes(x=x, y=y), size=5, shape=21, fill=color_patient,
             stroke=1) +
  geom_point(data=providers, aes(x=x, y=y), size=5, shape=22, fill=color_provider,
             stroke=1) +
  coord_cartesian(xlim=c(-0.8, 3.8), ylim=c(-0.7, 4.3), expand=FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill="transparent", color=NA))

# create the sticker
s <- hexSticker::sticker(
  subplot = p,
  s_x = 1,
  s_y = 0.82,
  s_width = 1,
  s_height = 1,
  package = "CareDensity",
  p_x = 1,
  p_y = 1.5,
  p_color = color_text,
  p_fontface = "bold",
  p_size = 16,
  p_family="sans",
  h_size = 1.5,
  h_fill = color_background,
  h_color = color_border,
  spotlight = FALSE,
  filename = "logo.png"
)
plot(s)

# alternative color palettes
color_patient <- "#2563EB"
color_provider <- "#14B8A6"
color_connect <- "#94A3B8"
color_border <- color_text <- "#172554"
color_background <- "#EAF8F7" #"#F8FAFC"

color_patient <- "#247BA0"
color_provider <- "#F4A261"
color_connect <- "#8DB3C7"
color_background <- "#F7FAFC"
color_border <- "#173F5F"
color_text <- "#173F5F"

color_patient <- "#4F46E5"
color_provider <- "#F97316"
color_connect <- "#A5B4FC"
color_border <- color_text <- "#312E81"
color_background <- "#FAFAF9"

color_patient <- "#1D4ED8"
color_provider <- "#10B981"
color_connect <- "#93C5B8"
color_border <- color_text <- "#172554"
color_background <- "#F8FAFC"
