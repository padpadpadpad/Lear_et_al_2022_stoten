# make map of Falmouth and sampling sites that is used in Figure 1 ####

# load packages
library(tidyverse)
library(ggmap)

# sampling sites latitude and longitude
points <- tibble(lat = c(50.152, 50.165, 50.181),
                 lon = c(-5.059, -5.083, -5.044),
                 label = c('Falmouth\nDocks', 'Falmouth\nWharf', 'Mylor Bank'))

# make bounding box
box <- make_bbox(lon, lat, data = points, f = .10)

# manually change the bounding box limits
box[1:4] <- c(-5.085, 50.145, -5.015, 50.187)
box

# make map
get_stamenmap(box, zoom = 14, maptype = "terrain") %>% 
  ggmap() +
  geom_point(aes(lon, lat), points, size = 6) +
  ggrepel::geom_label_repel(aes(lon, lat, label = label), points, label.padding = 0.4, nudge_x = 0.01, nudge_y = -0.001, size = MicrobioUoE::pts(14)) +
  theme_bw(base_size = 14) +
  labs(x = 'Longitude',
       y = 'Latitude') +
  annotate(y = 50.145, x = -5.015, label = 'Map made using ggmap, map tiles by Stamen Design, under CC BY 3.0', geom = 'label', size = MicrobioUoE::pts(8), hjust = 1, vjust = -0.03, label_padding = unit(0, 'lines'))

ggsave('Figures/Figure_1_map.pdf', last_plot(), height = 6, width = 6)
ggsave('Figures/Figure_1_map.png', last_plot(), height = 6, width = 6)
