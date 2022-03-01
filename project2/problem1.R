### Problem 1 ----
library(spatial)
library(MASS)
library(tidyverse)

## a) ----
# Plot cells point data.
cells <- ppinit("cells.dat")
cells.df <- data.frame(x = cells$x, y = cells$y)
cells$area
ggplot(cells.df) + geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), alpha = 0, color = 1) + theme_minimal()

# Plot redwood point data.
redwood <- ppinit("redwood.dat")
redwood.df <- data.frame(x = redwood$x, y = redwood$y)
redwood$area
ggplot(redwood.df) + geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = 0, xmax = 1, ymin = -1, ymax = 0), alpha = 0, color = 1) + theme_minimal()

# Plot pines point data.
pines <- ppinit("pines.dat")
pines.df <- data.frame(x = pines$x, y = pines$y)
ggplot(pines.df) + geom_point(aes(x = x, y = y)) + 
  geom_rect(aes(xmin = pines$area[1], xmax = pines$area[2], ymin = pines$area[3], ymax = pines$area[4]),
            alpha = 0, color = 1) + theme_minimal()
ggsave("figures/pines.pdf")
