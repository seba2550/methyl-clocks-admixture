# This script puts together most of the figures in the manuscript.
# The figures have been generated as individual plots in other scripts, saved as RDS files, and are then put into figure grids here

# setwd("~/Desktop/Capra Lab/Thesis Project/Aim_1/plots_RDS/")

library(ggpubr)
library(cowplot)
library(magick)
library(grid)

fig5b <- image_read_pdf("data/fig5b.pdf", density = 300)

# Turn into a grob (grid object)
fig5b_grob <- rasterGrob(as.raster(fig5b), interpolate = TRUE)

# Wrap grob in a cowplot ggdraw so it behaves like a ggplot object
panelB <- ggdraw() +
  draw_grob(fig5b_grob, 
            x = 0, y = 0,     # anchor at bottom-left
            width = 1, height = 1, 
            scale = 1.2)      # >1 enlarges, <1 shrinks

# grid.newpage()
# # Then, draw the grob
# grid.draw(fig5b_grob)

load_plots <- function(path = "data/") {
  files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
  plots <- lapply(files, readRDS)
  names(plots) <- gsub("\\.rds$", "", basename(files))
  plots
}

plots <- load_plots()

fig5a <- plots$fig5a
fig5a <- fig5a + ylab("% Clock CpGs with Differential Methylation in AFR")

fig5c <- plots$fig5c
fig5c <- fig5c + xlab("")

fig5d <- plots$fig5d

fig5e <- plots$fig5e
fig5e <- fig5e & theme(
  axis.title.x = element_blank()
)

fig5f <- plots$fig5f
fig5f <- fig5f + xlab("")

fig5g <- plots$fig5g

# Put together Figure 5
# Make a list in order
figure5_panels <- list(
  plots$fig5a,  # A
  panelB,       # B
  plots$fig5c,  # C
  plots$fig5d,  # D
  plots$fig5e,  # E
  plots$fig5f,  # F
  plots$fig5g   # G
)

# Combine into 3 rows
# -----------------------------
# We'll use 3 rows and auto-calculate columns (max 3 per row)
# ncol = ceiling(length(panels)/nrow) → roughly 3x3 layout
Figure5 <- ggarrange(
  plotlist = figure5_panels,
  labels = LETTERS[1:length(figure5_panels)],
  ncol = 3,
  nrow = 3
)


# Save PDFs
ggsave("manuscript/main_figs/Figure2.pdf", plots$fig2, width = 20, height = 11)
ggsave("manuscript/main_figs/Figure3.pdf", plots$fig3, width = 20, height = 11)
ggsave("manuscript/main_figs/Figure6.pdf", Figure5, width = 20, height = 15)
