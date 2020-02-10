# Creates the figures for the introduction
setwd("~/GitHubs/thesis_introduction")

suppressMessages(library(ggtree))


#
#      +--------x   A
#   ---+
#      |     +------- B
#      +-----+
#            +------- C
#
#  -+--+-----+--+---- (time)
#   t0 t1    t2 t3
#

ggtree::ggtree(
  ape::read.tree(text = "((C:2,B:2):1, A:2):1;"),
  ladderize = FALSE
) + ggtree::geom_tiplab(cex = 10) +
  ggtree::geom_rootedge() +
  ggtree::geom_treescale(width = 6, x = -2, fontsize = 0) +
  ggplot2::geom_vline(xintercept = -1, lty = "dotted") +
  ggplot2::geom_vline(xintercept = 0, lty = "dotted") +
  ggplot2::geom_vline(xintercept = 1, lty = "dotted") +
  ggplot2::geom_vline(xintercept = 2, lty = "dotted") +
  ggplot2::geom_vline(xintercept = 3, lty = "dotted") +
  ggplot2::geom_label(
    data = data.frame(x = seq(-1, 3), y = 0, label = paste0("t", seq(0,4))),
    aes(x = x, y = y, label = label),
    label.size = 0.0,
    nudge_x = 0,
    nudge_y = -0.18,
    cex = 10
  ) +
  ggplot2::geom_label(
    data = data.frame(x = 4, y = 0, label = "t"),
    aes(x = x, y = y, label = label),
    label.size = 0.0,
    nudge_x = 0.1,
    nudge_y = 0.0,
    cex = 10
  ) + ggplot2::ggsave(filename = "t_4.png", width = 7, height = 4)
