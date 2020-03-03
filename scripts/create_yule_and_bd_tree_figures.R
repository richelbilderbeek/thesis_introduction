# Creates the figures for the Yule and BD trees

library(pirouette)
library(ggplot2)
library(ggtree)

set.seed(314)
tree <- pirouette::create_yule_tree(n_taxa = 10)

grDevices::png(filename = "~/GitHubs/thesis_introduction/yule_tree.png", width = 800, height = 600)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 5)
ggtree::ggtree(tree, size = 2, col = "black")  + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

grDevices::png(filename = "~/GitHubs/thesis_introduction/yule_tree_nltt.png", width = 800, height = 600)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 5)
nLTT::nltts_plot(c(tree), col = "black", lwd = 7, main = "nLTT", cex.lab = 2.5, cex.main = 4, cex.axis = 2.5)
grDevices::dev.off()


set.seed(315)
tree <- pirouette::create_bd_tree(n_taxa = 10)

grDevices::png(filename = "~/GitHubs/thesis_introduction/bd_tree.png", width = 800, height = 600)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 5)
ggtree::ggtree(tree, size = 2, col = "black")  + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

grDevices::png(filename = "~/GitHubs/thesis_introduction/bd_tree_nltt.png", width = 800, height = 600)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 5)
nLTT::nltts_plot(c(tree), col = "black", lwd = 7, main = "nLTT", cex.lab = 2.5, cex.main = 4, cex.axis = 2.5)
grDevices::dev.off()

