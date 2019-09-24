# Creates the nLTT figures

library(pirouette)
library(ggplot2)

root_folder <- path.expand("~/GitHubs/thesis/introduction")
set.seed(314)

true_tree <- ape::read.tree(text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);")
twin_tree <- pirouette::create_twin_tree(
  phylogeny = true_tree,
  twinning_params = create_twinning_params(rng_seed = 314)
)

################################################################################
# Plotting
################################################################################
grDevices::png(filename = file.path(root_folder, "true_tree.png"), width = 800, height = 600)
# ape::plot.phylo(phylogeny, cex = 2.0, edge.width = 2.0)
ggtree::ggtree(true_tree, size = 2, col = "red") + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

grDevices::png(filename = file.path(root_folder, "twin_tree.png"), width = 800, height = 600)
# ape::plot.phylo(phylogeny, cex = 2.0, edge.width = 2.0)
ggtree::ggtree(twin_tree, size = 2, col = "blue") + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

grDevices::png(filename = file.path(root_folder, "nltt.png"), width = 800, height = 600)
par(mar = c(5.1, 4.1, 4.1, 2.1) + 5)
nLTT::nltts_plot(c(true_tree), col = "red", lwd = 7, main = "nLTTs", cex.lab = 2.5, cex.main = 4, cex.axis = 2.5)
nLTT::nltts_plot(c(twin_tree), replot = TRUE, col = "blue", lwd = 7)
grDevices::dev.off()

