# Creates the figures for the workflow image in the documentation

library(pirouette)
library(ggplot2)

dna_sequence_length <- 40
root_folder <- path.expand(
  paste0(
    "~/GitHubs/thesis_introduction/scripts/",
    dna_sequence_length
  )
)
set.seed(314)

phylogeny  <- ape::read.tree(text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);")
n_taxa <- ape::Ntip(phylogeny)
crown_age <- max(ape::branching.times(phylogeny))
################################################################################
# Use png here, to make Peregrine fail fast
################################################################################
grDevices::png(filename = file.path(root_folder, "phylogeny.png"), width = 800, height = 600)
# ape::plot.phylo(phylogeny, cex = 2.0, edge.width = 2.0)
ggtree::ggtree(phylogeny, size = 2) + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()


pir_params <- create_pir_params(
  alignment_params = create_alignment_params(
    root_sequence = pirouette::create_blocked_dna(length = dna_sequence_length),
    mutation_rate = 0.5 / crown_age
  ),
  twinning_params = create_twinning_params()
)

################################################################################
# Settings to run on Peregrine cluster
################################################################################
if (1 == 2) {
  pir_params$alignment_params$fasta_filename <- file.path(root_folder, "alignment.fasta")
  pir_params$experiments[[1]]$beast2_options$input_filename <- file.path(root_folder, "beast2_input.xml")
  pir_params$experiments[[1]]$beast2_options$output_log_filename <- file.path(root_folder, "beast2_output.log")
  pir_params$experiments[[1]]$beast2_options$output_trees_filenames <- file.path(root_folder, "beast2_output.trees")
  pir_params$experiments[[1]]$beast2_options$output_state_filename <- file.path(root_folder, "beast2_output.xml.state")
  pir_params$experiments[[1]]$errors_filename <- file.path(root_folder, "error.csv")
}
################################################################################

errors <- pir_run(phylogeny = phylogeny, pir_params = pir_params)

################################################################################
# trees
################################################################################

grDevices::png(filename = file.path(root_folder, "phylogeny.png"), width = 800, height = 600)
# ape::plot.phylo(phylogeny, cex = 2.0, edge.width = 2.0)
ggtree::ggtree(phylogeny, size = 2) + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()


grDevices::png(filename = file.path(root_folder, "phylogeny_twin.png"), width = 800, height = 600)
# ape::plot.phylo(ape::read.tree(pir_params$twinning_params$twin_tree_filename), cex = 2.0, edge.width = 2.0)
ggtree::ggtree(ape::read.tree(pir_params$twinning_params$twin_tree_filename), size = 2) + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()


################################################################################
# alignment
################################################################################

n_mutations <- count_n_mutations(
  ape::read.FASTA(pir_params$alignment_params$fasta_filename),
  root_sequence = pir_params$alignment_params$root_sequence
)
print(paste("n_mutations", n_mutations))
print(paste("n_mutations per taxon", n_mutations / n_taxa))

grDevices::png(filename = file.path(root_folder, "alignment.png"), width = 800, height = 300)
ape::image.DNAbin(
  ape::read.FASTA(file = pir_params$alignment_params$fasta_filename),
  grid = TRUE,
  show.bases = FALSE,
  legend = FALSE,
  cex.lab = 2.0,
  cex.axis = 2.0
)
grDevices::dev.off()

# Add annotation
fasta <- readLines(pir_params$alignment_params$fasta_filename)
fasta <- c(">", pir_params$alignment_params$root_sequence, fasta)
fasta_filename <- file.path(root_folder, "alignment_with_root.fas")
writeLines(text = fasta, con = fasta_filename)
png_filename <- file.path(root_folder, "alignment_with_root.png")
grDevices::png(filename = png_filename, width = 800, height = 300)
ape::image.DNAbin(
  ape::read.FASTA(file = fasta_filename),
  grid = TRUE,
  show.bases = FALSE,
  legend = FALSE,
  cex.lab = 2.0,
  cex.axis = 2.0
)
grDevices::dev.off()

grDevices::png(filename = file.path(root_folder, "alignment_twin.png"), width = 800, height = 300)
ape::image.DNAbin(
  ape::read.FASTA(file = pir_params$twinning_params$twin_alignment_filename),
  grid = TRUE,
  show.bases = FALSE,
  legend = FALSE,
  cex.lab = 2.0,
  cex.axis = 2.0
)
grDevices::dev.off()

################################################################################
# posteriors
################################################################################

grDevices::png(filename = file.path(root_folder, "densitree.png"), width = 1000, height = 800)
babette::plot_densitree(
  phylos = tracerer::parse_beast_trees(pir_params$experiments[[1]]$beast2_options$output_trees_filenames),
  alpha = 0.01,
  consensus = LETTERS[6:1],
  cex = 6.0,
  scaleX = TRUE,
  scale.bar = FALSE
)
grDevices::dev.off()

grDevices::png(filename = file.path(root_folder, "densitree_twin.png"), width = 1000, height = 800)
babette::plot_densitree(
  phylos = tracerer::parse_beast_trees(to_twin_filename(pir_params$experiments[[1]]$beast2_options$output_trees_filenames)),
  alpha = 0.01,
  consensus = LETTERS[6:1],
  cex = 6.0,
  scaleX = TRUE,
  scale.bar = FALSE
)
grDevices::dev.off()

################################################################################
# histogram of errors
################################################################################
df_errors <- data.frame(error = utils::read.csv(pir_params$experiments[[1]]$errors_filename)$x)
df_errors_twin <- data.frame(error = utils::read.csv(to_twin_filename(pir_params$experiments[[1]]$errors_filename))$x)

ggplot2::ggplot(
  df_errors,
  aes(x = error)
) + geom_histogram(binwidth = 0.01) +
  ggplot2::theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 30, face = "bold")
  ) + ggsave(file.path(root_folder, "errors.png"))

ggplot2::ggplot(
  df_errors_twin,
  aes(x = error)
) + geom_histogram(binwidth = 0.01) +
  ggplot2::theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 30, face = "bold")
  ) + ggsave(file.path(root_folder, "errors_twin.png"))

ggplot2::ggplot(
  df_errors,
  aes(x = "", y = error)
) + geom_violin() +
  xlab("") +
  scale_y_continuous(breaks = seq(0.0, 1.0, by = 0.02)) +
  ggsave(file.path(root_folder, "errors_violin.png"))

ggplot2::ggplot(
  df_errors_twin,
  aes(x = "", y = error)
) + geom_violin() +
  xlab("") +
  scale_y_continuous(breaks = seq(0.0, 1.0, by = 0.02)) +
  ggsave(file.path(root_folder, "errors_violin_twin.png"))


################################################################################
# Create ML tree from alignment
################################################################################
alignment <- ape::read.FASTA(pir_params$alignment_params$fasta_filename)

dm  <- phangorn::dist.ml(alignment)
tree_upgma  <- phangorn::upgma(dm)
tree_nj  <- phangorn::NJ(dm)
tree_fastme  <- ape::fastme.bal(dm)
tree_bionj  <- ape::bionj(dm)

grDevices::png(filename = file.path(root_folder, "phylogeny_upgma.png"), width = 800, height = 600)
ggtree::ggtree(tree_upgma, size = 2) + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

grDevices::png(filename = file.path(root_folder, "phylogeny_nj.png"), width = 800, height = 600)
ggtree::ggtree(tree_nj, size = 2) + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

grDevices::png(filename = file.path(root_folder, "phylogeny_fastme.png"), width = 800, height = 600)
ggtree::ggtree(tree_fastme, size = 2) + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

grDevices::png(filename = file.path(root_folder, "phylogeny_bionj.png"), width = 800, height = 600)
ggtree::ggtree(tree_bionj, size = 2) + ggtree::geom_tiplab(size = 16) + ggplot2::theme(plot.margin = unit(c(1,1,1,1), "cm"))
grDevices::dev.off()

