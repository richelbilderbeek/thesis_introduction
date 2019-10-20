# Creates the posteriors for humans

library(pirouette)
library(ggplot2)

root_folder <- path.expand("~/GitHubs/thesis_introduction")
setwd(root_folder)
set.seed(314)

# Read primates
nexus_file <- beastier::get_beast2_example_filename("Primates.nex")
alignment <- ape::as.DNAbin(ape::read.nexus.data(nexus_file))

# Create a FASTA filename
fasta_filename <- tempfile()
ape::write.FASTA(x = alignment, file = fasta_filename)
readLines(fasta_filename)

# Inference models
mcmc <- create_mcmc(chain_length = 1e6, store_every = 1e3)
n_states <- 1 + (mcmc$chain_length / mcmc$store_every)
halfway_index <- (mcmc$chain_length / mcmc$store_every) / 2


# Run babette with Yule model
inference_model_yule <- create_inference_model(
  tree_prior = create_yule_tree_prior(),
  mcmc = mcmc
)
out_yule <- bbt_run_from_model(
  fasta_filename = fasta_filename,
  inference_model = inference_model_yule
)

png("primates_yule.png")
babette::plot_densitree(
  out_yule[[1]][halfway_index:n_states],
  scaleX = TRUE,
  main = "Yule",
  type = "phylogram",
  width = 2,
  alpha = 10 / n_states
)
dev.off()

# Birth-Death
inference_model_bd <- create_inference_model(
  tree_prior = create_bd_tree_prior(),
  mcmc = mcmc
)
out_bd <- bbt_run_from_model(
  fasta_filename = fasta_filename,
  inference_model = inference_model_bd
)
png("primates_bd.png")
babette::plot_densitree(
  out_bd[[1]][halfway_index:n_states],
  scaleX = TRUE,
  main = "BD",
  type = "phylogram",
  width = 2,
  alpha = 10 / n_states
)
dev.off()

# Run babette with CEP
inference_model_cep <- create_inference_model(
  tree_prior = create_cep_tree_prior(),
  mcmc = mcmc
)
out_cep <- bbt_run_from_model(
  fasta_filename = fasta_filename,
  inference_model = inference_model_cep
)
png("primates_cep.png")
babette::plot_densitree(
  out_cep[[1]][halfway_index:n_states],
  scaleX = TRUE,
  main = "CEP",
  type = "phylogram",
  width = 2,
  alpha = 10 / n_states
)
dev.off()


if (1 == 2) {

  install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")
  require(PhyloOrchard)
  data(BinindaEmondsEtAl2007)
  BinindaEmondsEtAl2007 <- .compressTipLabel(BinindaEmondsEtAl2007)
  phangorn::densiTree(BinindaEmondsEtAl2007, type = "phylogram", col="red", main = "Hello")



  png("~/title_around.png")
  plot(1:10, main = "my_title")
  dev.off()

  png("~/title_missing.png")
  phangorn::densiTree(c(ape::rcoal(n = 10)), main = "my_title")
  dev.off()


  phangorn::densiTree(c(ape::rcoal(n = 10)), main = "Hello")
  ?plot
  ?phangorn::densiTree
}
