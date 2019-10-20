# Effect of species tree prior on inference

library(pirouette)
library(ggplot2)

set.seed(314)
root_folder <- path.expand("~/GitHubs/thesis_introduction")
setwd(root_folder)

true_trees <- TreeSim::sim.bd.taxa(n = 9, numbsim = 1, lambda = 10.0, mu = 9.0, complete = FALSE)
true_tree <- true_trees[[1]]
true_tree$tip.label <- paste0("t", seq(1, 10))
ape::plot.phylo(true_tree)

# Inference models
mcmc <- create_mcmc(chain_length = 1e6, store_every = 1e3)
n_states <- 1 + (mcmc$chain_length / mcmc$store_every)
halfway_index <- (mcmc$chain_length / mcmc$store_every) / 2

alignment_params <- create_test_alignment_params(
  root_sequence = create_blocked_dna(1000)
)

mrca_prior <- create_mrca_prior(
  mrca_distr = create_normal_distr(mean = 1.0, sigma = 0.0001),
  is_monophyletic = TRUE
)

pir_params_yule <- create_test_pir_params(
  alignment_params = alignment_params,
  experiments = list(
      create_test_gen_experiment(
      inference_model = create_inference_model(
        tree_prior = create_yule_tree_prior(),
        mcmc = mcmc,
        mrca_prior = mrca_prior
      )
    )
  )
)
pir_params_bd <- create_test_pir_params(
  alignment_params = alignment_params,
  experiments = list(
      create_test_gen_experiment(
      inference_model = create_inference_model(
        tree_prior = create_bd_tree_prior(),
        mcmc = mcmc,
        mrca_prior = mrca_prior
      )
    )
  )
)


out_yule <- pirouette::pir_run(
  phylogeny = true_tree,
  pir_params = pir_params_yule
)

out_bd <- pirouette::pir_run(
  phylogeny = true_tree,
  pir_params = pir_params_bd
)

png("bd_tree_yule_prior.png")
babette::plot_densitree(
  tracerer::parse_beast_trees(
    pir_params_yule$experiments[[1]]$beast2_options$output_trees_filenames
  )[halfway_index:n_states],
  alpha = 10 / n_states,
  type = "phylogram",
  consensus = sort(true_tree$tip.label)
)
dev.off()

png("bd_tree_bd_prior.png")
babette::plot_densitree(
  tracerer::parse_beast_trees(
    pir_params_bd$experiments[[1]]$beast2_options$output_trees_filenames
  )[halfway_index:n_states],
  alpha = 10 / n_states,
  type = "phylogram",
  consensus = sort(true_tree$tip.label)
)
dev.off()
