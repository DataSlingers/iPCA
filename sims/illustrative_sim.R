library(reshape2)
library(RColorBrewer)
library(R.matlab)
library(randomForest)
library(R.utils)
library(plotly)

library(foreach)
library(doParallel)

sourceDirectory("../functions/", modifiedOnly = F, recursive = F) # useful functions

#### ####
#### Illustrative Example ####
sim <- ipca_base_sim(nsim = 1, n = 200)[[1]]
sim_data <- sim$sim_data
truth <- sim$truth

# initialize list of Sighs
Sigh_ls <- list()

# inidividual pca
for (k in 1:length(sim_data)) {
  pca_name <- paste0("pca", k)
  Sigh_ls[[pca_name]] <- IndividualPCA(dat = sim_data, k = k)$Sig
}

# multfrob
Sigh_ls[["multfrob"]] <- FFmleMultFrob(dat = sim_data, lamDs = c(100, 1, 10000))$Sig

# initialize plot list
plt_ls <- list()
plt_ls[[1]] <- plot_ipca(Sig = truth$Sig_true, y = truth$y_true, show_plot = F,
                         plot_title = "Truth")$ipca_plot # plot truth
for (i in 1:length(Sigh_ls)) {
  plt_ls[[i+1]] <- plot_ipca(Sig = Sigh_ls[[i]], y = truth$y_true, show_plot = F,
                             plot_title = names(Sigh_ls)[i])$ipca_plot
}

# illustrative simulation plot
p_illustrative <- ggarrange(plotlist = plt_ls, nrow = 1, ncol = 5, labels = "AUTO")
p_illustrative
