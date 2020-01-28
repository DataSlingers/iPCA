library(reshape2)
library(RColorBrewer)
library(R.matlab)
library(randomForest)
library(R.utils)
library(plotly)
library(foreach)
library(doParallel)
library(SpatioTemporal)

sourceDirectory("../functions/", modifiedOnly = F, recursive = F) # useful functions

#### ####
#### iPCA Simulation with Different/Uncommon Sigmas ####
nsim <- 50
param_name <- "diff_idx"
params <- 1:3

avg_err_df <- as.data.frame(matrix(NA, nrow = 15, ncol = length(params)))
colnames(avg_err_df) <- params
rownames(avg_err_df) <- c("pca1", "pca2", "pca3", 
                          "concatenated", "concatenated_without",
                          "mfa", "mfa_without", "jive", "jive_without",
                          "addfrob", "addfrob_without",
                          "multfrob", "multfrob_without",
                          "l1", "l1_without")

for (param in params) {
  # make sims data
  sims <- ipca_model_diff_sigmas(nsim = nsim, diff_idx = param)
  
  metric_df <- data.frame(pca1 = NULL, pca2 = NULL, pca3 = NULL, 
                          concatenated = NULL, concatenated_without = NULL, 
                          mfa = NULL, mfa_without = NULL,
                          jive = NULL, jive_without = NULL,
                          addfrob = NULL, addfrob_without = NULL,
                          multfrob = NULL, multfrob_without = NULL, 
                          l1 = NULL, l1_without = NULL)
  
  for (i in 1:nsim) {
    sim <- sims[[i]]
    sim_data <- sim$sim_data
    sim_data_without <- sim$sim_data[-param]
    truth <- sim$truth
    Sig_true <- truth$Sig_true
    dim_U <- 2
    
    # initialize list for Sighs
    Sigh_ls <- list()
    
    # individual PCAs
    for (k in 1:length(sim_data)) {
      pca_name <- paste0("pca", k)
      Sigh_ls[[pca_name]] <- IndividualPCA(dat = sim_data, k = k)$Sig
      metric_df[i, pca_name] <- subspace_recovery(Sig = Sig_true, 
                                                  Sigh = Sigh_ls[[pca_name]], 
                                                  dim_U = dim_U)
    }
    
    # concatenated PCA
    Sigh_ls[["concatenated"]] <- ConcatenatedPCA(dat = sim_data)$Sig
    metric_df[i, "concatenated"] <- subspace_recovery(Sig = Sig_true,
                                                      Sigh = Sigh_ls[["concatenated"]],
                                                      dim_U = dim_U)
    Sigh_ls[["concatenated_without"]] <- ConcatenatedPCA(dat = sim_data_without)$Sig
    metric_df[i, "concatenated_without"] <- subspace_recovery(Sig = Sig_true,
                                                      Sigh = Sigh_ls[["concatenated_without"]],
                                                      dim_U = dim_U)
    
    # MFA
    Sigh_ls[["mfa"]] <- my_MFA(dat = sim_data)$U
    metric_df[i, "mfa"] <- subspace_recovery(Sig = Sig_true,
                                             Uh = Sigh_ls[["mfa"]],
                                             dim_U = dim_U)
    Sigh_ls[["mfa_without"]] <- my_MFA(dat = sim_data_without)$U
    metric_df[i, "mfa_without"] <- subspace_recovery(Sig = Sig_true,
                                             Uh = Sigh_ls[["mfa_without"]],
                                             dim_U = dim_U)
    
    # JIVE
    Sigh_ls[["jive"]] <- my_JIVE(dat = sim_data)$Sig
    metric_df[i, "jive"] <- subspace_recovery(Sig = Sig_true,
                                              Sigh = Sigh_ls[["jive"]],
                                              dim_U = dim_U)
    Sigh_ls[["jive_without"]] <- my_JIVE(dat = sim_data_without)$Sig
    metric_df[i, "jive_without"] <- subspace_recovery(Sig = Sig_true,
                                              Sigh = Sigh_ls[["jive_without"]],
                                              dim_U = dim_U)
    
    # Addfrob
    # to reduce computational time, make lam_grid smaller and choose fewer lams
    lams <- c(1e-4, 1e-2, 1, 100, 1000, 10000, 100000) # for additive penalties
    lam_grid <- expand.grid(lams, lams, lams, lams)
    choose_lambdas_ans <- choose_lambdas(dat = sim_data, lam_grid = lam_grid, 
                                         q = "addfrob", trcma = T, maxit = 10, 
                                         greedy.search = T, maxit.search = 1,
                                         seed = sample(x = 1:10000, 1))
    best_lambdas <- choose_lambdas_ans$best_lambdas
    Sigh_ls[["addfrob"]] <- FFmleAddFrob(dat = sim_data, 
                                         lamDs = best_lambdas[-1],
                                         lamS = best_lambdas[1])$Sig
    metric_df[i, "addfrob"] <- subspace_recovery(Sig = Sig_true,
                                                 Sigh = Sigh_ls[["addfrob"]],
                                                 dim_U = dim_U)
    lam_grid <- expand.grid(lams, lams, lams)
    choose_lambdas_ans <- choose_lambdas(dat = sim_data_without, lam_grid = lam_grid, 
                                         q = "addfrob", trcma = T, maxit = 10, 
                                         greedy.search = T, maxit.search = 1,
                                         seed = sample(x = 1:10000, 1))
    best_lambdas <- choose_lambdas_ans$best_lambdas
    Sigh_ls[["addfrob_without"]] <- FFmleAddFrob(dat = sim_data_without, 
                                         lamDs = best_lambdas[-1],
                                         lamS = best_lambdas[1])$Sig
    metric_df[i, "addfrob_without"] <- subspace_recovery(Sig = Sig_true,
                                                 Sigh = Sigh_ls[["addfrob_without"]],
                                                 dim_U = dim_U)
    
    # Multfrob
    # to reduce computational time, make lam_grid smaller and choose fewer lams
    lams <- c(1e-4, 1e-2, 1, 10, 100, 1000, 10000) # for multiplicative penalties
    lam_grid <- expand.grid(lams, lams, lams)
    choose_lambdas_ans <- choose_lambdas(dat = sim_data, lam_grid = lam_grid, 
                                         q = "multfrob", trcma = T, maxit = 10, 
                                         greedy.search = T, maxit.search = 1,
                                         seed = sample(x = 1:10000, 1))
    best_lambdas <- choose_lambdas_ans$best_lambdas
    Sigh_ls[["multfrob"]] <- FFmleMultFrob(dat = sim_data, 
                                           lamDs = best_lambdas)$Sig
    metric_df[i, "multfrob"] <- subspace_recovery(Sig = Sig_true,
                                                  Sigh = Sigh_ls[["multfrob"]],
                                                  dim_U = dim_U)
    lam_grid <- expand.grid(lams, lams)
    choose_lambdas_ans <- choose_lambdas(dat = sim_data_without, lam_grid = lam_grid, 
                                         q = "multfrob", trcma = T, maxit = 10, 
                                         greedy.search = T, maxit.search = 1,
                                         seed = sample(x = 1:10000, 1))
    best_lambdas <- choose_lambdas_ans$best_lambdas
    Sigh_ls[["multfrob_without"]] <- FFmleMultFrob(dat = sim_data_without, 
                                           lamDs = best_lambdas)$Sig
    metric_df[i, "multfrob_without"] <- subspace_recovery(Sig = Sig_true,
                                                  Sigh = Sigh_ls[["multfrob_without"]],
                                                  dim_U = dim_U)
    
    # L1
    # to reduce computational time, make lam_grid smaller and choose fewer lams
    lams <- c(1e-4, 1e-2, 1, 100, 1000, 10000, 100000) # for additive penalties
    lam_grid <- expand.grid(lams, lams, lams, lams)
    choose_lambdas_ans <- choose_lambdas(dat = sim_data, lam_grid = lam_grid, 
                                         q = "1_off", trcma = T, maxit = 10, 
                                         greedy.search = T, maxit.search = 1,
                                         seed = sample(x = 1:10000, 1))
    best_lambdas <- choose_lambdas_ans$best_lambdas
    Sigh_ls[["l1"]] <- FFmleGlasso(dat = sim_data, 
                                   lamDs = best_lambdas[-1],
                                   lamS = best_lambdas[1], 
                                   maxit = 1, pen_diag = F)$Sig
    metric_df[i, "l1"] <- subspace_recovery(Sig = Sig_true,
                                            Sigh = Sigh_ls[["l1"]],
                                            dim_U = dim_U)
    lam_grid <- expand.grid(lams, lams, lams)
    choose_lambdas_ans <- choose_lambdas(dat = sim_data_without, lam_grid = lam_grid, 
                                         q = "1_off", trcma = T, maxit = 10, 
                                         greedy.search = T, maxit.search = 1,
                                         seed = sample(x = 1:10000, 1))
    best_lambdas <- choose_lambdas_ans$best_lambdas
    Sigh_ls[["l1_without"]] <- FFmleGlasso(dat = sim_data_without, 
                                   lamDs = best_lambdas[-1],
                                   lamS = best_lambdas[1], 
                                   maxit = 1, pen_diag = F)$Sig
    metric_df[i, "l1_without"] <- subspace_recovery(Sig = Sig_true,
                                            Sigh = Sigh_ls[["l1_without"]],
                                            dim_U = dim_U)
  }
  
  avg_err_df[as.character(param)] <- colMeans(metric_df)
  
}

avg_err_df <- cbind(method = rownames(avg_err_df), avg_err_df)
plt_df <- melt(avg_err_df, id.vars = "method")
plt_df$variable <- as.numeric(as.character(plt_df$variable))
ggplot(plt_df) +
  aes(x = method, y = value, color = method, fill = method) +
  facet_grid(variable ~ ., scales = "free") +
  geom_bar(stat = "identity", position = "dodge", alpha = .75, width = .9)