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
#### CMF Simulation - vary noise in model ####
nsim <- 50
param_name <- "noise"
params <- c(.1, .5, 1, 1.5, 2)

avg_err_df <- as.data.frame(matrix(NA, nrow = 9, ncol = length(params)))
colnames(avg_err_df) <- params
rownames(avg_err_df) <- c("pca1", "pca2", "pca3", "concatenated", "cmf",
                          "mfa", "jive", "addfrob", "multfrob")

for (param in params) {
  # make sims data
  sims <- cmf_model(nsim = nsim, noise = param)
  
  metric_df <- data.frame(pca1 = NULL, pca2 = NULL, pca3 = NULL, 
                          concatenated = NULL, cmf = NULL, mfa = NULL,
                          jive = NULL, addfrob = NULL, multfrob = NULL)
  for (i in 1:nsim) {
    sim <- sims[[i]]
    sim_data <- sim$sim_data
    truth <- sim$truth
    U_true <- truth$U_true
    dim_U <- 2
    
    Sigh_ls <- list()
    
    # individual PCAs
    for (k in 1:length(sim_data)) {
      pca_name <- paste0("pca", k)
      Sigh_ls[[pca_name]] <- IndividualPCA(dat = sim_data, k = k)$Sig
      metric_df[i, pca_name] <- subspace_recovery(U = U_true, 
                                                  Sigh = Sigh_ls[[pca_name]], 
                                                  dim_U = dim_U)
    }
    
    # concatenated PCA
    Sigh_ls[["concatenated"]] <- ConcatenatedPCA(dat = sim_data)$Sig
    metric_df[i, "concatenated"] <- subspace_recovery(U = U_true,
                                                      Sigh = Sigh_ls[["concatenated"]],
                                                      dim_U = dim_U)
    
    # cmf
    Sigh_ls[["cmf"]] <- CMF(dat = sim_data, r = 2, lamS = 1, lamDs = c(1, 1, 1))$U
    metric_df[i, "cmf"] <- subspace_recovery(U = U_true,
                                             Uh = Sigh_ls[["cmf"]],
                                             dim_U = dim_U)
    
    # MFA
    Sigh_ls[["mfa"]] <- my_MFA(dat = sim_data)$U
    metric_df[i, "mfa"] <- subspace_recovery(U = U_true,
                                             Uh = Sigh_ls[["mfa"]],
                                             dim_U = dim_U)
    
    # JIVE
    Sigh_ls[["jive"]] <- my_JIVE(dat = sim_data)$Sig
    metric_df[i, "jive"] <- subspace_recovery(U = U_true,
                                              Sigh = Sigh_ls[["jive"]],
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
    metric_df[i, "addfrob"] <- subspace_recovery(U = U_true,
                                                 Sigh = Sigh_ls[["addfrob"]],
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
    metric_df[i, "multfrob"] <- subspace_recovery(U = U_true,
                                                  Sigh = Sigh_ls[["multfrob"]],
                                                  dim_U = dim_U)
  }
  
  avg_err_df[as.character(param)] <- colMeans(metric_df)
}

avg_err_df <- cbind(method = rownames(avg_err_df), avg_err_df)
plt_df <- melt(avg_err_df, id.vars = "method")
plt_df$variable <- as.numeric(as.character(plt_df$variable))
ggplot(plt_df) +
  aes(x = variable, y = value, color = method) +
  geom_line()
