library(reshape2)
library(RColorBrewer)
library(R.matlab)
library(simulator)
library(randomForest)
library(R.utils)
library(plotly)

library(foreach)
library(doParallel)

sourceDirectory("../funs/", modifiedOnly = F, recursive = F) # useful functions

#### ####
#### EM CV and Convergence Plots ####
mypath <- "./em_sim/"

# simulate small datasets
sim <- ipca_base_sim(n = 50, pks = c(60,70,70), snr = 1, nsim = 1) 

# only take two datasets for simplicity
sim_data <- sim[[1]]$sim_data[c(1,3)]
truth <- sim[[1]]$truth
truth$Deltks_true <- truth$Deltks_true[c(1,3)]

# possible choices for lambdas
lam_grid <- expand.grid(c(10^(-3), 10^(-2.5), 10^(-2), 10^(-1.5), 10^(-1), 10^(-.5), 10^0, 10^.5, 10^1, 10^1.5, 10^2),
                        c(10^0, 10^.5, 10^1, 10^1.5, 10^2, 10^2.5, 10^3))

trcma_times <- list()
for (trial in 1:10) { # trcma
  rand_seed <- sample(1:1000,1)
  
  ptm <- proc.time()
  ans <- choose_lambdas(dat = sim_dat, q = "multfrob", lam_grid = lam_grid,
                        prop_miss = 0.05, trcma = T, greedy.search = F, 
                        save.filepath = mypath, save.fileext = paste0("_trcma_",trial),
                        seed = rand_seed)
  trcma_times[[trial]] <- proc.time() - ptm
}

mcecm_times <- list()
for (trial in 1:10) { # full mcecm
  rand_seed <- sample(1:1000,1)

  ptm <- proc.time()
  ans <- choose_lambdas(dat = sim_dat, q = "multfrob", lam_grid = lam_grid,
                        prop_miss = 0.05, trcma = F, greedy.search = F,
                        save.filepath = mypath, save.fileext = paste0("_mcecm_maxit_",trial),
                        seed = rand_seed, maxit = 5)
  mcecm_times[[trial]] <- proc.time() - ptm
}

trcma_err_tab <- mcecm_err_tab <- matrix(NA, nrow = nrow(lam_grid), ncol = 10)
for (trial in 1:10) {
  trcma_ans <- readRDS(paste0(mypath, "lambda_errors_trcma_", trial, ".rds"))
  mcecm_ans <- readRDS(paste0(mypath, "lambda_errors_mcecm_maxit_", trial, ".rds"))
  trcma_err_tab[,trial] <- trcma_ans
  mcecm_err_tab[,trial] <- mcecm_ans
}
trcma_err <- rowMeans(trcma_err_tab)
mcecm_err <- rowMeans(mcecm_err_tab)
trcma_se <- apply(trcma_err_tab, 1, sd)/sqrt(10)
mcecm_se <- apply(mcecm_err_tab, 1, sd)/sqrt(10)

(trcma_best_lams <- as.numeric(lam_grid[which.min(trcma_err),]))
(mcecm_best_lams <- as.numeric(lam_grid[which.min(mcecm_err),]))

# cv plots
p_trcma <- p_mcecm <- list()
plt_idx <- 1
for (i in 2:1) {
  trcma_err_plot <- trcma_err[which(lam_grid[,i] == trcma_best_lams[i])]
  mcecm_err_plot <- mcecm_err[which(lam_grid[,i] == mcecm_best_lams[i])]
  trcma_se_plot <- trcma_se[which(lam_grid[,i] == trcma_best_lams[i])]
  mcecm_se_plot <- mcecm_se[which(lam_grid[,i] == mcecm_best_lams[i])]
  
  trcma_vary_lams <- lam_grid[lam_grid[,i] == trcma_best_lams[i], -i]
  mcecm_vary_lams <- lam_grid[lam_grid[,i] == mcecm_best_lams[i], -i]
  
  trcma_plot_df <- data.frame(Lambda = trcma_vary_lams,
                              Error = trcma_err_plot,
                              SE = trcma_se_plot)
  mcecm_plot_df <- data.frame(Lambda = mcecm_vary_lams,
                              Error = mcecm_err_plot,
                              SE = mcecm_se_plot)
  
  my_trcma_title <- ifelse(plt_idx==1, "One-Step Approximation", "")
  my_mcecm_title <- ifelse(i==2, "Full MCECM", "")
  my_xlab <- ifelse(i==2, expression(lambda[1]), expression(lambda[2]))
  my_width <- ifelse(i==2, .25, .15)
  if (i == 2) { my_breaks <- 10^(-3:2) } else { my_breaks <- 10^(0:3)}
  
  p_trcma[[plt_idx]] <- ggplot(trcma_plot_df) +
    aes(x = Lambda, y = Error) +
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = Error - SE, ymax = Error + SE), width = my_width, size = .3) +
    labs(x = my_xlab, y = "Imputation Error", title = my_trcma_title) +
    scale_x_continuous(trans = "log10", breaks = my_breaks)
  
  p_mcecm[[plt_idx]] <- ggplot(mcecm_plot_df) +
    aes(x = Lambda, y = Error) +
    geom_point() + geom_line() +
    geom_errorbar(aes(ymin = Error - SE, ymax = Error + SE), width = my_width, size = .3) +
    labs(x = my_xlab, y = "Imputation Error", title = my_mcecm_title) +
    scale_x_continuous(trans = "log10", breaks = my_breaks)
  plt_idx <- plt_idx + 1
}

p_cv <- c(p_trcma, p_mcecm)
p <- ggarrange(plotlist=p_cv, ncol=2, nrow=2, labels=c("A", "", "B", ""))
p
saveRDS(p, paste0(mypath, "em_cv_plots.rds"))
ggsave(paste0(mypath, "em_cv_plots.pdf"), p, width = 8, height = 8)

# find best lambdas according to subspace recovery metric -> for datasets 1&3: c(1e-2.5, 1e1)
dim_U <- 2
Sig <- truth$Sig_true
Sig_eigs <- eigs(Sig, dim_U, which = "LM"); U <- Sig_eigs$vectors
L <- nrow(lam_grid)
sr_metric <- rep(NA, L)
for (l in 1:L) {
  est <- FFmleMultFrob(dat = sim_dat, lamDs = as.numeric(lam_grid[l,]))
  Sigh <- est$Sig
  
  Sigh_eigs <- eigs(Sigh, dim_U, which = "LM")
  Uh <- Sigh_eigs$vectors
  
  # subspace recovery metric
  sr_metric[l] <- norm(U %*% solve(t(U) %*% U) %*% t(U) - 
                         Uh %*% solve(t(Uh) %*% Uh) %*% t(Uh), "F")^2
  saveRDS(sr_metric, paste0(mypath, "sr_metric.rds"))
}

# subspace recovery as vary lambdas individually
sr_metric <- readRDS(paste0(mypath, "sr_metric.rds"))
(sr_best_lams <- as.numeric(lam_grid[which.min(sr_metric),]))
p_sr <- list(); plt_idx <- 1
for (i in 2:1) {
  sr_plot <- sr_metric[which(lam_grid[,i] == sr_best_lams[i])]
  sr_vary_lams <- lam_grid[lam_grid[,i] == sr_best_lams[i], -i]
  sr_plot_df <- data.frame(Lambda = sr_vary_lams, Error = sr_plot)
  my_xlab <- ifelse(i==2, expression(lambda[1]), expression(lambda[2]))
  if (i == 2) { my_breaks <- 10^(-3:1) } else { my_breaks <- 10^(0:3)}
  my_title <- ifelse(plt_idx == 1, "Subspace Recovery Error Plots", "")
  p_sr[[plt_idx]] <- ggplot(sr_plot_df) +
    aes(x = Lambda, y = Error) + 
    geom_point() + geom_line(size = .75) +
    labs(x = my_xlab, y = "Subspace Recovery Error", title = my_title) + 
    scale_x_continuous(trans = "log10", breaks = my_breaks)
  plt_idx <- plt_idx + 1
}
p_true <- ggarrange(plotlist=p_sr, ncol=2, nrow=1, labels=c("A", ""))
p_true
saveRDS(p_true, paste0(mypath, "em_cv_plots_sr.rds"))
ggsave(paste0(mypath, "em_cv_plots_sr.pdf"), p_true, width = 8, height = 4)


## make convergence plots
Xs <- sim_dat
n <- nrow(Xs[[1]]); pks <- sapply(Xs, FUN = ncol); p <- sum(pks); K <- length(Xs); prop_miss <- 0.05
Xmiss <- Xs
for (k in 1:K) {
  Xk <- Xs[[k]]
  miss_idx <- sample(1:length(Xk), prop_miss*length(Xk), replace = F)
  Xmiss[[k]][miss_idx] <- NA
}

lamDs <- c(1,1)
trcma_ll <- TRCMAimpute_iPCA(x=Xmiss, lamDs=lamDs, q="multfrob", seed=1, trace=T)
mcecm_ll <- RMNimpute_iPCA(x=Xmiss, lamDs=lamDs, q="multfrob", seed=1, trace=T, maxit=20, thr.em = 1e-12)

trcma_ll_df <- data.frame(Algorithm = "One Step Approx.", 
                          Iteration = 1:length(trcma_ll$loglike), 
                          LL = trcma_ll$loglike)
mcecm_ll_df <- data.frame(Algorithm = "Full MCECM", 
                          Iteration = 1:length(mcecm_ll$loglike[-1]),
                          LL = mcecm_ll$loglike[-1])
conv_ll_df <- rbind(trcma_ll_df, mcecm_ll_df)
p_conv <- ggplot(conv_ll_df) + 
  aes(x = Iteration, y = LL, color = Algorithm) +
  scale_x_continuous(breaks = seq(0,15,by=1), limits = c(1,15)) +
  labs(x = "Iteration", y = "Log-Likelihood", 
       title = "Convergence of Imputation Algorithms") +
  geom_line(size = 1) + geom_point(size = 2)
p_conv
saveRDS(p_conv, paste0(mypath, "em_convergence_plot.rds"))
ggsave(paste0(mypath, "em_convergence_plot.pdf"), p_conv, width = 7, height = 4.5)
