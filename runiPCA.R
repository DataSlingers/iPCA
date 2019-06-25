## this is a driver script to run iPCA on real data

library(R.utils)

# load in necessary R functions to run iPCA
sourceDirectory("./functions/", modifiedOnly = F, recursive = F)

# load in data; should be a list of K data matrices with the samples as rows and features as columns
dat <- # read in data here

# set lambdas to search over for penalty parameter selection in covariance estimation step; can modify lams as needed
lams <- c(1e-4, 1e-2, 1, 100, 10000)

# choose penalty parameters for covariance estimation step using CV-like method
# note: set greedy.search = T for computational speedup
best_lambdas <- choose_lambdas(dat = dat, q = "multfrob", lams = lams, 
                               greedy.search = T)$best_lambdas

# run iPCA with the multiplicative Frobenius estimator
iPCA_out <- FFmleMultFrob(dat = dat, lamDs = best_lambdas)

# plot top 2 iPCs/joint patterns from iPCA
# note: there is also an optional parameter y; use if there is a response of interest and want to color points in PC plots by y
iPCA_plot <- plot_ipca(Sig = iPCA_out$Sig, pcs = 1:2)
# plot_ipca(Sig = iPCA_out$Sig, pcs = 1:2, y = y)  # uncomment and use if y is provided

