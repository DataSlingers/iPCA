# plot_ipca_varexplained: function to plot proportion of variance explained for iPCA
#   covariance estimators over the specified number of components
#
# inputs:
#   Xcell: K-cell of n x p_k data matrices
#   Sig: n x n matrix; estimate of Sigma or true Sigma
#   Delts: K-cell of p_k x p_k matrices; estimate of Deltas or true Deltas
#   plot_title: string of plot title
#   show_legend: T/F
#   point_size: size of pc points
#   line_size: width of line
#   text_size: text size
#   show_plot: T/F
#   save_plot: T/F
#   mypath: main directory for saving
#   file_ext: file extension for saving
# outputs:
#   ipca_varexplained_plot: variation explained plot

library(ggplot2)
library(gridExtra)
library(ggpubr) # http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization

plot_ipca_varexplained <- function(Xs, Sig, Delts,
                                   plot_title = "iPCA Proportion of Variance Explained", 
                                   show_legend = T, point_size = 1.5, line_size = 1, 
                                   show_plot = T, save_plot = T, mypath, file_ext = "") {
  
  mytheme <- theme(legend.text = element_text(family = "Helvetica", size = rel(1.25)),
                   strip.text = element_text(family = "Helvetica", size = rel(1.25), face = "bold"),
                   axis.title = element_text(family = "Helvetica", size = rel(1.25)), 
                   axis.text = element_text(family = "Helvetica", size = rel(1.25)), 
                   axis.line = element_line(size = 1,colour = "black"), 
                   axis.ticks = element_line(colour="black",size = rel(1)),
                   panel.grid.major = element_line(colour="grey90",size = rel(0.5)), 
                   panel.grid.minor = element_blank(), 
                   panel.background = element_rect(fill = "grey98"), 
                   legend.key = element_rect(fill = "grey98"), 
                   legend.title = element_text(family = "Helvetica", size = rel(1.25)), 
                   plot.title = element_text(face = "bold", size = rel(1.75),family = "Helvetica"))
  
  K <- length(Xs); n <- nrow(Xs[[1]]); pks <- sapply(Xs, FUN = ncol)
  
  # center the data
  Xs <- lapply(X = Xs, FUN = function(X) {return(scale(X, center = T, scale = F))})
  
  
  Sig_eig <- eigen(Sig)
  U <- Sig_eig$vectors
  d <- Sig_eig$values
  Uinv <- solve(U)
  
  props_U <- props_Uk <- matrix(0, nrow = n, ncol = K)
  for (k in 1:K) {
    # regular PCA
    X_svd <- svd(Xs[[k]])
    Uk <- X_svd$u; Dk <- diag(X_svd$d); Vk <- X_svd$v
    
    Deltk_V <- eigen(Delts[[k]])$vectors
    
    tv <- (norm(Xs[[k]], 'F'))^2 # total variance
    
    # plot proportion of variance of X explained per matrix
    for (l in 1:min(n,pks[k])) {
      
      proj_U <- t(U[,1:l]) %*% Xs[[k]] %*% Deltk_V[,1:l]
      var_U <- norm(proj_U, 'F')^2
      
      prop_var_U <- var_U / tv
      props_U[l,k] <- prop_var_U

      proj_Uk <- Uk[,1:l] %*% t(Uk[,1:l]) # usual PCA
      var_Uk <- norm(proj_Uk %*% Xs[[k]], 'F')^2
      prop_var_Uk <- var_Uk / tv
      props_Uk[l,k] <- prop_var_Uk
      
    }
  }
  
  props_U_long <- cbind(reshape2::melt(data.frame(m = 1:n, props_U), id = "m"), var_ex = "iPCA (Joint PCs)")
  props_Uk_long <- cbind(reshape2::melt(data.frame(m = 1:n, props_Uk), id = "m"), var_ex = "PCA (Individual PCs)")
  
  df_long <- rbind(props_U_long, props_Uk_long)
  df_long[df_long == 0] <- NA
  
  p1 <- ggplot(df_long) +
    aes(x = m, y = value, group = var_ex) +
    facet_grid(~variable) +
    # geom_point(size = point_size, aes(color = var_ex)) +
    geom_line(size = line_size, aes(color = var_ex)) +
    labs(x = "Number of Components", y = "Proportion of Variance Explained", color = "") +
    labs(title = plot_title) +
    ylim(0-1e-5,1+1e-5) +
    mytheme
  
  if (show_legend == F) {
    p1 <- p1 + guides(color = F, shape = F)
  }
  
  if (show_plot) {
  #   annotate_figure(p1,
  #                   top = text_grob(plot_title, color = "black", face = "bold", size = 16))
    print(p1)
  }
  
  if (save_plot) {
    ggsave(paste0(mypath,"ipca_prop_var_explained", file_ext, ".pdf"), p1, width = 11, height = 6)
    saveRDS(p1, file = paste0(mypath,"ipca_prop_var_explained", file_ext, ".rds"))
  }
  
  return(list(ipca_varexplained_plot = p1, df = df_long))
}
