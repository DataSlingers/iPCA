# plot_ipca: given Sigma or Sigma estimate, plot ipca scores
#
# inputs:
#   Sig: n x n matrix; estimate of Sigma or true Sigma
#   U: n x r matrix; directly gives scores, e.g. output of CMF()
#   y: n-length vector of class labels 1, 2, 3, .... (factors)
#   plot_title: string of plot title
#   show_legend: T/F
#   point_size: size of pc points
#   text_size: text size
#   show_plot: T/F
#   pcs: e.g. c(1,2) = show pcs 1 and 2
# outputs:
#   ipca_plot: PC1 vs PC2
#   ipca_scores: all ipca scores PC1, ...., PCn

library(ggplot2)
library(RColorBrewer)
library(gridExtra)
# see library(ggpubr) # http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization

plot_ipca <- function(Sig, U, y, plot_title, show_legend = F, point_size = 1, text_size = 12, show_plot = T, pcs) {
  
  mytheme <- theme(legend.text = element_text(family = "Helvetica", size = rel(1)), 
                   axis.title = element_text(family = "Helvetica", size = rel(1)), 
                   axis.text = element_text(family = "Helvetica", size = rel(1)), 
                   axis.line = element_line(size = 1,colour = "black"), 
                   axis.ticks = element_line(colour="black",size = rel(1)),
                   panel.grid.major = element_line(colour="grey90",size = rel(0.5)), 
                   panel.grid.minor = element_blank(), 
                   panel.background = element_rect(fill = "grey98"), 
                   legend.key = element_rect(fill = "grey98"), 
                   legend.title = element_text(family = "Helvetica", size = rel(1)), 
                   plot.title = element_text(face = "bold", size = rel(1.25),family = "Helvetica"))
  
  # myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
  # sc <- scale_colour_gradientn(colours = myPalette(100)) # for continuous y labels
  # sc <- scale_colour_viridis(discrete = F, option = "plasma")
  sc <- scale_colour_viridis(discrete = F, option = "plasma", begin = 0, end = 0.95)
  
  if (!(missing(U))) {
    Sig_V <- U
  }else {
    Sig_eig <- eigen(Sig)
    Sig_V <- Sig_eig$vectors; Sig_d <- Sig_eig$values
  }
  
  if (!missing(pcs)) { # rearrange depending on desired pcs
    Sig_V <- cbind(Sig_V[,pcs], Sig_V[,-pcs]) 
  }
  
  if (!missing(y)) {
    p1 <- ggplot(data = data.frame(Sig_V, y = y)) +
      aes(x = X1, y = X2, color = y) +
      geom_point(size = point_size) +
      labs(x = "PC1", y = "PC2", color = "Class") +
      mytheme 
    if (!is.factor(y)) {
      p1 <- p1 + sc
    }
  }else {
    p1 <- ggplot(data = data.frame(Sig_V)) +
      aes(x = X1, y = X2) +
      geom_point(size = point_size) +
      labs(x = "PC1", y = "PC2", color = "Class") +
      mytheme 
  }
  
  if (!(missing(plot_title))) {
    p1 <- p1 + labs(title = plot_title)
  }
  
  if (show_legend == F) {
    p1 <- p1 + guides(color = F, shape = F)
  }
  
  if (show_plot) {
    print(p1)
  }
  
  return(list(ipca_plot = p1, ipca_scores = Sig_V))
}