# helper function for rf metric for iPCA (only works for categorical Y right now)
rf_metric <- function(dat, split, Y, trials = 1) {
  
  metric <- rep(NA, l = trials)
  
  for (trial in 1:trials) {
    Xtr = dat[split,]
    Xts = dat[-split,]
    Ytr = Y[split]
    Yts = Y[-split]
    
    rf_fit <- randomForest(x = Xtr, y = Ytr,
                           xtest = Xts, ytest = Yts)
    
    if (is.factor(Y)) {
      metric[trial] <- rf_fit$test$err.rate[length(rf_fit$test$err.rate)] # classification accuracy
      # con_mat <- rf_fit$test$confusion[,-ncol(rf_fit$test$confusion)]
      # metric[trial] <- sum(diag(con_mat)) / sum(sum(con_mat)) # classification accuracy; same as above
    }else {
      metric[trial] <- sqrt(rf_fit$test$mse[length(rf_fit$test$mse)]) # RMSE
    }
  }
  
  return(mean(metric))
}
