approximate_R2 <- function(auc, prev, n = 1000000){
  
  # define mu as a function of the C-statistic
  mu   <- sqrt(2) * stats::qnorm(auc)
  
  # sigmain <- sqrt(2)*qnorm(auc)
  #  mu<-0.5*(2*prev-1)*(sigmain^2)+log(prev/(1-prev))
  
  
  # simulate large sample linear prediction based on two normals 
  # for non-eventsN(0, 1), events and N(mu, 1)
  
  LP <- c(stats::rnorm(prev*n,  mean=0, sd=1), stats::rnorm((1-prev)*n, mean=mu, sd=1))
  y <- c(rep(0, prev*n), rep(1, (1-prev)*n))
  
  # Fit a logistic regression with LP as covariate;
  # this is essentially a calibration model, and the intercept and    
  # slope estimate will ensure the outcome proportion is accounted 
  # for, without changing C-statistic
  
  fit <- rms::lrm(y~LP)
  
  max_R2 <- function(prev){
    1-(prev^prev*(1-prev)^(1-prev))^2
  }
  return(list(R2.nagelkerke = as.numeric(fit$stats['R2']), 
              R2.coxsnell = as.numeric(fit$stats['R2']) * max_R2(prev)))   
}