expected_s_n <- function(n, S, mean_eta, variance_eta,  p, c, n.predictors, nsim = 1000, nval = 20000, parallel = TRUE){
  
  set.seed(2022)
  
  # Find beta that corresponds to that variance
  
  beta    <- rep(1, n.predictors)
  beta    <- beta * sqrt(variance_eta/sum(beta^2))
  sigma   <- diag(1, n.predictors)
  
  xval    <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)
  
  
  if (parallel==TRUE) {
       cores <- parallel::detectCores()
       cl    <- parallel::makeCluster(cores[1]-2)} else
       cl    <- parallel::makeCluster(2)
  
  doParallel::registerDoParallel(cl)
  
  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`
  
  cs <- NULL
  i  <- 0
  

  a <- foreach::foreach(i=1:nsim,  .packages=c('mvtnorm','RcppNumerical', 'ggplot2' )) %dopar% {
    
    set.seed(i)
    
    invlogit<-function(x) 1/(1+exp(-x))
    
    x     <- mvtnorm::rmvnorm(round(n), rep(0, n.predictors), sigma = sigma )
    y     <- stats::rbinom(round(n),  1, invlogit(mean_eta + x%*%beta))
    yval  <- stats::rbinom(nval, 1, invlogit(mean_eta + xval%*%beta))
    
    
    #a       <- fastglm(cbind(1,x), y, family=binomial())
    a        <- RcppNumerical::fastLR(cbind(1,x), y)
    eta_est  <- cbind(1, xval) %*% as.vector(a$coef)
    
    #fit      <- fastglm(cbind(1,eta_est), yval, family=binomial()) 
    fit      <- RcppNumerical::fastLR(cbind(1,eta_est), yval)
    cs[i]    <- fit$coef[2]
    cs[i]
  
  }
  
  parallel::stopCluster(cl)

  cs <- unlist(a)

  graphics::hist(cs, main = paste("CS=",round(mean(cs,na.rm=TRUE)/0.0025)*0.0025, "N=",n))
  c(round(mean(cs,na.rm=TRUE)/0.0025)*0.0025, sqrt(stats::var(cs)/nsim))
  
 
}