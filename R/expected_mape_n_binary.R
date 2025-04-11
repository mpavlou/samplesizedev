expected_mape_n_binary <- function(n, MAPE, mean_eta, variance_eta,  p, c, beta, n.predictors, nsim = 1000, nval = 25000, parallel = TRUE){

  set.seed(2022)

  # Find beta that corresponds to that variance

  # beta    <- rep(1, n.predictors)
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

  mape <- NULL
  i    <- 0


  a <- foreach::foreach(i=1:nsim,  .packages=c('mvtnorm','RcppNumerical', 'ggplot2' )) %dopar% {

    set.seed(i)

    invlogit<-function(x) 1/(1+exp(-x))

    x     <- mvtnorm::rmvnorm(round(n), rep(0, n.predictors), sigma = sigma )
    y     <- stats::rbinom(round(n),  1, invlogit(mean_eta + x%*%beta))
    yval  <- stats::rbinom(nval, 1, invlogit(mean_eta + xval%*%beta))

    p_true   <- as.vector(invlogit(mean_eta + xval%*%beta))

    #a       <- fastglm(cbind(1,x), y, family=binomial())
    a        <- RcppNumerical::fastLR(cbind(1,x), y)
    eta_est  <- cbind(1, xval) %*% as.vector(a$coef)
    p_est    <- as.vector(invlogit(eta_est))


    #fit      <- fastglm(cbind(1,eta_est), yval, family=binomial())
    fit      <- RcppNumerical::fastLR(cbind(1,eta_est), yval)
    mape[i]  <- mean(abs(p_true-p_est))
    mape[i]


  }

  parallel::stopCluster(cl)

  mape <- unlist(a)

  #graphics::hist(cs, main = paste("CS=",round(mean(cs,na.rm=TRUE)/0.0025)*0.0025, "N=",n))
  df        <- data.frame(mape)
  df        <- stats::na.omit(df)
  mape_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = mape), size=12) +
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("N = ", n,", Scaled MAPE = ", round(mean(mape/phi,na.rm=TRUE),4), ", Expected MAPE = ", round(mean(mape,na.rm=TRUE)/0.0025)*0.0025, ", SD(MAPE) = ", round(sqrt(stats::var(mape)),3))) +
    ggplot2::geom_vline( ggplot2::aes(xintercept = mean(mape, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
    ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
    ggplot2::xlab("MAPE")

  print(mape_plot)
  c(round(median(mape,na.rm=TRUE)/0.00025)*0.00025, sqrt(var(mape)/nsim))
  c(round(median(mape,na.rm=TRUE)/phi/1)*1, round(median(mape,na.rm=TRUE)/0.00025)*0.00025, sqrt(var(mape)/nsim))


}
