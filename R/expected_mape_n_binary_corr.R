expected_mape_n_binary_corr <- function(n, MAPE, mean_eta, variance_eta,  p, c, beta, n.predictors, nsim = 1000, nval = 25000, parallel = TRUE, cor0=0, cor1=0){

  set.seed(2022)

  # Find beta that corresponds to that variance
  if (cor0 == 0 & cor1 == 0)

    sigma   <- diag(1, n.predictors) else {

      n.noise <- length(beta[beta==0])
      n.true  <- n.predictors-n.noise
      sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
      sigma[1:n.true, 1:n.true] <- cor0
      if (n.noise>0) {
        sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
      diag(sigma) <- 1
    }


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
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("N = ", n, ", Expected MAPE = ", round(mean(mape,na.rm=TRUE)/0.0001)*0.0001, ", SD(MAPE) = ", round(sqrt(stats::var(mape)),4))) +
    ggplot2::geom_vline( ggplot2::aes(xintercept = mean(mape, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
    ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
    ggplot2::xlab("MAPE")

  print(mape_plot)
  c(round(mean(mape,na.rm=TRUE)/0.0001)*0.0001, sqrt(var(mape)/nsim))

}
