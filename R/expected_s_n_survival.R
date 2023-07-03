expected_s_n_survival <- function(n, S, variance_eta,  p, c, n.predictors, beta, nsim = 1000, nval = 25000, parallel = TRUE){

  set.seed(2022)


  # Find beta that corresponds to that variance

  #beta     <- rep(1, n.predictors)
  beta     <- beta * sqrt(variance_eta/sum(beta^2))
  beta     <- as.vector(beta)
  sigma    <- diag(1, n.predictors)
  lambda   <- 1
  p.censor <- 1 - p

  xval    <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)


  if (parallel==TRUE) {
       cores <- parallel::detectCores()
       cl    <- parallel::makeCluster(cores[1])} else
       cl    <- parallel::makeCluster(2)

  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`
  `%do%`    <- foreach::`%do%`

  cs   <- NULL
  i    <- 0



a <- foreach::foreach(i = 1: nsim, .packages = c('mvtnorm','RcppNumerical', 'ggplot2', 'survival')) %dopar% {

    set.seed(i)

    # Predictors and Outcome

    sigma               <- diag(1, n.predictors)
    x                   <- mvtnorm::rmvnorm(n, rep(0, n.predictors), sigma = sigma)
    eta                 <- x %*% as.vector(beta)
    names(x)            <- paste("X", seq(1:n.predictors), sep = "")
    u                   <- stats::runif(n)
    t                   <- -log(u)/((lambda) * exp(eta) )
    censor              <- rep(1,n)
    cutoff              <- stats::quantile(t, 1 - p.censor)
    censor[t > cutoff]  <- 0
    status              <- censor
    eventtime           <- t

    data                <- data.frame(eventtime, status, x)

    # Formula

    xvars              <- paste("X", seq(1 : n.predictors), sep = "")
    measurevar         <- "survival::Surv(eventtime, status)"
    formula            <- stats::as.formula(paste(measurevar, paste(xvars, collapse=" + "), sep=" ~ "))
    fit                <- survival::coxph(formula, data = data, x = TRUE, y = TRUE)


    # Validation

    xval               <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)
    eta_val            <- xval %*% as.vector(beta)
    names(xval)        <- paste("x", seq(1:n.predictors), sep = "")
    u                  <- stats::runif(nval)
    t                  <- -log(u)/((lambda) * exp(eta_val) )
    censor             <- rep(1,nval)
    censor[t > cutoff] <- 0
    status             <- censor
    eventtime          <- t
    eta_est            <- as.matrix(xval) %*% as.vector(fit$coef)
    dataval            <- data.frame(eventtime, status, eta_est)

    # Performance measures

    cs[i]          <- survival::coxph(Surv(eventtime, status) ~ eta_est, data=dataval)$coef

    c(cs[i])

  }

  parallel::stopCluster(cl)

 cs <- unlist(a)


 df        <- data.frame(cs)
 df        <- stats::na.omit(df)
 cs_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = cs), size=12) +
   ggplot2::geom_density() +  ggplot2::ggtitle(paste("N = ", n, ", Expected CS = ", round(mean(cs,na.rm=TRUE)/0.0025)*0.0025, ", SD(CS) = ", round(sqrt(stats::var(cs)),3))) +
   ggplot2::geom_vline( ggplot2::aes(xintercept = mean(cs, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
   ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
   ggplot2::xlab("Calibration Slope")

 if ( abs(mean(cs, na.rm=TRUE)- 0.9) > 0.005)   cs_plot <-  cs_plot + ggplot2::geom_vline( ggplot2::aes(xintercept = 0.9), color="red", linetype ="dashed", size = 1)
print(cs_plot)

  # graphics::hist(cs, main = paste("CS=", round( mean(cs, na.rm = TRUE) / 0.0025) * 0.0025, "N=",n))
   c(round(mean(cs, na.rm = TRUE) / 0.0025) * 0.0025, sqrt(stats::var(cs) / nsim))


}
