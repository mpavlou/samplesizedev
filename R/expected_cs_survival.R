
#' Calculate the expected calibration slope and mean for a given  sample size
#'
#' @description
#' This function calculates the expected calibration slope given key data and model characteristics
#' (proportion of censoring at time point of interest, C-index and number of predictors).
#' @param n (numeric) The sample size
#' @param p (numeric) The anticipated proportion of events at the time-point of interest
#' @param c (numeric) The C-index
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000)
#' @param nval (numeric) Size of validation data
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#'
#' @return df: the expected calibration slope
#'
#' @examples
#' # expected_cs_survival(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = FALSE)
#'
#' # Prefer parallel computing with >2 cores that ensure faster running
#' # expected_cs_survival(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = TRUE)

#' @seealso
#' samplesizedev_survival

#'
#'
expected_cs_survival <- function(n, p, c, n.predictors, beta=rep(1/n.predictors, n.predictors), nsim = 1000, nval = 25000, parallel=TRUE, method="MLE"){

  # Find mean and variance of for Normal linear predictor

  set.seed(2022)
  p.censor <- 1 - p

  #n       <- 50000
  lambda  <- 1

  if (c<0.8) variance <-  find_sigma_quick(c)  else variance <-  find_sigma(c)

  # check
  # Size - big data

  #beta      <-  rep(1, n.predictors)
  f         <-  sqrt(variance/sum(beta^2))
  beta      <-  f * beta

  sigma    <- diag(1, length(beta))
  x        <- mvtnorm::rmvnorm(nval, rep(0, length(beta)), sigma = sigma)
  eta      <- x %*% as.vector(beta)


  u        <- stats::runif(nval)
  t        <- -log(u)/((lambda) * exp(eta) )

  censor   <- rep(1,nval)
  cutoff   <- stats::quantile(t, 1- p.censor)
  censor[t > cutoff]=0
  ptrue <- mean(censor) ;ptrue; cutoff

  data <- data.frame(t, censor, eta)


  if (parallel==TRUE) {
    cores <- parallel::detectCores()
    cl    <- parallel::makeCluster(cores[1])} else
    cl    <- parallel::makeCluster(2)

  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`
  `%do%`    <- foreach::`%do%`

  cs   <- NULL
  i    <- 0

  xval           <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)

  a<- foreach::foreach(i = 1:nsim, .packages=c('mvtnorm','RcppNumerical', 'ggplot2', 'survival', 'stats' )) %dopar% {

    set.seed(i)

    # Predictors and Outcome

    sigma                <- diag(1, length(beta))
    x                    <- mvtnorm::rmvnorm(n, rep(0, length(beta)), sigma = sigma)
    eta                  <- x %*% as.vector(beta)
    names(x)             <- paste("x", seq(1:n.predictors), sep = "")
    u                    <- stats::runif(n)
    t                    <- -log(u)/((lambda) * exp(eta) )
    censor               <- rep(1,n)
    cutoff               <- stats::quantile(t, 1- p.censor)
    censor[t > cutoff]   <- 0
    status               <- censor
    eventtime            <- t
    data                 <- data.frame(eventtime, status, x)

    # Formula

    xvars       <- paste("X", seq(1:n.predictors), sep = "")
    measurevar  <- "Surv(eventtime, status)"
    formula     <- stats::as.formula(paste(measurevar, paste(xvars, collapse=" + "), sep=" ~ "))
    fit         <- survival::coxph(formula, data=data, x=TRUE, y=TRUE)
    fit


    # validation

    nval               <- nval
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

    cs[i] <- survival::coxph(Surv(eventtime, status) ~ eta_est, data=dataval)$coef

    c(cs[i])

  }

  parallel::stopCluster(cl)


  b      <- matrix(unlist(a), byrow=TRUE, nrow=nsim)
  cs     <- b[,1]


  df        <- data.frame(cs)
  df        <- stats::na.omit(df)
  cs_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = cs), size=12) +
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("Expected CS = ",round(mean(cs,na.rm=TRUE)/0.0025)*0.0025)) +
    ggplot2::geom_vline( ggplot2::aes(xintercept = mean(cs, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
    ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
    ggplot2::xlab("Calibration Slope")

  if ( abs(mean(cs, na.rm=TRUE)- 0.9) > 0.005)   cs_plot <-  cs_plot + ggplot2::geom_vline( ggplot2::aes(xintercept = 0.9), color="red", linetype ="dashed", size = 1)

  print(cs_plot)
  set.seed(2022)

  # xval      <- mvtnorm::rmvnorm(200000, rep(0,n.predictors), sigma = sigma)
  # yval      <- stats::rbinom(200000, 1, invlogit(mean + xval %*% beta))
  # prev      <- mean(yval)
  # cstat     <- quickcstat(yval, invlogit(mean + xval %*% beta))

  df        <- data.frame(n, round(mean(cs, na.rm = TRUE)/0.0025) * 0.0025, round(sqrt(stats::var(cs,na.rm = TRUE)), 4), round(p, 2), round(c, 2 ) )
  names(df) <- c("N", "Expected CS", "SD(CS)", "Proportion of events", "C-index")

  #performance <- df[,-3]
  performance <- df


  performance

}


# system.time(a <- expected_cs_survival(n = 241, p = 0.1, c = 0.7, n.predictors = 12, nsim = 1000, nval = 10000, parallel = TRUE) )
# a
#
# system.time(a <- expected_cs_survival(n = 376, p = 0.5, c = 0.7, n.predictors = 12, nsim = 1000, nval = 10000, parallel = TRUE) )
# a
#
# system.time(a <- expected_cs_survival(n = 228, p = 0.5, c = 0.75, n.predictors = 12, nsim = 1000, nval = 10000, parallel = TRUE) )
# a
#
# system.time(a <- expected_cs_survival(n = 350, p = 0.5, c = 0.75, n.predictors = 12, nsim = 1000, nval = 10000, parallel = TRUE) )
# a



