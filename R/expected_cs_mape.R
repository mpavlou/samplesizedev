
#' Calculate the expected calibration slope and mean absolute prediction error (MAPE) for a given  sample size
#'
#' @description
#' This function calculates the expected calibration slope and MAPE given key data and model characteristics
#' (outcome prevalence, C-statistic and number of predictors). It takes approximately 15 seconds to test_allrun.
#' @param n (numeric) The sample size
#' @param p (numeric) The anticipated outcome prevalence
#' @param c (numeric) The C-statistic
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000)
#' @param nval (numeric) Size of validation data
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#'
#' @return df: the expected calibration slope and mape
#' @export
#'
#' @examples
#' # expected_cs_mape(n = 1000, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = FALSE)
#'
#' # Fast run with parallel computing
#' expected_cs_mape(n = 1000, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = TRUE)

#' @seealso
#' samplesizedev

#'
#'
expected_cs_mape <- function(n, p, c, n.predictors, nsim = 1000, nval = 25000, parallel=TRUE){

  # Find mean and variance of for Normal linear predictor

  set.seed(2022)

  mean_var     <- find_mu_sigma(p,c)
  mean         <- mean_var[1]
  variance     <- mean_var[2]

  # Find beta that corresponds to that variance

  beta    <- rep(1, n.predictors)
  beta    <- beta * sqrt(mean_var[2]/sum(beta^2))
  sigma   <- diag(1, n.predictors)

  xval    <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)


  if (parallel==TRUE) {
    cores <- parallel::detectCores()
    cl    <- parallel::makeCluster(cores[1]-2)} else
    cl    <- parallel::makeCluster(2)

  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  cs   <- NULL
  mape <- NULL
  i    <- 0

  a<- foreach::foreach(i = 1:nsim, .packages=c('mvtnorm','RcppNumerical', 'ggplot2' )) %dopar% {

    set.seed(i)
    invlogit <- function(x) 1/(1+exp(-x))

    x        <- mvtnorm::rmvnorm(round(n), rep(0, n.predictors), sigma = sigma )
    y        <- stats::rbinom(round(n),  1, invlogit(mean + x%*%beta))
    yval     <- stats::rbinom(nval, 1, invlogit(mean + xval%*%beta))
    p_true   <- as.vector(invlogit(mean + xval%*%beta))

    a        <- RcppNumerical::fastLR(cbind(1,x), y)
    eta_est  <- cbind(1, xval) %*% as.vector(a$coef)
    p_est    <- as.vector(invlogit(eta_est))

    fit      <- RcppNumerical::fastLR(cbind(1,eta_est), yval )
    cs[i]    <- fit$coef[2]
    mape[i]  <- mean(abs(p_true-p_est))

    c(cs[i],mape[i])

  }

  parallel::stopCluster(cl)


  b      <- matrix(unlist(a), byrow=TRUE, nrow=nsim)
  cs     <- b[,1]
  mape   <- b[,2]


  df        <- data.frame(cs)
  df        <- stats::na.omit(df)
  cs_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = cs), size=12) +
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("Expected CS = ",round(stats::median(cs,na.rm=TRUE)/0.005)*0.005)) +
    ggplot2::geom_vline( ggplot2::aes(xintercept=stats::median(cs, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) + ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13))


  df        <- data.frame(mape)
  df        <- stats::na.omit(df)
  mape_plot <- ggplot2::ggplot(df,  ggplot2::aes(x = mape), size=12) +
    ggplot2::geom_density() + ggplot2::ggtitle(paste("Expected MAPE = ", round(stats::median(mape,na.rm=TRUE),3), sep = "")) +
    ggplot2::geom_vline( ggplot2::aes(xintercept=stats::median(mape, na.rm = TRUE)), color="blue", linetype = "dashed", size=1) + ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13))

  figure  <- ggpubr::ggarrange(cs_plot, mape_plot,
                       ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

  cs_mape_plot <- ggpubr::annotate_figure(figure,
                                  top = ggpubr::text_grob(paste("Distribution of the Calibration Slope and MAPE\n","N=", n, ", Prevalence=", p, ", C-stat=",c,sep=""),
                                                  color = "black", face = "bold", size = 13))

  print(cs_mape_plot)

  set.seed(2022)

  xval      <- mvtnorm::rmvnorm(200000, rep(0,n.predictors), sigma = sigma)
  yval      <- stats::rbinom(200000, 1, invlogit(mean + xval %*% beta))
  prev      <- mean(yval)
  cstat     <- quickcstat(yval, invlogit(mean + xval %*% beta))

  df        <- data.frame(n, round(stats::median(cs, na.rm = TRUE)/0.005) * 0.005, round(sqrt(stats::var(cs,na.rm = TRUE)/nsim), 4), round(stats::median(mape, na.rm = TRUE),4), round(prev, 3), round(cstat, 3) )
  names(df) <- c("N", "Expected CS", "MCE(CS)", "Expected MAPE", "Prevalence", "C-Statistic")

  performance <- df[,-3]

  performance

}
