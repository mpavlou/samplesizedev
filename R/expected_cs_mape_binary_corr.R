
#' Calculate the expected calibration slope and mean absolute prediction error (MAPE) for a given  sample size
#'
#' @description
#' This function calculates the expected calibration slope and MAPE given key data and model characteristics
#' (outcome prevalence, C-statistic and number of predictors). It takes approximately 15 seconds to test_allrun.
#' @param n (numeric) The sample size
#' @param p (numeric) The anticipated outcome prevalence
#' @param c (numeric) The C-statistic
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param cor0 (numeric) correlation between true predictors
#' @param cor1 (numeric) Correlation between noise predictors
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000)
#' @param nval (numeric) Size of validation data
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#' @param method (character) the fitting method. "MLE" is the default and currently only option, but others will be added in future versions
#' @param parallel (numeric) relative strength of predictor variables (same length as n_predictors)
#' @param beta (numeric) Strength of predictors (same length as n.predictors)

#'
#' @return df: the expected calibration slope and mape
#'
#' @examples
#' # expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = FALSE)
#'
#' # Prefer parallel computing with >2 cores that ensure faster running
#' # expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 100, parallel = TRUE)

#' @seealso
#' samplesizedev_binary samplesizedev_survival

#'
#'
expected_cs_mape_binary_corr <- function(n, p, c, beta = rep(1/n.predictors, n.predictors), n.predictors, nsim = 1000, nval = 25000, cor0=0, cor1=0, method ="MLE", parallel=TRUE){

  # Find mean and variance of for Normal linear predictor

  set.seed(2022)

  mean_var     <- find_mu_sigma(p, c, tol = 0.0001)
  mean         <- mean_var[1]
  variance     <- mean_var[2]

  # Find beta that corresponds to that variance

  if (cor0==0 & cor1 ==0) {

    beta    <- beta * sqrt(mean_var[2]/sum(beta^2))
    sigma   <- diag(1, n.predictors)} else

      {

    beta  <- adjust_multiplier_correlated(c=c, mean = mean, beta = beta, n.predictors = n.predictors, cor0=cor0, cor1=cor1)
    beta

      n.noise <- length(beta[beta==0])
      n.true  <- n.predictors-n.noise


      # Specify correlation matrix
        sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
        sigma[1:n.true, 1:n.true] <- cor0
        if (n.noise>0) {
        sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
        diag(sigma) <- 1

      }


  xval    <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)
  #yval     <- stats::rbinom(nval, 1, invlogit(mean + xval%*%beta))


  if (parallel==TRUE) {
    cores <- parallel::detectCores()
    cl    <- parallel::makeCluster(cores[1]-2)} else
    cl    <- parallel::makeCluster(2)

  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  cs        <- NULL
  mape      <- NULL
  cstat_est <- NULL
  i    <- 0

  if (method== "MLE") {

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

    fit      <- RcppNumerical::fastLR(cbind(1,eta_est), yval, start = c(0,0.9) )
    cs[i]    <- fit$coef[2]
    mape[i]  <- mean(abs(p_true-p_est))
    #NEW
     # cstat_est[i] <- quickcstat(yval, p_est)
    cstat_est[i] <- pROC::roc(as.vector(yval), as.vector(eta_est), quiet=TRUE)$auc

    c(cs[i],mape[i], cstat_est[i])

  }

  } else if (method == "LSF")

  {

    bootsf<-function(data,n=100){
      #first column outcome
      cal<-NULL
      for (j in 1:n){
        bs <- sample(nrow(data), replace=T)
        databs=data[bs,]
        xvarsbs=databs[,-1];ybs<-databs[,1]
        fitbs <- speedglm::speedglm(ybs~xvarsbs, family=binomial())

          eta_est    <- as.matrix(cbind(1,data[,-1]))%*%coef(fitbs)
          fitcal     <- speedglm::speedglm(data[,1]~eta_est, family=binomial())
          cal[j]     <- as.vector(stats::coef(fitcal)[2])

      }
      return(stats::median(cal,na.rm=TRUE))
    }


    a<- foreach::foreach(i = 1:nsim, .packages=c('mvtnorm','RcppNumerical', 'ggplot2', 'speedglm' )) %dopar% {

      set.seed(i)
      invlogit <- function(x) 1/(1+exp(-x))

      x        <- mvtnorm::rmvnorm(round(n), rep(0, n.predictors), sigma = sigma )
      y        <- stats::rbinom(round(n),  1, invlogit(mean + x%*%beta))
      yval     <- stats::rbinom(nval, 1, invlogit(mean + xval%*%beta))
      p_true   <- as.vector(invlogit(mean + xval%*%beta))

      a         <- RcppNumerical::fastLR(cbind(1,x), y)
      datasf    <- cbind(y, x)
      sf        <- bootsf(datasf, 100)
      betasf    <- c(1,rep(sf,n.predictors))*a$coef
      off       <- speedglm::speedglm(y~1,offset=cbind(1,x)%*%betasf,family=binomial())
      betasf[1] <- betasf[1]+stats::coef(off)
      eta_est   <- cbind(1, xval)%*%betasf
      p_est     <- as.vector(invlogit(eta_est))

      fit      <- RcppNumerical::fastLR(cbind(1,eta_est), yval )
      cs[i]    <- fit$coef[2]
      mape[i]  <- mean(abs(p_true-p_est))
      #NEW
      cstat_est[i] <-quickcstat(yval, p_est)

      c(cs[i],mape[i], cstat_est[i])
    }
    }


  parallel::stopCluster(cl)


  b      <- matrix(unlist(a), byrow=TRUE, nrow=nsim)
  cs     <- b[,1]
  mape   <- b[,2]
  #NEW
  cstat_est  <- b[,3]


  df        <- data.frame(cs)
  df        <- stats::na.omit(df)
  cs_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = cs), size=12) +
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("Mean Calibration Slope = ",round(mean(cs,na.rm=TRUE)/0.0025)*0.0025)) +
    ggplot2::geom_vline( ggplot2::aes(xintercept = mean(cs, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
    ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
    ggplot2::xlab("Calibration Slope")

  if ( abs(mean(cs, na.rm=TRUE)- 0.9) > 0.005)   cs_plot <-  cs_plot + ggplot2::geom_vline( ggplot2::aes(xintercept = 0.9), color="red", linetype ="dashed", size = 1)



  df        <- data.frame(mape)
  df        <- stats::na.omit(df)
  mape_plot <- ggplot2::ggplot(df,  ggplot2::aes(x = mape), size=12) +
    ggplot2::geom_density() + ggplot2::ggtitle(paste("Mean MAPE = ", round(mean(mape,na.rm=TRUE),3), sep = "")) +
    ggplot2::geom_vline( ggplot2::aes(xintercept=mean(mape, na.rm = TRUE)), color="blue", linetype = "dashed", size=1) + ggplot2::ylab("Density") +
    ggplot2::theme(text =  ggplot2::element_text(size = 13)) +     ggplot2::xlab("MAPE")

  figure  <- ggpubr::ggarrange(cs_plot, mape_plot,
                       ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

  cs_mape_plot <- ggpubr::annotate_figure(figure,
                                  top = ggpubr::text_grob(paste("Distribution of the Calibration Slope and MAPE\n","N=", n, ", Prevalence=", p, ", C-stat=",c,sep=""),
                                                  color = "black", face = "bold", size = 13)) + ggplot2::xlab("MAPE")

  print(cs_mape_plot)

  set.seed(2022)

  xval      <- mvtnorm::rmvnorm(2000000, rep(0,n.predictors), sigma = sigma)
  yval      <- stats::rbinom(2000000, 1, invlogit(mean + xval %*% beta))
  prev      <- mean(yval)
  cstat     <- quickcstat(yval, invlogit(mean + xval %*% beta))

  df        <- data.frame(n, ceiling(mean(cs, na.rm = TRUE)/0.0025) * 0.0025,
                             round(sqrt(stats::var(cs,na.rm = TRUE)), 4),
                             round(sqrt( mean( ((cs-1)^2), na.rm=TRUE) ), 4),
                             round(mean(ifelse( (cs < 0.8), 1, 0),na.rm=TRUE), 2),
                             round(mean(mape, na.rm = TRUE),4),
                             round(sqrt(stats::var(mape,na.rm = TRUE)), 4),
                             round(prev, 4),
                             round(cstat, 4 ),
                             n.predictors)
  names(df) <- c("N", "Mean_CS", "SD_CS", "RMSD_CS", "Pr(CS<0.8)", "Mean_MAPE",  "SD_MAPE", "Prev.", "C-Stat.", " # Predictors")

  performance <- df[,-3]
  performance <- df

  performance

  #long format for simulations
  performance <- data.frame(n, n.predictors,  round(cstat,3), round(prev, 3), cs, mape, cstat_est)
  names(performance) = c("n", "npred", "cstat","prev","cs_est","mape_est","cstat_est")
  performance

}

#expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, beta = c(0.5,0.3,0.3,0.15,0.15, rep(0,7)), n.predictors = 12, nsim = 2000, parallel = TRUE)


#expected_cs_mape_binary_corr(n = 530, p = 0.2, c = 0.85, beta = c(0.5,0.3,0.3,0.15,0.15, rep(0,7)), n.predictors = 12, nsim = 2000, parallel = TRUE, cor0=0.5, cor1=0.05)


# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 2000, parallel = TRUE)
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, beta= c(0.5,0.3,0.2,0.1,0.1,
# rep(0,5)), nsim = 2000, parallel = TRUE)
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, beta= c(0.9,0.1 ,0,0,0,
# rep(0,5)), nsim = 2000, parallel = TRUE)



