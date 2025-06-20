
#' Calculate the expected calibration slope and mean absolute prediction error (MAPE) for a given  sample size
#'
#' @description
#' This function calculates the expected calibration slope and MAPE given key data and model characteristics
#' (outcome prevalence, C-statistic and number of predictors). It takes approximately 15 seconds to run.
#' @param n (numeric) The sample size
#' @param p (numeric) The anticipated outcome prevalence
#' @param c (numeric) The C-statistic
#' @param n.predictors (numeric) The number of candidate predictor variables
#' @param nsim (numeric) The number of simulations (at least 500, default value 1000)
#' @param nval (numeric) Size of validation data
#' @param parallel (logical) parallel processing to speed up computations (default=TRUE)
#' @param method (character) the fitting method. "MLE" is the default and currently only option, but others will be added in future versions
#' @param parallel (numeric) relative strength of predictor variables (same length as n_predictors)
#' @param beta (numeric) Strength of predictors (same length as n.predictors)
#' @param long (logical) Extract all simulations instead of just averages
#' @param approx (logical) Extract all simulations instead of just averages



#'
#' @return a data frame df with elements:
#'             theinut sample size
#'             the expected calibration slope (mean_CS)
#'             the standard deviation of the CS (sd_CS)
#'             the probability of obtaining a miscalibrated model with calibration slope <0.8 (Pr(CS<0.8))
#'             the expected MAPE (MAPE)
#'             the standard deviation of the expected MAPE (sd_MAPE)
#'             the expected optimism in R square Nagelgerke (optimism_R2_Nag)
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
#'
expected_cs_mape_binary <- function(n, p, c, n.predictors, beta, nsim = 1000, nval = 25000, method ="MLE", parallel = TRUE, long = FALSE, approx=FALSE){

  # Find mean and variance of for Normal linear predictor
  # beta=rep(1/n.predictors, n.predictors)
  if (p>0.5) p  <- 1-p

  # set.seed(2022)

  mean_var     <- find_mu_sigma(p,c, tol = 0.00001)
  mean         <- mean_var[1]
  variance     <- mean_var[2]

  # Find beta that corresponds to that variance

  beta    <- beta * sqrt(mean_var[2]/sum(beta^2))
  sigma   <- diag(1, n.predictors)

  set.seed(2022)
  xval    <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)
  # yval     <- stats::rbinom(nval, 1, invlogit(mean + xval%*%beta))

  # True R2
  MaxR2      <- 1-(((p^(p))*((1-p)^(1-p)))^2)
  ncalc      <- 500000
  x          <- mvtnorm::rmvnorm(ncalc, rep(0, n.predictors), sigma = sigma )
  eta        <- mean+x%*% beta
  y          <- stats::rbinom(ncalc,  1, invlogit(eta))


  a          <- RcppNumerical::fastLR(cbind(1,x), y)
  L1         <- a$loglikelihood
  L0         <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
  LR         <- -2*(L0-L1)
  r2_cs_true <- 1 - exp(-LR/ncalc)
  r2_cs_true

  if (approx==TRUE) {

  ncalc      <- 100000
  x          <- mvtnorm::rmvnorm(ncalc, rep(0, n.predictors), sigma = sigma )
  eta        <- mean+x%*% beta
  y          <- stats::rbinom(ncalc,  1, invlogit(eta))
  fit <- glm(y~x, family=binomial())
  varbeta <- vcov(fit)
  sampbeta <- rmvnorm(nsim, c(mean, beta), sigma = varbeta*ncalc/n)
  }

  # data.calc <- data.frame(y,x)
  #
  # system.time(fit <- glm(y ~ ., data = data.calc, family = 'binomial'))
  #
  # LR      <- -2 * (as.numeric(logLik(glm(y ~ 1, data = data.calc,
  #                                        family = binomial(link = "logit")))) -
  #                    as.numeric(logLik(fit)))
  # r2_cs_true <- 1 - exp(-LR/ncalc)
  #
  # r2_cs_true
  #
  # system.time(fit <- glm(y~eta, data = data.calc, family="binomial"))
  #
  #
  # LR      <- -2 * (as.numeric(logLik(glm(y ~ 1, data = data.calc,
  #                                        family = binomial(link = "logit")))) -
  #                    as.numeric(logLik(fit)))
  # r2_cs_true <- 1 - exp(-LR/ncalc)
  #
  # r2_cs_true


  if (parallel==TRUE) {
    cores <- parallel::detectCores()
    cl    <- parallel::makeCluster(cores[1]-2)} else
      cl    <- parallel::makeCluster(2)

  doParallel::registerDoParallel(cl)

  `%dopar%` <- foreach::`%dopar%`
  `%do%` <- foreach::`%do%`

  cs            <- NULL
  mape          <- NULL
  opt           <- NULL
  r2_app        <- NULL
  heuristic     <- NULL
  ave_pred_risk <- NULL
  cest          <- NULL

  i    <- 0

  if (method== "MLE") {



    a<- foreach::foreach(i = 1:nsim, .packages=c('mvtnorm','RcppNumerical', 'ggplot2' )) %dopar% {

      set.seed(i)
      invlogit <- function(x) 1/(1+exp(-x))

      # Approximation of the C-statistic (Large n)

      quickcstat <- function(y, pred, seed=1){
        #set.seed(seed)
        casepred=pred[y == 1]
        conpred=pred[y == 0]

        if (length(conpred)>length(casepred)){
          conpred=conpred[sample(length(conpred),length(casepred),replace=FALSE)]
          auc.true=sum(casepred>conpred)/length(casepred)} else
          {
            casepred=casepred[sample(length(casepred),length(conpred),replace=FALSE)]
            auc.true=sum(casepred>conpred)/length(conpred)
          }

        return(auc.true)
      }

      x        <- mvtnorm::rmvnorm(round(n), rep(0, n.predictors), sigma = sigma )
      y        <- stats::rbinom(round(n),  1, invlogit(mean + x%*%beta))
      yval     <- stats::rbinom(nval, 1, invlogit(mean + xval%*%beta))
      p_true   <- as.vector(invlogit(mean + xval%*%beta))


      r2_cs_app <- NA

      a <- NULL
      if (approx==TRUE){
        a$coef <- sampbeta[i,]

      } else

        {a        <- RcppNumerical::fastLR(cbind(1,x), y)
        L1        <- a$loglikelihood
        # a0      <- RcppNumerical::fastLR(as.matrix(rep(1,n)), y)
        # L0      <- a0$loglikelihood
        L0        <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
        LR        <- -2*(L0-L1)
        r2_cs_app <- 1 - exp(-LR/n)
        }


      # a        <- RcppNumerical::fastLR(cbind(1,x), y)
      # L1        <- a$loglikelihood
      # # a0      <- RcppNumerical::fastLR(as.matrix(rep(1,n)), y)
      # # L0      <- a0$loglikelihood
      # L0        <- sum(y*log(mean(y)) + (1-y)*log(1-mean(y)))
      # LR        <- -2*(L0-L1)
      # r2_cs_app <- 1 - exp(-LR/n)


      eta_est  <- cbind(1, xval) %*% as.vector(a$coef)
      p_est    <- as.vector(invlogit(eta_est))


      fit              <- RcppNumerical::fastLR(cbind(1,eta_est), yval, start = c(0,0.9) )
      cs[i]            <- fit$coef[2]
      mape[i]          <- mean(abs(p_true-p_est))
      opt[i]           <- r2_cs_app/MaxR2 - r2_cs_true/MaxR2
      heuristic[i]     <- 1 - n.predictors/LR
      r2_app[i]        <- r2_cs_app
      ave_pred_risk[i] <- mean(p_est)
      cest[i]          <- quickcstat(yval, p_est)


      c(cs[i], mape[i], opt[i], heuristic[i], r2_app[i],  ave_pred_risk[i], cest[i])

    }

  } else if (method == "LSF") {

     quickcstat <- function(y, pred, seed=1){
        #set.seed(seed)
        casepred=pred[y == 1]
        conpred=pred[y == 0]

        if (length(conpred)>length(casepred)){
          conpred=conpred[sample(length(conpred),length(casepred),replace=FALSE)]
          auc.true=sum(casepred>conpred)/length(casepred)} else
          {
            casepred=casepred[sample(length(casepred),length(conpred),replace=FALSE)]
            auc.true=sum(casepred>conpred)/length(conpred)
          }

        return(auc.true)
      }

    bootsf <- function(data,n=100){
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
      opt[i]           <- NA
      heuristic[i]     <- NA
      r2_app[i]        <- NA
      ave_pred_risk[i] <- mean(p_est)
      cest[i]          <- quickcstat(yval, p_est)


      c(cs[i], mape[i], opt[i], heuristic[i], r2_app[i],  ave_pred_risk[i], cest[i])

    }
  }


  parallel::stopCluster(cl)


  b             <- matrix(unlist(a), byrow=TRUE, nrow=nsim)
  cs            <- b[,1]
  mape          <- b[,2]
  opt           <- b[,3]
  heuristic     <- b[,4]
  r2_app        <- b[,5]
  ave_pred_risk <- b[,6]
  cest          <- b[,7]


  df        <- data.frame(cs)
  df        <- stats::na.omit(df)
  cs_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = cs), size=12) +
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("Mean Calibration Slope = ",round(mean(cs,na.rm=TRUE),3))) +
    ggplot2::geom_vline( ggplot2::aes(xintercept = mean(cs, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
    ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
    ggplot2::xlab("Calibration Slope") + ggplot2::theme_bw()+  ggplot2::theme(legend.position="bottom")

  if ( abs(mean(cs, na.rm=TRUE)- 0.9) > 0.005)   cs_plot <-  cs_plot + ggplot2::geom_vline( ggplot2::aes(xintercept = 0.9), color="red", linetype ="dashed", size = 1)



  df        <- data.frame(mape)
  df        <- stats::na.omit(df)
  mape_plot <- ggplot2::ggplot(df,  ggplot2::aes(x = mape), size=12) +
    ggplot2::geom_density() + ggplot2::ggtitle(paste("Mean MAPE = ", round(mean(mape,na.rm=TRUE),3), sep = "")) +
    ggplot2::geom_vline( ggplot2::aes(xintercept=mean(mape, na.rm = TRUE)), color="blue", linetype = "dashed", size=1) + ggplot2::ylab("Density") +
    ggplot2::theme(text =  ggplot2::element_text(size = 13)) +ggplot2::xlab("MAPE") +
    ggplot2::theme_bw()+ ggplot2::theme(legend.position="bottom")

  figure  <- ggpubr::ggarrange(cs_plot, mape_plot,
                               ncol = 2, nrow = 1, common.legend = TRUE, legend="bottom")

  cs_mape_plot <- ggpubr::annotate_figure(figure,
                                          top = ggpubr::text_grob(paste("Distribution of the Calibration Slope and MAPE\n","N=", n, ", Prevalence=", p, ", C-stat=",c,sep=""),
                                                                  color = "black", face = "bold", size = 13)) + ggplot2::xlab("MAPE")

  print(cs_mape_plot)

  set.seed(2022)

  xval      <- mvtnorm::rmvnorm(500000, rep(0,n.predictors), sigma = sigma)
  yval      <- stats::rbinom(500000, 1, invlogit(mean + xval %*% beta))
  prev      <- mean(yval)
  cstat     <- quickcstat(yval, invlogit(mean + xval %*% beta))

  # A <- 2*p*(1-p)*qnorm(c)^2
  # app <- sqrt(1/(A*n)+2/(n-2) )

  options("scipen"=100, "digits"=4)

  df        <- data.frame(round(n),
                          round(prev, 2),
                          round(cstat, 2 ), n.predictors,
                          c(0),
                          round(mean(cs, na.rm = TRUE),3),
                          round(sqrt(stats::var(cs,na.rm = TRUE)), 4),
                          # round(sqrt( mean( ((cs-1)^2), na.rm=TRUE) ), 4),
                          round(mean(ifelse( (cs > 0.85 & cs <1.15), 1, 0),na.rm=TRUE), 2),
                          round(stats::median(mape, na.rm = TRUE),4),
                          round(sqrt(stats::var(mape,na.rm = TRUE)), 4),
                          round(mean(opt, na.rm = TRUE),3),
                          round(mean(cest, na.rm = TRUE),3),
                          round(sqrt(stats::var(cest,na.rm = TRUE)), 4),
                          round(sqrt(var(ave_pred_risk, na.rm = TRUE)),3),
                          round(median(cs, na.rm = TRUE),3))
                          # round(mean(heuristic, na.rm = TRUE),3),
                          # round(r2_cs_true,4),
                          # round(mean(r2_app,na.rm=TRUE)*mean(cs, na.rm = TRUE),4),

  names(df) <- c("n","True prevalence", "True c-statistic", "Number of predictors","---------------------------",  "Mean_calibration_slope", "SD(CS)", "Pr(0.85<CS<1.15)", "Mean_MAPE",  "SD(MAPE)", "Optimism_R2_Nag", "Mean_AUC", "SD(AUC)", "SD(Average Predicted Risk)", "Median CS")


  # names(df) <- c("n", "mean_CS", "sd_CS", "Pr(CS<0.8)", "mean_MAPE",  "sd_MAPE", "optimism_R2_Nag", "heuristic_SF", "r2_true", "r2_app/cs", "prevalence", "c-statistic", " # predictors")

  performance <- df[,-3]
  performance <- df

  options("scipen"=100, "digits"=4)

  if (long == FALSE)  t(performance) else

    {
    b        <- cbind(n, p, n.predictors, b)
    b        <- data.frame(b)
    names(b) <- c("n", "phi", "p", "cs", "mape", "opt_r2_nag", "heur_cs", "r2_apparent", "average_risk", "cstat")
    b
  }

}

# expected_cs_mape_binary(n = 530, phi = 0.2, c = 0.85, p=10, nsim = 2000, parallel = TRUE)
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, beta= c(0.5,0.3,0.2,0.1,0.1,
# rep(0,5)), nsim = 2000, parallel = TRUE)
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, beta= c(0.9,0.1 ,0,0,0,
# rep(0,5)), nsim = 2000, parallel = TRUE)


