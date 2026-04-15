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
#' @param gamma (numeric) relative strength of predictor variables (same length as n_predictors)
#' @param beta (numeric) Strength of predictors (same length as n.predictors)
#' @param long (logical) Extract all simulations instead of just averages

#' @return a data frame df with elements:
#'             the input sample size
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
expected_cs_mape_binary <- function(n, p, c, n.predictors, beta, nsim = 1000, nval = 25000, method ="MLE", parallel = TRUE, long = FALSE, approx=FALSE, threshold, individual_predicted_probability){

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
  p_true          <- as.vector(invlogit(mean + x%*%beta))

  if (length(individual_predicted_probability)==0) individual_predicted_probability=median(p_true)

  # index of closest value
  idx                  <- which.min(abs( p_true - individual_predicted_probability))
  x_ipp                <- x[idx, ]
  p_ipp_true           <- as.vector(invlogit(mean + x_ipp%*%beta))

  individual_quantile  <- round(mean(p_true <= p_ipp_true, na.rm=TRUE), 2)


  # quantiles you want
  quant_grid <- seq(0.1, 0.9, by = 0.01)

  ord <- order(eta)
  eta_sorted <- eta[ord]
  x_sorted   <- x[ord, , drop = FALSE]

  eta_q <- quantile(eta, probs = quant_grid)

  idx <- findInterval(eta_q, eta_sorted)
  idx <- pmax(1, pmin(idx, length(eta_sorted)))

  x_quantile_all <- x_sorted[idx, , drop = FALSE]

  # corresponding true probabilities at those points
  p_quantile_true_all <- as.vector(invlogit(mean + x_quantile_all %*% beta))
#
#   i_individual    <- which.min(abs(quant_grid - individual_predicted_probability))
#
#   #reduce to one
#   x_quantile      <- x_quantile_all[ i_individual, , drop = FALSE]
#   p_quantile_true <- p_quantile_true_all[ i_individual]

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
  fit        <- glm(y~x, family=binomial())
  varbeta    <- vcov(fit)
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

  cs             <- NULL
  mape           <- NULL
  opt            <- NULL
  r2_app         <- NULL
  heuristic      <- NULL
  ave_pred_risk  <- NULL
  cest           <- NULL
  brier          <- NULL
  sens           <- NULL
  nb             <- NULL
  p_quantile     <- NULL
  p_quantile_all <- matrix(NA, nrow = nsim, ncol = length(quant_grid))

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

      x          <- mvtnorm::rmvnorm(round(n), rep(0, n.predictors), sigma = sigma )
      y          <- stats::rbinom(round(n),  1, invlogit(mean + x%*%beta))
      yval       <- stats::rbinom(nval, 1, invlogit(mean + xval%*%beta))
      p_true     <- as.vector(invlogit(mean + xval%*%beta))


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

      # Set decision threshold
      # threshold <- p
      # Classify predictions
      pred_class <- ifelse(p_est >= threshold, 1, 0)

      # Calculate TP, FP, and N
      TP <- sum(pred_class == 1 & yval == 1)
      FP <- sum(pred_class == 1 & yval == 0)

      fit              <- RcppNumerical::fastLR(cbind(1,eta_est), yval, start = c(0,0.9) )
      cs[i]            <- fit$coef[2]
      mape[i]          <- mean(abs(p_true-p_est))
      brier[i]         <- mean((p_est-yval)^2)
      opt[i]           <- r2_cs_app/MaxR2 - r2_cs_true/MaxR2
      heuristic[i]     <- 1 - n.predictors/LR
      r2_app[i]        <- r2_cs_app
      ave_pred_risk[i] <- mean(p_est)
      cest[i]          <- quickcstat(yval, p_est)
      nb[i]            <- (TP / nval) - (FP / nval) * (threshold / (1 - threshold))  #net benefit
      sens[i]          <- (TP / sum(yval))
      p_quantile[i]         <- as.vector(invlogit(c(1,x_ipp)%*%as.vector(a$coef)))
      p_quantile_all_i <- as.vector(invlogit(cbind(1, x_quantile_all) %*% as.vector(a$coef)))


      c(cs[i], mape[i], opt[i], heuristic[i], r2_app[i],  ave_pred_risk[i], cest[i], brier[i], sens[i], nb[i],   p_quantile[i],  p_quantile_all_i)


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
      brier[i]         <- mean((p_est-yval)^2)
      opt[i]           <- NA
      heuristic[i]     <- NA
      r2_app[i]        <- NA
      ave_pred_risk[i] <- mean(p_est)
      cest[i]          <- quickcstat(yval, p_est)
      p_quantile[i]         <- as.vector(invlogit(c(1,x_ipp)%*%betasf))
      p_quantile_all_i <- as.vector(invlogit(cbind(1, x_quantile_all) %*% as.vector(betasf)))

      # Set decision threshold
      # threshold <- p
      # Classify predictions
      pred_class <- ifelse(p_est >= threshold, 1, 0)

      # Calculate TP, FP, and N
      TP <- sum(pred_class == 1 & yval == 1)
      FP <- sum(pred_class == 1 & yval == 0)

      nb[i]           <- (TP / nval) - (FP / nval) * (threshold / (1 - threshold))  #net benefit
      sens[i]         <- (TP / sum(yval))


      c(cs[i], mape[i], opt[i], heuristic[i], r2_app[i],  ave_pred_risk[i], cest[i], brier[i], sens[i], nb[i],   p_quantile[i],  p_quantile_all_i)

    }
  }  else if (method == "ridge" | method=="lasso") {

    a<- foreach::foreach(i = 1:nsim,
                         .packages=c('mvtnorm','RcppNumerical', 'ggplot2','glmnet','foreach')) %dopar% {

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

      mod_penal_ave_foreach <- function(x, y, bn=50, method=method, f = nfolds/(nfolds-1)-1,
                                        parallel=TRUE, nfolds=10, boot=TRUE){

        if (method=="lasso") a=1
        if (method=="ridge") a=0

        pred        <- NULL ; beta.boot   <- NULL
        lambda.boot <- NULL ; cv.boot     <- NULL
        yval        <- y
        xval        <- x
        npred       <- ncol(x)
        data        <- cbind(x,y)

        fit <- glmnet(x, y, type.measure = "deviance", family = "binomial", alpha = a, standardize = TRUE, parallel = T)

        # Obtain the default sequence of lambda, and include some smaller values too
        lambdaseq <- fit$lambda
        lambdaseq <- unique(c(lambdaseq, seq( lambdaseq[length(lambdaseq)], lambdaseq[length(lambdaseq)]/20,
                                              -lambdaseq[length(lambdaseq)]/20 )))

        # cvd is the cross-validated deviance obtained in each of the 'over-sampled datasets'
        cvd <- NULL

        cv.boot <- foreach(i=1:bn, .combine='rbind', .packages=c('MASS','RcppNumerical', 'pROC', 'brglm2', 'detectseparation',
                                                                 'RcppNumerical', 'logistf', 'glmnet')) %dopar% {


                                                                   withwarnings <- function(expr) {

                                                                     val        <- NULL
                                                                     myWarnings <- NULL
                                                                     wHandler   <- function(w) {
                                                                       myWarnings <<- c(myWarnings, w$message)
                                                                       invokeRestart("muffleWarning")
                                                                     }

                                                                     myError  <- NULL
                                                                     eHandler <- function(e) {
                                                                       myError <<- e$message
                                                                       NULL
                                                                     }
                                                                     val <- tryCatch(withCallingHandlers(expr, warning = wHandler), error = eHandler)
                                                                     list(value = val, warnings = myWarnings, error=myError)
                                                                   }

                                                                   bs     <- sample(nrow(data), replace=T)
                                                                   databs <- data[bs,]

                                                                   if (boot==FALSE) databs<-data

                                                                   if (f==0) databs  <- data else
                                                                   {

                                                                     bs2     <- sample(nrow(data)*f, replace=T)
                                                                     databs2 <- data[bs2,]
                                                                     databs  <- rbind(databs,databs2)
                                                                   }

                                                                   xbs   <- databs[,1:npred]
                                                                   ybs   <- databs[,npred+1]

                                                                   fitbs <- withwarnings(cv.glmnet(xbs, ybs, type.measure = "deviance", family = "binomial", alpha = a, nfolds = nfolds,
                                                                                                   standardize = TRUE, parallel = parallel, lambda = lambdaseq))


                                                                   if ( length(fitbs$error) == 0 ) {
                                                                     fitbs  <- fitbs$value
                                                                     cvd    <- fitbs$cvm

                                                                   } else

                                                                   {
                                                                     cvd <- cvd
                                                                   }

                                                                   rbind(cvd)

                                                                 }

        b <- cbind(lambdaseq, colMeans(cv.boot, na.rm = TRUE))

        # Selected the value of lambda that minimises the cross-validated deviance

        lambda.boot <- b[order(b[,2]),] [1,1]


        # Obtain the solution using ridge/lasso at the selected value of lambda for the *original data*

        fit <- glmnet(x, y, type.measure = "deviance", family = "binomial", alpha = a, standardize = TRUE,
                      parallel = T, lambda = lambdaseq)

        beta.boot   <- as.vector(coef(fit, s = lambda.boot) )

        return(list("beta.boot" = beta.boot, "lambda.boot" = lambda.boot))
      }

      x        <- mvtnorm::rmvnorm(round(n), rep(0, n.predictors), sigma = sigma )
      y        <- stats::rbinom(round(n),  1, invlogit(mean + x%*%beta))
      yval     <- stats::rbinom(nval, 1, invlogit(mean + xval%*%beta))
      p_true   <- as.vector(invlogit(mean + xval%*%beta))

      pen       <-  mod_penal_ave_foreach(x=x, y=y, method=method, bn=5, nfolds=10)

      eta_est  <- cbind(1, xval) %*% as.vector(pen$beta.boot)
      p_est    <- as.vector(invlogit(eta_est))

      # Set decision threshold
      # threshold <- p
      # Classify predictions
      pred_class <- ifelse(p_est >= threshold, 1, 0)

      # Calculate TP, FP, and N
      TP <- sum(pred_class == 1 & yval == 1)
      FP <- sum(pred_class == 1 & yval == 0)

      fit              <- RcppNumerical::fastLR(cbind(1,eta_est), yval, start = c(0,0.9) )
      cs[i]            <- fit$coef[2]
      mape[i]          <- mean(abs(p_true-p_est))
      brier[i]         <- mean((p_est-yval)^2)
      opt[i]           <- NA
      heuristic[i]     <- NA
      r2_app[i]        <- NA
      ave_pred_risk[i] <- mean(p_est)
      cest[i]          <- quickcstat(yval, p_est)
      nb[i]           <- (TP / nval) - (FP / nval) * (threshold / (1 - threshold))  #net benefit
      sens[i]         <- (TP / sum(yval))
      p_quantile[i]         <- as.vector(invlogit(c(1,x_ipp)%*%as.vector(pen$beta.boot)))
      p_quantile_all_i <- as.vector(invlogit(cbind(1, x_quantile_all) %*% as.vector(pen$beta.boot)))


      c(cs[i], mape[i], opt[i], heuristic[i], r2_app[i],  ave_pred_risk[i], cest[i], brier[i], sens[i], nb[i],   p_quantile[i],  p_quantile_all_i)

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
  brier         <- b[,8]
  sens          <- b[,9]
  nb            <- b[,10]
  p_quantile    <- b[,11]

  res_mat <- do.call(rbind, a)

  p_quantile_all <- res_mat[, (12:ncol(res_mat))]


  df        <- data.frame(cs)
  df        <- stats::na.omit(df)
  cs_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = cs), size=12) +
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("Median CS = ",round(median(cs,na.rm=TRUE),2), ", Pr(0.85<CS<1.15) = ", round(mean(ifelse( (cs > 0.85 & cs <1.15), 1, 0),na.rm=TRUE), 2) )) +
    ggplot2::geom_vline( ggplot2::aes(xintercept = mean(cs, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
    ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
    ggplot2::xlab("Calibration Slope") + ggplot2::theme_bw()+  ggplot2::theme(legend.position="bottom")+
    ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10))+
    ggplot2::theme(plot.title =  ggplot2::element_text(size = 10)) +
    ggplot2::coord_cartesian(xlim = c(0.6, 1.6)) +
    ggplot2::scale_x_continuous(
      breaks = seq(0.6, 1.4, by = 0.1))

  if ( abs(mean(cs, na.rm=TRUE)- 0.9) > 0.005)   cs_plot <-  cs_plot + ggplot2::geom_vline( ggplot2::aes(xintercept = 0.9), color="red", linetype ="dashed", size = 1)


  df        <- data.frame(mape)
  df        <- stats::na.omit(df)
  mape_plot <- ggplot2::ggplot(df,  ggplot2::aes(x = mape), size=12) +
    ggplot2::geom_density() + ggplot2::ggtitle(paste("Median MAPE = ", round(median(mape,na.rm=TRUE),3), sep = "")) +
    ggplot2::geom_vline( ggplot2::aes(xintercept=mean(mape, na.rm = TRUE)), color="blue", linetype = "dashed", size=1) + ggplot2::ylab("Density") +
    ggplot2::theme(text =  ggplot2::element_text(size = 12)) +ggplot2::xlab("MAPE") + ggplot2::theme_bw()+
    ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10))+
    ggplot2::theme(plot.title =  ggplot2::element_text(size = 10))

  # Predicted Probability

  df        <- data.frame(p_quantile)
  df        <- stats::na.omit(df)

  # subset for shaded region
  df_shade <- df[df$p_quantile > threshold, ]

  prob_above <- mean(p_quantile > threshold, na.rm = TRUE)

  dens <- stats::density(p_quantile, na.rm = TRUE)
  df_dens <- data.frame(x = dens$x, y = dens$y)
  # df_shade <- df_dens[df_dens$x > threshold, ]

  # shading rule


  if (threshold < p_ipp_true) {
    # shade BELOW threshold
    prob_region <- mean(p_quantile < threshold, na.rm = TRUE)
    df_shade <- df_dens[df_dens$x < threshold, ]
    prob_label <- "P(IPP < threshold)"
  } else {
    # shade ABOVE threshold
    prob_region <- mean(p_quantile > threshold, na.rm = TRUE)
    df_shade <- df_dens[df_dens$x > threshold, ]
    prob_label <- "P(IPP > threshold)"
  }

  p_plot <- ggplot2::ggplot(df_dens, ggplot2::aes(x = x, y = y)) +

    # full density curve
    ggplot2::geom_line() +

    # shaded region (adaptive)
    ggplot2::geom_area(
      data = df_shade,
      fill = "red",
      alpha = 0.4
    ) +

    ggplot2::ggtitle(
      paste(
        "Sampling distribution of IPP=", round(p_ipp_true, 2),
        "\nMedian IPP = ", round(median(p_quantile, na.rm = TRUE), 2),
        "\n95% CI IPP = (",
        round(stats::quantile(p_quantile, probs = 0.025), 2), ", ",
        round(stats::quantile(p_quantile, probs = 0.975), 2), ") ; Width = ",
        round(
          stats::quantile(p_quantile, probs = 0.975) -
            stats::quantile(p_quantile, probs = 0.025), 2
        ),
        "\n", prob_label, " = ", round(prob_region, 2),
        sep = ""
      )
    ) +

    # median
    ggplot2::geom_vline(
      xintercept = median(p_quantile, na.rm = TRUE),
      color = "blue",
      linetype = "dashed",
      linewidth = 1
    ) +

    # threshold
    ggplot2::geom_vline(
      xintercept = threshold,
      color = "black",
      linetype = "dotted",
      linewidth = 1
    ) +
    ggplot2::ylab("Density") +
    ggplot2::xlab("Predicted Probability") +
    ggplot2::theme(text = ggplot2::element_text(size = 10)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(size = 10)
    ) +
    ggplot2::coord_cartesian(
      xlim = c(0, min(stats::quantile(p_quantile, probs = 0.999), 1))
    )

  # True probability distribution
  df        <- data.frame(p_true)
  df        <- stats::na.omit(df)
  p_true_plot <- ggplot2::ggplot(df,  ggplot2::aes(x = p_true), size=10) +
    ggplot2::geom_density() + ggplot2::ggtitle(paste("Distribution of true probabilities")) +
    ggplot2::geom_vline( ggplot2::aes(xintercept=  p_ipp_true), color="blue", linetype = "dashed", size=1) + ggplot2::ylab("Density") +
    ggplot2::theme(text =  ggplot2::element_text(size = 10)) +ggplot2::xlab("True Probabilities") +
    ggplot2::theme_bw()+ ggplot2::theme(legend.position="bottom")+
    ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10))+
    ggplot2::theme(plot.title =  ggplot2::element_text(size = 10))+
    ggplot2::coord_cartesian(
      xlim = c(
        max(0, stats::quantile(p_true, probs = 0.0001)),
        min(stats::quantile(p_true, probs = 0.9999), 1)
      )
    )

  # Get max density to scale placement nicely
  dens <- stats::density(p_true, na.rm = TRUE)
  df_dens <- data.frame(x = dens$x, y = dens$y)
  peak_y <- max(df_dens$y)

  y_anchor <- peak_y * 0.7
  y_offset <- peak_y * 0.12

  x_range <- max(df_dens$x) - min(df_dens$x)
  x_offset <- 0.15 * x_range

  p_true_plot <- p_true_plot +

    ggplot2::annotate(
      "segment",
      x = p_ipp_true + x_offset,
      y = y_anchor,
      xend = p_ipp_true,
      yend = y_anchor - y_offset,
      arrow = grid::arrow(length = grid::unit(0.25, "cm")),
      color = "red"
    ) +

    ggplot2::annotate(
      "text",
      x = p_ipp_true + x_offset,
      y = y_anchor + peak_y * 0.05,
      label = paste0(
        "IPP=", round(p_ipp_true, 3),
        "\n(", individual_quantile * 100, "- percentile)"
      ),
      color = "red",
      hjust = 0,
      size = 3
    )


  # All quantiles plot

  p_median <- apply(p_quantile_all, 2, median)
  p_q25    <- apply(p_quantile_all, 2, quantile, 0.25)
  p_q75    <- apply(p_quantile_all, 2, quantile, 0.75)

  df <- data.frame(
    q = quant_grid,
    med = p_median,
    lo = p_q25,
    hi = p_q75,
    true =   p_quantile_true_all
  )

  df <- stats::na.omit(df)

  p_plot_quantiles_all <- ggplot2::ggplot(df, ggplot2::aes(x = q, y = med)) +
    ggplot2::geom_line(color = "blue", linewidth = 1.2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      width = 0.005,
      color = "blue"
    ) +     ggplot2::theme_bw()+
    ggplot2::labs(x = "Percentile of distribution of true probs", y = "Median (95% CI) of IPPs") +
    ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10)) +
    ggplot2::theme(plot.title =  ggplot2::element_text(size = 10))


  df <- data.frame(
    q = quant_grid,
    med = p_median,
    lo = p_q25,
    hi = p_q75,
    true =   p_quantile_true_all
  )

  df <- stats::na.omit(df)

  p_plot_quantiles_all <- ggplot2::ggplot(df, ggplot2::aes(x = q)) +

    # median estimated curve
    ggplot2::geom_line(
      ggplot2::aes(y = med),
      color = "blue",
      linewidth = 1
    ) + ggplot2::geom_line(
      ggplot2::aes(y = true),
      color = "red",
      linewidth = 1,
      linetype = "dashed"
    ) +

    # vertical IQR bars
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lo, ymax = hi),
      width = 0.005,
      color = "blue",
      alpha = 0.5
    ) +
    ggplot2::labs(
      x = "Percentile of distribution of true probs",
      y = "Median (95% CI) of IPPs")+
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text=ggplot2::element_text(size=10),
                   axis.title=ggplot2::element_text(size=10)) +
    ggplot2::theme(plot.title =  ggplot2::element_text(size = 10))



  # ---- CS + MAPE ----
  figure1 <- ggpubr::ggarrange(
    cs_plot, mape_plot,
    ncol = 2, nrow = 1,
    common.legend = TRUE,
    legend = "bottom"
  )

  cs_mape_plot <- ggpubr::annotate_figure(
    figure1,
    top = ggpubr::text_grob(
      sprintf(
        "Sampling Distribution of CS and MAPE (nsims=%s) \nMethod=%s, N=%s, Prevalence=%s, C-stat=%s, No Predictors=%s",
        nsim, method, n, p, c, n.predictors
      ),
      color = "black", face = "bold", size = 10
    )
  )

  # ---- Probability plots ----
  figure2 <- ggpubr::ggarrange(
    p_true_plot, p_plot, p_plot_quantiles_all,
    ncol = 3, nrow = 1, widths= c(3,4,4)
  )

  prob_plot <- ggpubr::annotate_figure(
    figure2,
    top = ggpubr::text_grob(
      sprintf("\n Uncertainty/Stability of Individual Predicted Probabilities (IPP)"),
      color = "black", face = "bold", size = 10))

  final_plot <- ggpubr::ggarrange(
    cs_mape_plot, prob_plot,
    ncol = 1, nrow = 2
  )

  print(final_plot)

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
                          round(median(cest, na.rm = TRUE),3),
                          round(sqrt(stats::var(cest,na.rm = TRUE)), 4),
                          round(sqrt(var(ave_pred_risk, na.rm = TRUE)),3),
                          round(median(cs, na.rm = TRUE),3),
                          round(median(brier, na.rm = TRUE),3),
                          round(median(sens, na.rm = TRUE),3),
                          round(median(nb, na.rm = TRUE),3),
                          round(median(p_quantile, na.rm = TRUE),3),
                          round(sqrt(stats::var(p_quantile,na.rm = TRUE)), 3)
  )

  # round(mean(heuristic, na.rm = TRUE),3),
  # round(r2_cs_true,4),
  # round(mean(r2_app,na.rm=TRUE)*mean(cs, na.rm = TRUE),4),

  names(df) <- c("n","True prevalence", "True c-statistic", "Number of predictors","---------------------------",
                 "Mean_calibration_slope", "SD(CS)", "Pr(0.85<CS<1.15)", "Mean_MAPE",  "SD(MAPE)", "Optimism_R2_Nag",
                 "Mean_AUC", "SD(AUC)", "SD(Average Predicted Risk)", "Median CS", "Brier", "Sensitivity", "NB",
                 "Individual Predicted risk", "SD(IPP)")


  # names(df) <- c("n", "mean_CS", "sd_CS", "Pr(CS<0.8)", "mean_MAPE",  "sd_MAPE", "optimism_R2_Nag", "heuristic_SF", "r2_true", "r2_app/cs", "prevalence", "c-statistic", " # predictors")

  performance <- df[,-3]
  performance <- df

  options("scipen"=100, "digits"=4)

  if (long == FALSE)  t(performance) else

  {
    b        <- cbind(n, p, n.predictors, b)
    b        <- data.frame(b)
    names(b) <- c("n", "phi", "p", "cs", "mape", "opt_r2_nag", "heur_cs", "r2_apparent", "average_risk", "cstat", "brier", "sens", "nb", "  p_quantile")
    b
  }

}


# expected_cs_mape_binary(n = 530, phi = 0.2, c = 0.85, p=10, nsim = 2000, parallel = TRUE)
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, beta= c(0.5,0.3,0.2,0.1,0.1,
# rep(0,5)), nsim = 2000, parallel = TRUE)
# expected_cs_mape_binary(n = 530, p = 0.2, c = 0.85, n.predictors = 10, beta= c(0.9,0.1 ,0,0,0,
# rep(0,5)), nsim = 2000, parallel = TRUE)
