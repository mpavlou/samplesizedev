#  Functions for the simulations in Paper
# 'Evaluation of sample size requirements for the development of prediction
#  models for binary outcomes'

#################################################################################################################

# Function to calculate the expected calibration slope and MAPE

# Inputs/Parameters

# n (numeric) The sample size
# p (numeric) The anticipated outcome prevalence
# c (numeric) The C-statistic
# n.predictors (numeric) The number of candidate predictor variables
# cor0 (numeric) correlation between true predictors
# cor1 (numeric) Correlation between noise predictors
# nsim (numeric) The number of simulations (at least 500, default value 1000)
# nval (numeric) Size of validation data
# parallel (logical) parallel processing to speed up computations (default=TRUE)
# method (character) the fitting method. "MLE" is the default and currently only option, but others will be added in future versions
# parallel (numeric) relative strength of predictor variables (same length as n_predictors)
# beta (numeric) Strength of predictors (same length as n.predictors)

# Output:
# df: the expected calibration slope and MAPE

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


  # Short format for presentation

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

  # Long format for simulations
  performance <- data.frame(n, n.predictors,  round(cstat,3), round(prev, 3), cs, mape, cstat_est)
  names(performance) = c("n", "npred", "cstat","prev","cs_est","mape_est","cstat_est")
  performance

}

# Examples

# expected_cs_mape_binary_corr(n = 530, p = 0.2, c = 0.85, n.predictors = 10, nsim = 1000, parallel = TRUE)

#################################################################################################################


# Function to calculate Sample size required to develop a risk prediction model based on
# CALIBRATION SLOPE

# Inputs
# S (numeric)      : The target expected calibration slope
# p (numeric)      : The anticipated outcome prevalence
# c (numeric)      : The anticipated C-statistic
# n.predictors     :(numeric) The number of candidate predictor variables
# nsim (numeric)   : The number of simulations (at least 500, default value 1000 to ensure small simulation error)
# nval (numeric)   : Size of validation data (at least 10000 )
# parallel(logical):  parallel processing to speed up computations (default=TRUE)

# Outputs: The required sample size to achieve required calibration slope
# (Actual by Simulation and RvS-2 for comparison)


samplesizedev_binary_s_corr <- function(S, p, c, n.predictors, beta=rep(1/n.predictors, n.predictors), nval = 25000, nsim = 1000, parallel = TRUE, cor0=0.1, cor1=0.05){

  set.seed(2022)

  mean_var     <- find_mu_sigma(p, c, tol = 0.00001)
  mean_eta         <- mean_var[1]
  variance_eta     <- mean_var[2]

  # Find beta that corresponds to that variance

  if (cor0==0 & cor1 ==0) {

    betan    <- beta * sqrt(mean_var[2]/sum(beta^2))
    sigma   <- diag(1, n.predictors)} else

    {

      betan  <- adjust_multiplier_correlated(c=c, mean =  mean_eta, beta = beta, n.predictors = n.predictors, cor0=cor0, cor1=cor1)

      n.noise <- length(beta[beta==0])
      n.true  <- n.predictors-n.noise


  # Specify correlation matrix
      sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
      sigma[1:n.true, 1:n.true] <- cor0
      if (n.noise>0) {
        sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
      diag(sigma) <- 1

    }


  r2   <- as.numeric(approximate_R2(c, p, n = 200000)[2])

  n_init <- round((n.predictors)/ ((S-1)*log(1-r2/S)))


  if (c<=0.7)            {inflation_f   <- 1.3 ; min.opt  <- n_init*0.7}
  if (c>0.7  & c<=0.8)   {inflation_f   <- 1.6  ; min.opt <- n_init}
  if (c>0.8  & c<=0.85)  {inflation_f   <- 2.1    ; min.opt <- n_init}
  if (c>0.85 & c<=0.9)   {inflation_f   <- 2.8  ; min.opt <- n_init}

  max.opt  <- inflation_f*n_init

  tol = ceiling(round(n_init/200)/5) * 5

  print("Optimisation Starting ~ 1 min left...")
  s_est <- function(n, nsim=nsim){

    s <-  expected_s_n_binary_corr(n, S = S, mean_eta = mean_eta, variance_eta = variance_eta,  beta = betan, p = p, c = c, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel=parallel, cor0=cor0, cor1=cor1)
    s[1] - S
  }

  n <- bisection(s_est, min.opt, max.opt, tol = tol, nsim = nsim)

  tol = ceiling(round(n/200)/5) * 5
  n <- ceiling(n/tol)*tol

  size        <- NULL
  size$rvs1   <- as.vector(n_init)
  size$actual <- as.vector(round(n))
  size$correct_to_nearest <- as.vector(tol)

  size

}

#################################################################################################################


# Function to calculated Sample size required to develop a risk prediction model based on MAPE

# Inputs
# MAPE (numeric)    : The target expected MAPE
# p (numeric)      : The anticipated outcome prevalence
# c (numeric)      : The anticipated C-statistic
# n.predictors     :(numeric) The number of candidate predictor variables
# nsim (numeric)   : The number of simulations (at least 500, default value 1000 to ensure small simulation error)
# nval (numeric)   : Size of validation data (at least 10000 )
# parallel(logical):  parallel processing to speed up computations (default=TRUE)

# Outputs: The required sample size to achieve required MAPE (Actual by Simulation and RvS-2 for comparison)


samplesizedev_binary_mape_corr <- function(MAPE, p, c,  n.predictors, beta, nval = 25000, nsim = 1000, parallel = TRUE, cor0=0, cor1=0){

  set.seed(2022)

  mean_var         <- find_mu_sigma(p,c)
  mean_eta         <- mean_var[1]
  variance_eta     <- mean_var[2]

  # Find beta that corresponds to that variance

  if (cor0==0 & cor1 ==0) {

    betan    <- beta * sqrt(mean_var[2]/sum(beta^2))
    sigma   <- diag(1, n.predictors)} else

    {

      betan  <- adjust_multiplier_correlated(c=c, mean = mean_eta, beta = beta, n.predictors = n.predictors, cor0=cor0, cor1=cor1)

      n.noise <- length(beta[beta==0])
      n.true  <- n.predictors-n.noise


  # Specify correlation matrix
      sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
      sigma[1:n.true, 1:n.true] <- cor0
      if (n.noise>0) {
        sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
      diag(sigma) <- 1

    }

  n_init <- exp((-0.508 + 0.259 * log(p) + 0.504 * log(n.predictors) - log(MAPE))/0.544) ;

  min.opt = round(n_init*0.5)
  max.opt = round(n_init*1.5)


  tol = ceiling(round(n_init/200)/5) * 5

  print("Optimisation Starting ~ 1 min left...")

  mape_est <- function(n, nsim=nsim){

    mape <-  expected_mape_n_binary_corr(n, MAPE = MAPE, mean_eta = mean_eta, variance_eta = variance_eta, beta=betan,  p = p, c = c, n.predictors = n.predictors, nval = nval, nsim = nsim, parallel = parallel, cor0=cor0)

    MAPE-mape[1]

  }

  n <- bisection_mape(mape_est, MAPE=MAPE, min.opt, max.opt, tol = tol, nsim = nsim)
  tol = ceiling(round(n/200)/5) * 5
  n <- ceiling(n/tol)*tol

  size        <- NULL
  size$rvs2   <- as.vector(round(n_init))
  size$actual <- as.vector(round(n))

  size

}


#######################################################################################################

# Function to calculate the expected calibration slope  for a given  sample size
# Feed into the function for sample size calculations based on shrinkage (calibration slope)
# Inputs same as above
# Outputs: mean calibration slope

expected_s_n_binary_corr <- function(n, S, mean_eta, variance_eta,  p, c, n.predictors, beta, nsim = 1000, nval = 25000, parallel = TRUE, cor0, cor1){

  set.seed(2022)

  # Specify correlation matrix

  if (cor0==0 & cor1 ==0) {
    sigma   <- diag(1, n.predictors)} else

    {

      n.noise <- length(beta[beta==0])
      n.true  <- n.predictors-n.noise

      sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
      sigma[1:n.true, 1:n.true] <- cor0
      if (n.noise>0) {
        sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
      diag(sigma) <- 1

    }

  xval    <- mvtnorm::rmvnorm(nval, rep(0, n.predictors), sigma = sigma)
  #yval  <- stats::rbinom(nval, 1, invlogit(mean_eta + xval%*%beta))


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

  #graphics::hist(cs, main = paste("CS=",round(mean(cs,na.rm=TRUE)/0.0025)*0.0025, "N=",n))

  yval  <- stats::rbinom(nval, 1, invlogit(mean_eta + xval%*%beta))


  df        <- data.frame(cs)
  df        <- stats::na.omit(df)
  cs_plot   <- ggplot2:: ggplot(df,  ggplot2::aes(x = cs), size=12) +
    ggplot2::geom_density() +  ggplot2::ggtitle(paste("N = ", round(n), "p = ", mean(yval), ", Expected CS = ", round(mean(cs,na.rm=TRUE)/0.0025)*0.0025, ", SD(CS) = ", round(sqrt(stats::var(cs,na.rm=TRUE)),3))) +
    ggplot2::geom_vline( ggplot2::aes(xintercept = mean(cs, na.rm = TRUE)), color="blue", linetype ="dashed", size = 1) +
    ggplot2::ylab("Density") +  ggplot2::theme(text =  ggplot2::element_text(size = 13)) +
    ggplot2::xlab("Calibration Slope")

  if ( abs(mean(cs, na.rm=TRUE)- 0.9) > 0.005)   cs_plot <-  cs_plot + ggplot2::geom_vline( ggplot2::aes(xintercept = 0.9), color="red", linetype ="dashed", size = 1)
  print(cs_plot)
  c(round(mean(cs,na.rm=TRUE)/0.0025)*0.0025, sqrt(stats::var(cs)/nsim))


}


#######################################################################################################


# Function to calculate the expected calibration slope  for a given  sample size
# Feed into the function for sample size calculations based on MAPE
# Inputs same as above
# Outputs: mean MAPE

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


#######################################################################################################

# Bisection method for sample size calculation (shrinkage)
# a, b: starting values
# f: funciton evaluated

bisection <- function(f, a, b, iter = 10, tol = ceiling(round(a/200)/5) * 5, nsim = 1000) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.

  tol = ceiling(round(a/200)/5) * 5

  nsim1 <- nsim

  if (nsim>=1000) divide <- 2 else divide <- 1

  nsim  <- round(nsim/divide)

  fa <- f(a, nsim = nsim)
  fb <- f(b, nsim = nsim)
  if (!(fa < 0) && (fb > 0)) {
    stop('signs of f(a) and f(b) differ')
  } else if ((fa > 0) && (fb < 0)) {
    stop('signs of f(a) and f(b) differ')
  }

  for (k in 1:iter) {

    if (k <= divide ) nsim <- nsim1/divide*k
    # a <- round(a/tol)*tol
    # b <- round(b/tol)*tol

    c  <- (a + b) / 2 # Calculate midpoint
    fc <- f(c, nsim = nsim)

    #print(c(k, a, b, c, fa, fb, fc ))
    #print(c(k))


    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
    if (  ((abs(fc) <= 0.0025) || ((b - a) / 2) < tol)   &  (k >=2 )) {
      #if (  abs(fc) <= 0.0025 &  (k >=2 )) {

      return(c)
    }

    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(fc) == sign(fa),
           a <- c,
           b <- c)

    ifelse(sign(fc) == sign(fa),
           fa <- fc,
           fb <- fc)

  }
  # If the max number of iterations is reached and no root has been found,
  # return message and end function.
  print('Too many iterations')
}

#############################################################################


# Bisection method for sample size calculation (MAPE)
# a, b: starting values
# f: funciton evaluated

bisection_mape <- function(f, a, b, MAPE = 0.0001, iter = 10, tol = ceiling(round(a/200)/5) * 5, nsim = 1000) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.

  tol = ceiling(round(a/200)/5) * 5

  nsim1 <- nsim

  if (nsim>=1000) divide <- 2 else divide <- 1
  nsim  <- round(nsim/divide)

  fa <- f(a, nsim = nsim)
  fb <- f(b, nsim = nsim)
  if (!(fa < 0) && (fb > 0)) {
    stop('signs of f(a) and f(b) differ')
  } else if ((fa > 0) && (fb < 0)) {
    stop('signs of f(a) and f(b) differ')
  }

  for (k in 1:iter) {

    if (k <= divide ) nsim <- nsim1/divide*k
    # a <- round(a/tol)*tol
    # b <- round(b/tol)*tol

    c  <- (a + b) / 2 # Calculate midpoint
    fc <- f(c, nsim = nsim)

    #print(c(k, a, b, c, fa, fb, fc ))
    #print(c(k))


    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
    if (  ((abs(fc) <= MAPE/200) || ((b - a) / 2) < tol)   &  (k >=2 )) {
      # if (  abs(fc) <= MAPE/200 &  (k >=2 )) {
      return(c)
    }

    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(fc) == sign(fa),
           a <- c,
           b <- c)

    ifelse(sign(fc) == sign(fa),
           fa <- fc,
           fb <- fc)

  }
  # If the max number of iterations is reached and no root has been found,
  # return message and end function.
  print('Too many iterations')
}


########################################################################################

# Find mean and variance of linear predictor for given prevalence
# and C-statistic (using numerical integration)
# Part of functions calculating expected shrinkage and MAPE


invlogit <- function(x) 1/(1+exp(-x))


find_mu_sigma <- function(target.prev, target.c, min.opt = c(-10,0), max.opt = c(0,5), tol = 0.00001){

  pcfun <- function(x){

    mean      <- x[1]
    variance  <- x[2]

    f1 = function(x) {
      stats::integrate(function(y) {stats::dnorm(x, mean = mean, sd = sqrt(variance)) * stats::dnorm(y, mean = mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1) },
                       -Inf, x)$value
    }

    num = stats::integrate(Vectorize(f1), -Inf, Inf)$value

    f2 = function(x) {
      stats::integrate(function(y) {stats::dnorm(x, mean = mean, sd = sqrt(variance)) * stats::dnorm(y, mean = mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1) * (1  + exp(y)) ^ (-1) },
                       -Inf, Inf)$value
    }

    denom <- stats::integrate(Vectorize(f2), -Inf, Inf)$value

    c      <- num/denom
    f3     <- function(x) stats::dnorm(x, mean=mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1)

    prev   <- stats::integrate(f3, - Inf, Inf, subdivisions = 1000L)$value
    abs( c - target.c)^2 + abs(prev - target.prev )^2
  }

  out      <- stats::optim(par=c(-2.65,0.1), pcfun, c(min.opt, max.opt, tol = tol))$par
  out

  N        <- 2000000
  lp       <- stats::rnorm(N, mean = out[1], sd = sqrt(out[2]))
  p        <- (1 + exp(-lp)) ^ (-1)
  y        <- stats::rbinom(N, 1, prob = p)
  prev     <- mean(y)
  #c        <- as.vector(pROC::roc(y, lp, quiet = TRUE)$auc)
  c        <- quickcstat(y, lp)
  c(out[1], out[2], prev, c)
}

# Check
# round(find_mu_sigma(0.1, 0.65, tol=0.00001),4)

########################################################################################

# Find coefficients for desired c-statistic and prevalence (correlated predictors)
# Part of functions calculating expected shrinkage and MAPE

adjust_multiplier_correlated <- function(c, mean, beta, n.predictors, min.opt = variance/3, max.opt = variance*3, tol=0.01, cor0, cor1){

  N <- 3000000

  n.noise <- length(beta[beta==0])
  n.true  <- n.predictors - n.noise

  # mean_var <- find_mu_sigma(p, c, tol=0.001)

  # beta_new <- find_multiplier(sqrt(variance), beta=beta, n.true=n.true, n.noise=n.noise)

  # Specify correlation matrix
  sigma <- matrix(0, nrow = n.predictors,  ncol = n.predictors)
  sigma[1:n.true, 1:n.true] <- cor0
  if (n.noise>0) {
    sigma[(n.true+1):n.predictors, (n.true+1):n.predictors] <- cor1}
  diag(sigma) <- 1

  x     <- mvtnorm::rmvnorm(round(N), rep(0,nrow(sigma)), sigma = sigma )

  cfun <- function(adjust){

    eta   <- mean + x %*% as.matrix(beta) * adjust
    y     <- stats::rbinom(N,  1, invlogit(eta))
    cest  <- quickcstat(y,eta)

    abs(cest- c)
  }
  adjust <- optimize(cfun, c(0, 7, tol = 0.0005))$minimum
  adjust

  beta_new <- beta*adjust
  beta_new

  # Check beta's are ok
  #eta <- mean + x%*%as.matrix(beta_new)
  #y     <- rbinom(N,  1, invlogit(eta))
  #mean(y)
  #quickcstat(y, eta)

}

########################################################################################

# Riey's function (Statistics in Medicine) to approximate R2 from AUC and prevalence)

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



########################################################################################

# Approximation of the C-statistic (Large n)
# Part of various functions above

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

########################################################################################





