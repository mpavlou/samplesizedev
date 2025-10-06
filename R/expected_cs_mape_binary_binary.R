# Adjust coefficients for correlated predictors
adjust_multiplier_binary <- function(c, p, beta = c(-0.5,-0.3,0.3,0.15,0.15), type="binary", n.true, n.noise, N=300000){

  # beta = c(-0.5,-0.3,0.3,0.15,0.15)
  # N <- 3000000

  beta <- c(beta,rep(0.1,n.true-5), rep(0,n.noise))
  x    <- rmvbin(round(N), margprob=c(0.7, 0.7, 0.3, 0.4, 0.5, rep(0.5,n.true-5), rep(0.2,n.noise)))

  cfun <- function(b0_adjust){

    b0     <- b0_adjust[1]
    adjust <- b0_adjust[2]

    eta   <- b0 + x %*% as.matrix(beta) * adjust
    y     <- rbinom(N,  1, invlogit(eta))
    cest  <- quickcstat(y,eta)
    prev  <- mean(y)

    abs(cest - c) + abs(prev-p)
  }


  run <- round(optim(par=c(-2.65,1), cfun, c(min.opt=c(-7,0.3),  max.opt=c(6,13), tol = 0.0005))$par, 4)

  b0      <-  run[1]
  adjust  <-  run[2]

  #check
  # x    <- rmvbin(round(N), margprob=c(0.7, 0.7, 0.3, 0.4, 0.5, rep(0.5,n.true-5), rep(0.2,n.noise)))

  # eta <- b0 + x %*% as.matrix(beta) * adjust
  # y     <- rbinom(N,  1, invlogit(eta))
  # mean(y)
  # quickcstat(y, eta)

  c(b0, beta*adjust)

}

#mean_var


### Expected shrinkage Binary predictors

expected_s_binary <- function(n, S, p, c, beta = c(-0.5,-0.3,0.3,0.15,0.15) , n.true=5, n.noise=5, r2=0, nsim=50, nval=50000){

  #  beta = c(-0.5,-0.3,0.3,0.15,0.15)

  n.predictors <- n.true + n.noise

  set.seed(2022)

  beta <- adjust_multiplier_binary(c, p, beta = beta, n.true = n.true, n.noise = n.noise)

  cs           <- NULL
  mape         <- NULL

  #xval    <- rmvbin(1000000, margprob=c(0.7, 0.7, 0.3, 0.4, 0.5, rep(0.2,n.noise)))
  #yval    <- rbinom(1000000, 1, invlogit(cbind(1,xval)%*%beta))
  #p_true  <- invlogit(cbind(1,xval)%*%beta)
  # p_est <- mean(yval)
  # c_est <-quickcstat(yval,p_true)
  # p_est
  # c_est

  xval    <-rmvbin(round(nval), margprob=c(0.7, 0.7, 0.3, 0.4, 0.5, rep(0.5,n.true-5), rep(0.2,n.noise)))
  yval    <- rbinom(nval, 1, invlogit(cbind(1,xval) %*% beta))


  for (i in 1: nsim){

    set.seed(i)


    x     <- rmvbin(n, margprob=c(0.7, 0.7, 0.3, 0.4, 0.5, rep(0.5,n.true-5), rep(0.2,n.noise)))
    y     <- rbinom(n,  1, invlogit(cbind(1,x)%*%beta))

    yval    <- rbinom(nval, 1, invlogit(cbind(1,xval) %*% beta))
    p_true <- as.vector(invlogit(cbind(1,xval) %*% beta))

    #a       <- fastglm(cbind(1,x), y, family=binomial())
    a        <- fastLR(cbind(1,x), y)
    eta_est  <- cbind(1, xval) %*% as.vector(a$coef)
    p_est    <- as.vector(invlogit(eta_est))


    #fit      <- fastglm(cbind(1,eta_est), yval, family=binomial())
    fit      <- fastLR(cbind(1,eta_est), yval )
    cs[i]    <- fit$coef[2]
    mape[i]  <- mean(abs(p_true-p_est))

  }
  #hist(cs, main = paste("CS=",round(median(cs,na.rm=TRUE)/0.0025)*0.0025, "N=",n))
  df <- data.frame(cs)
  pl <- ggplot(df, aes(x=cs)) +
    geom_density() + ggtitle(paste("Target C = ", c, ",  Target Prevalence = ", p,"\n", "N = ",n, ",  Expected CS = ",round(median(cs,na.rm=TRUE)/0.005)*0.005,
                                   ", Expected MAPE=", round(median(mape,na.rm=TRUE),4), sep="")) +
    geom_vline(aes(xintercept=median(cs)), color="blue", linetype="dashed", size=1)
  print(pl)
  df<- c(round(median(cs,na.rm=TRUE)/0.005)*0.005, sqrt(var(cs,na.rm=TRUE)), round(median(mape,na.rm=TRUE),4) )


  xval    <-rmvbin(2000000, margprob=c(0.7, 0.7, 0.3, 0.4, 0.5, rep(0.5, n.true-5), rep(0.2, n.noise)))
  yval    <- rbinom(2000000, 1, invlogit(cbind(1,xval)%*%beta))

  set.seed(2022)

  prev      <- mean(yval)
  cstat     <- quickcstat(yval, invlogit(cbind(1,xval) %*% beta))

  df        <- data.frame(n, ceiling(mean(cs, na.rm = TRUE)/0.0025) * 0.0025,
                          round(sqrt(stats::var(cs,na.rm = TRUE)), 4),
                          round(sqrt( mean( ((cs-1)^2), na.rm=TRUE) ), 4),
                          round(mean(ifelse( (cs < 0.8), 1, 0),na.rm=TRUE), 2),
                          round(mean(mape, na.rm = TRUE),4),
                          round(sqrt(stats::var(mape,na.rm = TRUE)), 4),
                          round(prev, 3),
                          round(cstat, 3 ),
                          n.predictors)
  names(df) <- c("N", "Mean_CS", "SD_CS", "RMSD_CS", "Pr(CS<0.8)", "Mean_MAPE",  "SD_MAPE", "Prev.", "C-Stat.", " # Predictors")

  performance <- df[,-3]
  performance <- df

  performance

}
