##################################################################
# Find mean and variance of linear predictor for given prevalence
# and C-statistic

invlogit <- function(x) 1/(1+exp(-x))

logit <-function(x) log(x/(1-x))


find_mu_sigma <- function(target.prev, target.c, min.opt = c(-10,0), max.opt = c(0.02,5), tol = 0.00001){


    pcfun <- function(x){

    #target.prev = 0.2; target.c=0.7; min.opt = c(-7,0.5); max.opt = c(0,14)

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

    f3     <- function(x) stats::dnorm(x, mean=mean, sd = sqrt(variance)) * (1  + exp(-x)) ^ (-1)

    c      <- num/denom
    prev   <- stats::integrate(f3, - Inf, Inf, subdivisions = 1000L)$value

    abs( c - target.c)^2 + abs(prev - target.prev )^2
  }

  if (target.c>0.65) {
    sigma_c <- sqrt(2) * stats::qnorm(target.c)
    mu      <- 0.5 * (2 * target.prev - 1) * (sigma_c^2) + log(target.prev / (1 - target.prev))
  out      <- stats::optim(par=c(mu,0.15), pcfun, c(min.opt, max.opt, tol = tol))$par} else

    {
    sigma_c <- sqrt(2) * stats::qnorm(target.c)
    mu      <- 0.5 * (2 * target.prev - 1) * (sigma_c^2) + log(target.prev / (1 - target.prev))
    sigma   <- sqrt((sigma_c^2) * (1 + target.prev * (1 - target.prev) * (sigma_c^2)))
    out     <- c(mu, sigma^2)
    }

  N        <- 500000
  #Better 2000000 to check
  lp       <- stats::rnorm(N, mean = out[1], sd = sqrt(out[2]))
  p        <- (1 + exp(-lp)) ^ (-1)
  y        <- stats::rbinom(N, 1, prob = p)
  prev     <- mean(y)
  c        <- quickcstat(y, lp)
  c(out[1], out[2], prev, c)
}

# Check
    # round(find_mu_sigma(0.02, 0.675, tol=0.00001),4)
  # round(find_mu_sigma(0.1, 0.6, tol=0.00001),4)

 # round(find_mu_sigma(0.05, 0.7, tol=0.00001),4)




#####################################################


#########################################
# Find scale and variance of linear predictor for given prevalence and C-statistic

#
# find_lambda_sigma <- function(target.prev, target.c, time=5, min.opt = c(0.001,0.02), max.opt = c(0.5,0.5),tol=0.0001){
#
#   pcfunsurv <- function(x){
#
#     #target.prev = 0.2; target.c=0.7; min.opt = c(-7,0.5); max.opt = c(0,14)
#
#     lambda    <- x[1]
#     variance  <- x[2]
#
#     f3     <- function(y) stats::dnorm(y, mean = 0, sd = sqrt(variance)) * (1 - (exp(-time * lambda)) ^ (exp(y)) )
#    #  f3     <- function(y) stats::dnorm(y, mean = 0, sd = sqrt(variance)) * (1  + exp(-y)) ^ (-1)
#
#     prev   <- stats::integrate(f3, - Inf, Inf, subdivisions = 1000L)$value ;prev
#
#     f1 = function(y) {
#       stats::integrate(function(z) {stats::dnorm(y, mean = 0, sd = sqrt(variance)) * stats::dnorm(z, mean = 0, sd = sqrt(variance)) * (1 - (exp(-time * lambda)) ^ exp(y) ) * (exp(-(time-1) * lambda) ^ exp(z)) },
#                        -Inf, y)$value
#     }
#
#     num = stats::integrate(Vectorize(f1), -Inf, Inf)$value; num
#
#     f2 = function(y) {
#       stats::integrate(function(z) {stats::dnorm(y, mean = 0, sd = sqrt(variance)) * stats::dnorm(z, mean = 0, sd = sqrt(variance)) * (1 - (exp(-time * lambda)) ^ exp(y) ) * (exp(-(time-1) * lambda) ^ exp(z)) },
#                        -Inf, Inf)$value
#     }
#
#     denom <- stats::integrate(Vectorize(f2), -Inf, Inf)$value ;denom
#
#     c      <- num/denom ;c
#
#      # c <- NULL
#      # t <- c(1, 2, 5, 10, 20,30, 40, 100, 150, 200, 1000)
#      # for (time in t) {
#      #   c <- c(c, stats::integrate(Vectorize(f1), -Inf, Inf)$value / stats::integrate(Vectorize(f2), -Inf, Inf)$value)
#      # }
#      #
#      # mean(c)
#      #
#      # c <- weighted.mean(c, w = pweibull(t, shape=1, scale=1/lambda) )
#      # c
#
#     abs( c - target.c)^2 + abs(prev - target.prev )^2
#   }
#
#   out      <- stats::optim(par=c(0.01, 0.5), pcfunsurv, c(min.opt, max.opt, tol = tol))$par
#   out
#
#   N        <- 1000000
#   lp       <- stats::rnorm(N, mean = 0, sd = sqrt(out[2]))
#   p        <- (1 - (exp(-time * out[1])) ^ exp(lp) )
#   # y        <- stats::rbinom(N, 1, prob = p)
#   # prev     <- mean(y)
#   # c        <- as.vector(pROC::roc(y, lp, quiet = TRUE)$auc)
#   # c         <-quickcstat(y, lp)
#   c(out[1], out[2], mean(p))
# }
#
# # Check
#
# out <- round(find_lambda_sigma (0.2, 0.7, time=2, tol=0.00001),4)
# out
#
# n      <- 10000000
# thresh <- 0.001
#
# u      <- runif(n)
# lp     <- stats::rnorm(n, mean = 0, sd = sqrt(out[2]))
# t      <- -log(u)/((out[1]) * exp(lp) )
# censor <- rep(1,n)
# censor[t > quantile(t, p = thresh)]=0
# mean(censor)
# quickcindex(censor, t, lp, seed=1)
#
# n      <- 10000
# thresh <- 0.6
#
# u      <- runif(n)
# lp     <- stats::rnorm(n, mean = 0, sd = sqrt(out[2]))
# t      <- -log(u)/((out[1]) * exp(lp) )
# censor <- rep(1,n)
# censor[t > quantile(t, p = thresh)]=0
# mean(censor)
#
# data<-data.frame(t, censor, lp)
#
# a <- cindex(Surv(t, censor) ~ lp, data=data)
# a
#
#
#
#
#
#
#
# thresh <- 0
#
# u      <- runif(n)
# lp     <- stats::rnorm(n, mean = 0, sd = sqrt(out[2]))
# t      <- -log(u)/((out[1]) * exp(lp) )
# censor <- rep(1,n)
# censor[t > quantile(t, p = thresh)]=0
# mean(censor)
#
#
# quickcindex_nocensoring(t,lp)
#
#

find_sigma <- function(target.c,  min.opt = 0.1, max.opt = 15, tol=0.0001){

  n       <- 50000
  lambda  <- 1

  # Build formula
  # measurevar <- "Surv(t, censor)"
  etavar     <- "eta"
  formula    <- stats::as.formula  (paste("Surv(t, censor)", etavar, sep=" ~ "))


  cfun <- function(x) {
    set.seed(2)
    u      <- stats::runif(n)
    eta    <- stats::rnorm(n, mean = 0, sd = sqrt(x))
    t      <- -log(u)/((lambda) * exp(eta) )
    #c      <- 1 - concordance(t ~ eta)$concordance; c
    censor <- rep(1, n)
    data   <- data.frame(t, censor, eta)

    big <- survival::coxph(survival::Surv(t, censor)~eta, data = data, x = TRUE, y = TRUE)

    A1  <- pec::cindex(list("Cox X1"= big), formula = formula, data = data, eval.times = 1)

    c <- A1$AppCindex[["Cox X1"]]

    abs(c - target.c)
  }

  out      <- stats::optimize(cfun, lower=min.opt, upper=max.opt, tol = tol)$minimum
  out
}

# system.time(x <- find_sigma(target.c=0.65, min.opt=0.2, max.opt=15))
#  x
#  n <- 100000
#  u      <- runif(n)
#  eta    <- stats::rnorm(n, mean = 0, sd = sqrt(x))
#  t      <- -log(u)/((lambda) * exp(eta) ) ;mean(t)
#  gamma  <- 1
#  t      <-  (-log(u)/((lambda*1) * exp(eta) ))^(1/gamma) ;mean(t)
#
#  censor <- rep(1,n)
#  c      <- 1 - concordance(t ~ eta)$concordance; c
#  # Build formula
#  measurevar <- "Surv(t, censor)"
#  etavar     <- "eta"
#  formula    <- as.formula  (paste( measurevar, etavar, sep=" ~ "))
#
#  data <- data.frame(t, censor, eta)
#
# #
#  big <- coxph(formula, data=data, x=TRUE, y=TRUE)
#
#  A1  <- pec::cindex(list("Cox X1"= big), formula=Surv(t, censor)~eta, data=data, eval.times=1)
#
#  cest <- A1$AppCindex[["Cox X1"]]
#  cest

 #####################################################

 find_sigma_quick <- function(target.c,  min.opt = 0.1, max.opt = 15, tol=0.0001){

   n       <- 50000
   lambda  <- 1

   # Build formula
   measurevar <- "Surv(t, censor)"
   etavar     <- "eta"
   formula    <- as.formula  (paste( measurevar, etavar, sep=" ~ "))


   cfun <- function(x) {
     set.seed(2)
     u      <- runif(n)
     eta    <- stats::rnorm(n, mean = 0, sd = sqrt(x))
     t      <- -log(u)/((lambda) * exp(eta) )
     c      <- 1 - survival::concordance(t ~ eta)$concordance; c
     abs(c - target.c)
   }

   out      <- stats::optimize(cfun, lower=min.opt, upper=max.opt, tol = tol)$minimum
   out
 }


# system.time(x <- find_sigma_quick(target.c=0.65, min.opt=0.1, max.opt=15))
#  x
#  n <-100000
#  u      <- runif(n)
#  eta    <- stats::rnorm(n, mean = 0, sd = sqrt(x))
#  t      <- -log(u)/((lambda) * exp(eta) )
#  censor <- rep(1,n)
#  c      <- 1 - concordance(t ~ eta)$concordance; c
#  # Build formula
#  measurevar <- "Surv(t, censor)"
#  etavar     <- "eta"
#  formula    <- as.formula  (paste( measurevar, etavar, sep=" ~ "))
#
#  data <- data.frame(t, censor, eta)
#
# #
#  big <- coxph(formula, data=data, x=TRUE, y=TRUE)
#
#  A1  <- pec::cindex(list("Cox X1"= big), formula=Surv(t, censor)~eta, data=data, eval.times=1)
#
#  cest <- A1$AppCindex[["Cox X1"]]
#  cest





find_mu_sigma_sim <- function(target.prev, target.c, min.opt = c(-10,0), max.opt = c(4,5), tol = 0.00001){

  pcfun <- function(x){

    #target.prev = 0.2; target.c=0.7; min.opt = c(-7,0.5); max.opt = c(0,14)

    mean      <- x[1]
    variance  <- x[2]

    N        <- 200000
    lp       <- stats::rnorm(N, mean = mean, sd = sqrt(variance))
    p        <- (1 + exp(-lp)) ^ (-1)
    y        <- stats::rbinom(N, 1, prob = p)
    prev     <- mean(y)
    c        <- quickcstat(y, lp)
    abs(c - target.c) + abs(p-target.prev)
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
   # round(find_mu_sigma(0.17, 0.63, tol=0.00001),4)
   # round(find_mu_sigma_sim(0.17, 0.63, tol=0.001),4)




Hist <- function (time, event, entry = NULL, id = NULL, cens.code = "0",
          addInitialState = FALSE)
{
  cens.code <- as.character(cens.code[[1]])
  if (is.matrix(time))
    time <- data.frame(time)
  if (inherits(time, "list")) {
    if (length(time) != 2 || length(time[[1]]) != length(time[[2]]))
      stop("Argument time has a wrong format")
    time <- data.frame(time)
  }
  if (is.data.frame(time)) {
    cens.type <- "intervalCensored"
    L <- time[[1]]
    R <- time[[2]]
    N <- length(L)
    stopifnot(is.numeric(L))
    stopifnot(is.numeric(R))
    wrong <- L > R
    wrong[is.na(R)] <- FALSE
    stopifnot(all(wrong == FALSE))
    rm(wrong)
    status <- rep(2, N)
    status[L == R] <- 1
    status[is.infinite(R) | is.na(R) | (L != R & as.character(R) ==
                                          cens.code)] <- 0
    R[status == 0] <- Inf
  }
  else {
    stopifnot(is.numeric(time))
    cens.type <- "rightCensored"
    N <- length(time)
    status <- rep(1, N)
  }
  if (is.null(entry))
    entry.type <- ""
  else {
    if (is.matrix(entry))
      entry <- data.frame(entry)
    if (inherits(entry, "list")) {
      if (length(entry) != 2 || length(entry[[1]]) != length(entry[[2]]))
        stop("Argument entry has a wrong format")
      entry <- data.frame(entry)
    }
    if (is.data.frame(entry)) {
      entry.type <- "intervalCensored"
      U <- entry[[1]]
      V <- entry[[2]]
      stopifnot(is.numeric(U))
      stopifnot(is.numeric(V))
      stopifnot(all(!is.na(U)) | all(!is.na(V)))
    }
    else {
      stopifnot(is.numeric(entry))
      if (is.null(id))
        entry.type <- "leftTruncated"
      else entry.type <- "exact"
    }
  }
  if (cens.type == "intervalCensored") {
    if (entry.type == "intervalCensored") {
      stopifnot(all(V <= L))
    }
    else {
      stopifnot(is.null(entry) || all(entry <= L))
    }
  }
  else {
    if (entry.type == "intervalCensored") {
      stopifnot(all(V <= time))
    }
    else {
      stopifnot(is.null(entry) || all(entry <= time))
    }
  }
  if (missing(event)) {
    model <- "onejump"
    event <- rep(1, N)
    warning("Argument event is missing:\nassume observations of a survival model\nand only one event per subject")
  }
  else {
    if (is.matrix(event))
      event <- data.frame(event)
    if ((is.vector(event) & !(inherits(event, "list"))) ||
        is.factor(event))
      stopifnot(length(event) == N)
    if (inherits(event, "list")) {
      if (length(event) != 2 || length(event[[1]]) != length(event[[2]]))
        stop("Argument event has a wrong format")
      event <- data.frame(event)
    }
    if (!is.data.frame(event)) {
      if (is.null(id)) {
        model <- "onejump"
        if (is.logical(event))
          event <- as.numeric(event)
        status[is.na(event) | is.infinite(event) | as.character(event) ==
                 cens.code] <- 0
      }
      else {
        stopifnot(is.numeric(id) || is.factor(id))
        model <- "multi.states"
        if (cens.type == "intervalCensored") {
          stop("Dont know the order of transitions for interval censored observations.")
        }
        else {
          if (addInitialState == TRUE) {
            time <- c(rep(0, length(unique(id))), time)
            if (is.factor(event)) {
              event <- factor(c(rep("initial",
                                    length(unique(id))), as.character(event)),
                              levels = c("initial", levels(event)))
            }
            else {
              stopifnot(match("initial", unique(event),
                              nomatch = 0) == 0)
              event <- c(rep("initial", length(unique(id))),
                         event)
            }
            id <- c(unique(id), id)
          }
          sorted <- order(id, time)
          time <- time[sorted]
          id <- id[sorted]
          event <- event[sorted]
          if (length(unique(id)) != sum(time == 0))
            stop("There are ", length(unique(id)),
                 " different individuals (id's), but the state at time 0 is available for ",
                 sum(time == 0), " id's.")
          initialState <- event[time == 0]
          last.id <- c(diff(id) != 0, 1)
          first.id <- c(1, diff(id) != 0)
          from <- factor(event[last.id != 1])
          to <- factor(event[first.id != 1])
          id <- id[time != 0]
          time <- time[time != 0]
          status <- rep(1, length(to))
          status[is.na(to) | is.infinite(to) | as.character(to) ==
                   cens.code] <- 0
        }
      }
    }
    else {
      model <- "multi.states"
      from <- event[[1]]
      to <- event[[2]]
      status[is.na(to) | is.infinite(to) | as.character(to) ==
               cens.code] <- 0
      if (length(unique(from)) == 1) {
        model <- "onejump"
        event <- to
        if (is.logical(to))
          to <- as.numeric(to)
        status[is.na(to) | is.infinite(to) | as.character(event) ==
                 cens.code] <- 0
      }
    }
  }
  if (all(status == 1))
    cens.type <- "uncensored"
  if (model == "onejump") {
    if (is.factor(event)) {
      event <- factor(event)
      states <- levels(event)
    }
    else {
      states <- sort(as.character(unique(event)))
    }
    states <- as.character(states[states != cens.code])
    if (length(states) > 1)
      model <- "competing.risks"
    else model <- "survival"
    if (cens.type == "intervalCensored") {
      if (model == "survival") {
        if (entry.type == "intervalCensored")
          history <- cbind(U = U, V = V, L = L, R = R,
                           status = status)
        else history <- cbind(entry = entry, L = L, R = R,
                              status = status)
      }
      else {
        if (entry.type == "intervalCensored")
          history <- cbind(U = U, V = V, L = L, R = R,
                           status = status, event = as.integer(factor(event,
                                                                      levels = c(states, cens.code))))
        else history <- cbind(entry = entry, L = L, R = R,
                              status = status, event = as.integer(factor(event,
                                                                         levels = c(states, cens.code))))
      }
    }
    else {
      if (model == "survival") {
        if (entry.type == "intervalCensored")
          history <- cbind(U = U, V = V, time = time,
                           status = status)
        else history <- cbind(entry = entry, time = time,
                              status = status)
      }
      else {
        if (entry.type == "intervalCensored")
          history <- cbind(U = U, V = V, time = time,
                           status = status, event = as.integer(factor(event,
                                                                      levels = c(states, cens.code))))
        else {
          history <- cbind(entry = entry, time = time,
                           status = status, event = as.integer(factor(event,
                                                                      levels = c(states, cens.code))))
        }
      }
    }
  }
  else {
    if (any(as.character(from) == as.character(to)))
      stop("Data contain transitions from state x to state x")
    eventISfactor <- as.numeric(is.factor(from)) + as.numeric(is.factor(to))
    if (eventISfactor == 1)
      stop("Components of event have different classes")
    if (eventISfactor == 2)
      states <- unique(c(levels(from), levels(to)))
    else states <- as.character(unique(c(from, to)))
    states <- as.character(states[states != cens.code])
    if (cens.code %in% levels(from)) {
      stop(paste("The Cens.code", cens.code, " identifies censored data, but is found amoung the `from' state of some transitions"))
    }
    if (cens.type == "intervalCensored") {
      if (entry.type == "intervalCensored")
        history <- cbind(U = U, V = V, L = L, R = R,
                         status = status, from = as.integer(factor(from,
                                                                   levels = c(states, cens.code))), to = as.integer(factor(to,
                                                                                                                           levels = c(states, cens.code))))
      else {
        history <- cbind(entry = entry, L = L, R = R,
                         status = status, from = as.integer(factor(from,
                                                                   levels = c(states, cens.code))), to = as.integer(factor(to,
                                                                                                                           levels = c(states, cens.code))))
      }
    }
    else {
      if (entry.type == "intervalCensored")
        history <- cbind(U = U, V = V, time = time, status = status,
                         from = as.integer(factor(from, levels = c(states,
                                                                   cens.code))), to = as.integer(factor(to,
                                                                                                        levels = c(states, cens.code))))
      else {
        history <- cbind(entry = entry, time = time,
                         status = status, from = as.integer(factor(from,
                                                                   levels = c(states, cens.code))), to = as.integer(factor(to,
                                                                                                                           levels = c(states, cens.code))))
      }
    }
  }
  if (!is.null(id))
    history <- cbind(history, id)
  rownames(history) <- NULL
  class(history) <- c("Hist")
  attr(history, "states") <- states
  attr(history, "cens.type") <- cens.type
  attr(history, "cens.code") <- as.character(cens.code)
  attr(history, "model") <- model
  attr(history, "entry.type") <- entry.type
  history
}



