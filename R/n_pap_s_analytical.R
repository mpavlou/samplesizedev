n_pap_s_analytical <- function(c, p, n.predictors, l_s, u_s, PAP_s = PAP_s, min.opt = 0.1, max.opt = 0.99){

  r2   <- as.numeric(approximate_R2(c, p, n = 1000000)[2])

  prob <- function(S){

    n <- (n.predictors)/ ((S-1)*log(1-r2/S));
    se<- sqrt(S/ (2*p*(1-p)*qnorm(c)^2 *n) + 2*S^2/(n-2))

    zz <- 1-(stats::pnorm( (l_s - S)/se) +pnorm( (S - u_s )/se) )
    abs(zz - PAP_s)

  }

  s_est       <- stats::optimize(prob, c(min.opt, max.opt, tol = 0.001))$minimum
  n           <- (n.predictors)/ ((s_est-1)*log(1-r2/s_est))

  se          <- sqrt(s_est/ (2*p*(1-p)*qnorm(c)^2 *n) + 2*s_est^2/(n-2))

  probability <- 1-(stats::pnorm( ( l_s - s_est)/se) +pnorm( (s_est - u_s )/se) )

  c(s_est, round(n),se,probability)

}
