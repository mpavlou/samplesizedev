find_n_prap_s <- function(c, p, mean, variance, n.predictors, r2, l_s=0.85, u_s=1.15, PAP_s = 0.8, min.opt = 0.1, max.opt = 0.99){


  if (c>0.75) c_adj     <- adjusted_c_mu_sigma(mean, variance, n.predictors, p) else c_adj = c(c,r2)
  if (c>0.75) c_adj_var <- adjusted_c_mu_sigma(mean, variance, 1, p) else c_adj_var = c(c,r2)


  prob <- function(S){
    n  <- (n.predictors)/ ((S-1)*log(1-  c_adj[2]/S));
    se <- sqrt(S/ (2*p*(1-p)*stats::qnorm(c_adj_var[1])^2 *n) + 2*S^2/(n-2))

    zz <- 1-(stats::pnorm( (l_s-S)/se) +pnorm( (S-u_s)/se) )
    abs(zz - PAP_s)
  }

  s_est <- stats::optimize(prob, c(min.opt, max.opt, tol = 0.0001))$minimum
  n <- (n.predictors)/ ((s_est-1)*log(1 - c_adj[2]/s_est))
  se<- sqrt(s_est/ (2 * p * (1-p) * stats::qnorm(c_adj_var[1])^2 *n) + 2 * s_est^2/(n-2))
  probability <- 1-(stats::pnorm( (l_s-s_est)/se) + stats::pnorm( (s_est-u_s)/se) )

  # if (c>0.8  & c<=0.85 )       {inflation_f   <- 1.2  ; n <- n*inflation_f }
  # if (c>0.85  & c<=0.9 )       {inflation_f   <- 1.3  ; n <- n*inflation_f }

  c(round(n), s_est, se,probability)

}

# find_n_prap_s(c=c, p=0.2, mean=mean_eta, variance=variance_eta, r2=r2, n.predictors=10, l_s=0.85, u_s=1.15, PAP_s = 0.8, min.opt = 0.1, max.opt = 0.99)


