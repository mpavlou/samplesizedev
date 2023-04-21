# Function for bisection method

bisection <- function(f, a, b, iter = 1000, tol = 1e-7, nsim=200) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  
  nsim1 <- nsim
  nsim  <- nsim/4
  
  fa <- f(a, nsim = nsim)
  fb <- f(b, nsim = nsim)
  if (!(fa < 0) && (fb > 0)) {
    stop('signs of f(a) and f(b) differ')
  } else if ((fa > 0) && (fb < 0)) {
    stop('signs of f(a) and f(b) differ')
  }
  
  for (k in 1:iter) {
    
    if (k<=4) nsim <- nsim1/4*k
    print(c(k))
    a <- round(a/10)*10
    b <- round(b/10)*10
    
    c  <- (a + b) / 2 # Calculate midpoint
    fc <- f(c, nsim=nsim)
    
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if (  ((fc == 0) || ((b - a) / 2) < tol)   & (k>=4)) {
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