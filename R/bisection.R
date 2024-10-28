# Function for bisection method

bisection <- function(f, a, b, iter = 10, tol = ceiling(round(a/200)/5) * 5, nsim = 1000) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  ##ch

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

    c  <- round((a*3 + b*2) / 5) # Calculate midpoint
    fc <- f(c, nsim = nsim)

    #print(c(k, a, b, c, fa, fb, fc ))
    #print(c(k))



    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.

    if (nsim<=200){
    if ( ((abs(fc) <= 0.007) || ((b - a) / 2) < tol)   &  (k >=1 )) {
      #if (  abs(fc) <= 0.0025 &  (k >=2 )) {
      return(c)
    }
    }

    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
   if (  ((abs(fc) <= 0.0025) || ((b - a) / 2) < tol)   &  (k >=1 )) {
    #if (  abs(fc) <= 0.0025 &  (k >=2 )) {
      return(c)
   }

    if (  ((abs(fc) <= 0.005) || ((b - a) / 2) < tol)   &  (k >=7) ) {
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
