
# Function for bisection method

bisection <- function(f, a, b, tol, iter = 10, tol_n = ceiling(round(a/100)/5) * 5, nsim = 1000) {

  tol_n = round(a/50)

  fa <- f(a, nsim = nsim)
  if (abs(fa) <= tol) return(a)

  fb <- f(b, nsim = nsim)
  if (abs(fb) <= tol) return(b)


  if (fa > 0 & fb > 0) { a<-a*0.8;   fa <- f(a, nsim = nsim)}
  if (fa < 0 & fb < 0) { b<-b*1.2;   fb <- f(b, nsim = nsim)}


  if (!(fa < 0) && (fb > 0)) {
    stop('signs of f(a) and f(b) do not differ')
  } else if ((fa > 0) && (fb < 0)) {
    stop('signs of f(a) and f(b) do not differ')
  }


  for (k in 1:iter) {

    # if (k <= divide ) nsim <- nsim1/divide*k
    # a <- round(a/tol)*tol
    # b <- round(b/tol)*tol

    c  <- round((a*3 + b*2) / 5) # Calculate midpoint
    fc <- f(c, nsim = nsim)

    #print(c(k, a, b, c, fa, fb, fc ))
    #print(c(k))



    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.

    if (nsim<=200){
    if ( ((abs(fc) <= 0.01) || ((b - a) / 2) < tol_n)   &  (k >=1 )) {
      #if (  abs(fc) <= 0.0025 &  (k >=2 )) {
      return(c)
    }
    }

    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
   if (  ((abs(fc) <= tol) || ((b - a) / 2) < tol_n)   &  (k >=1 )) {
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
