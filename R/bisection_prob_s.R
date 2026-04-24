# Function for bisection method

bisection_prob_s <- function(f, a, b, iter = 20, tol = ceiling(round(a/200)/5) * 5, nsim = 1000) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  ##ch

  tol_n = round(a/50)

  c<-round((a+b)/2)

  fc <- f(c, nsim = nsim)
  if (abs(fc) <= 0.005) return(c)

  fa <- f(a, nsim = nsim)
  if (abs(fa) <= 0.0051) return(a)

  fb <- f(b, nsim = nsim)
  if (abs(fb) <= 0.0051) return(b)


  if (fa > 0 & fb > 0) {
    a<-a*0.95  ;  fa <- f(a, nsim = nsim) ; {if (abs(fa) <= 0.0051) return(a)}}
  if (fa < 0 & fb < 0) {
    b<-b*1.05  ;  fb <- f(b, nsim = nsim);  {if (abs(fb) <= 0.0051) return(b)}}


  if (!(fa < 0) && (fb > 0)) {
    stop('signs of f(a) and f(b) do not differ')
  } else if ((fa > 0) && (fb < 0)) {
    stop('signs of f(a) and f(b) do not differ')
  }


  for (k in 1:iter) {

    c  <- round((a*3 + b*2) / 5) # Calculate midpoint
    fc <- f(c, nsim = nsim)


    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.

    if (nsim<=200){
    if ( ((abs(fc) <= 0.1) || ((b - a) / 2) < tol_n)   &  (k >=1 )) {
      #if (  abs(fc) <= 0.0025 &  (k >=2 )) {
      return(c)
    }
    }

    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the
    # function and return the root.
   if (  ((abs(fc) <= 0.0051) || ((b - a) / 2) < tol_n)   &  (k >=1 )) {
    #if (  abs(fc) <= 0.0025 &  (k >=2 )) {
      return(c)
   }

    if (  ((abs(fc) <= 0.006) || ((b - a) / 2) < tol_n)   &  (k >=7) ) {
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
