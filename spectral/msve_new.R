library(fftwtools)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("lag.cpp")

mSVEfft <- function (A, b, method = "bartlett")
{
  n <- nrow(A) # A must be centered matrix
  p <- ncol(A)
  w <- as.vector(lag(1:b, n = n, b = b, method = method)) # calculate lags
  w <- c(1, w[1:(n-1)], 0, w[(n-1):1])  # starting steps from FFT paper
  w <- Re(fftw_r2c(w))            
  FF <- matrix(0, ncol = p, nrow = 2*n)  
  FF[1:n,] <- A    
  if(p > 1)  # multivariate
  {
    FF <- mvfftw_r2c (FF)        
    FF <- FF * matrix(w, nrow = 2*n, ncol = p) 
    FF <- mvfftw_c2r(FF) / (2* n ) 
    return ((t(A) %*% FF[1:n, ]) / n )    
  } else if(p == 1)  ##univariate calls
  {
    FF <- fftw_r2c (FF)        
    FF <- FF * matrix(w, nrow = 2*n, ncol = p) 
    FF <- fftw_c2r(FF) / (2* n ) 
    return ((t(A) %*% FF[1:n]) / n )
  }              

}


