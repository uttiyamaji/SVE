require(fftwtools)

# function to calculate SVE 

SVE = function (A, b, method = "Bartlett"){
  n = nrow(A)
  p = ncol(A)
  
  w = lag_weights(kernel = method, n, b)
  
  w = c(1, w[1:(n-1)], 0, w[(n-1):1])
  
  w = Re(fftw_r2c(w))            
  
  FF = matrix(0, ncol = p, nrow = 2*n)  
  FF[1:n,] = A                      
  
  FF = mvfftw_r2c (FF)        
  FF = FF * matrix(w, nrow = 2*n, ncol = p ) 
  FF = mvfftw_c2r(FF) / (2* n ) 

  return ((t(A) %*% FF[1:n, ]) / n )
}










# to get the weights
lag_weights = function(kernel = c("Truncated", "Bartlett", "Parzen", "Tukey - Hanning", "Quadratic Spectral"),
                   dimN, b){
  w = numeric(dimN)
  switch (kernel,
          Truncated = {
            w [1:b] = 1
          },
          Bartlett = {
            w [1:b] = 1 - (seq(1, b)/b)
          }, 
          Parzen = {
            seq1 = (seq(1, floor(b/2))) / b
            seq2 = (seq(floor(b/2) + 1, b)) / b
            w[1: length(seq1)] = 1 - 6*seq1^2 + 6*seq1 ^3
            w[(length(seq1)+1) : b] = 2*(1- seq2)^3
          }, 
          "Tukey - Hanning" = {
            w [1:b] = (1 + cos(pi*((seq(1,b)) / b))) / 2
          }, 
          "Quadratic Spectral" = {
            aa = pi*((seq (1, dimN))/b)/5
            w = 1/(12*aa ^2) * (sin(6*aa) / (6*aa) - cos(6*aa))
          })
  return (w)
}


