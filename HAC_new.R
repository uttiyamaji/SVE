

# function to calculate SVE 
HAC2 <- function (mcond , method , bw)
{
  require ( fftwtools )
  
  dimmcond <- dim( mcond )
  Nlen <- dimmcond [1]
  qlen <- dimmcond [2]
  
  ww = numeric(Nlen)
  ww <- kweightsHAC ( kernel = method , Nlen , bw)
  
  ww <- c(1, ww [1:( Nlen -1)], 0, ww [( Nlen -1) :1])
  
  ww <- Re( fftw_r2c(ww))            ### change
  
  FF = matrix(0,ncol = qlen, nrow = 2*Nlen)  ### change
  FF[1:Nlen,] <- mcond                      ### change
  
  FF <- mvfftw_r2c (FF)          ### change
  FF <- FF * matrix ( ww , nrow = 2*Nlen, ncol = qlen ) ### change
  FF <-  mvfftw_c2r(FF) / (2* Nlen )  ### change
  FF <- FF [1: Nlen ,]
  return ((t( mcond ) %*% FF) / Nlen )
}










# to get the weights
kweightsHAC <- function ( kernel = c(" Truncated ", " Bartlett ", " Parzen ", "Tukey - Hanning ", " Quadratic Spectral "),
                          dimN , bw){
  ww <- numeric ( dimN )
  switch (kernel ,
          Truncated = {
            ww [1: bw] <- 1
          },
          Bartlett = {
            ww [1: bw] <- 1 - ( seq (1, bw) / (bw))
          }, Parzen = {
            seq1 <- (seq (1, floor (bw/2))) / bw
            seq2 <- (seq( floor (bw/2) +1, bw)) / bw
            ww [1: length ( seq1 )] <- 1 - 6* seq1 ^2 + 6* seq1 ^3
            ww [( length ( seq1 )+1) :bw] <- 2*(1- seq2 )^3
          }, "Tukey - Hanning" = {
            ww [1: bw] <- (1 + cos(pi*(( seq (1, bw))/bw))) / 2
          }, "Quadratic Spectral" = {
            aa <- pi*(( seq (1, dimN ))/bw)/5
            ww <- 1/(12*aa ^2) * (sin (6*aa) / (6*aa) - cos (6*aa))
          })
  return (ww)
}


