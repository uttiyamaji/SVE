library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)
require(fftwtools)
##########################################################
### Fast HAC from paper #####
##########################################################
kweightsHAC <- function ( kernel = c("Truncated", "Bartlett", "Parzen",
                                     "Tukey - Hanning", "Quadratic Spectral"),dimN , bw)
{
  ww <- numeric ( dimN )
  switch (kernel ,
          "Truncated" = {
            ww [1: bw] <- 1
          },
          "Bartlett" = {
            ww [1: bw] <- 1 - ( seq (1, bw) / (bw))
          }, "Parzen" = {
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


HAC.new <- function(mcond , method , bw)
{
  dimmcond <- dim(mcond)
  Nlen <- dimmcond [1]
  qlen <- dimmcond [2]
  ww <- kweightsHAC(kernel = method, Nlen, bw)
  ww <- c(1, ww[1:(Nlen-1)], 0, ww[(Nlen-1):1])
  ww <- Re(fftw(ww))
  FF <- rbind(mcond, matrix(0, Nlen, qlen))
  FF <- mvfftw(FF)
  FF <- FF * matrix(rep(ww, qlen), ncol=qlen)
  FF <- Re(mvfftw(FF, inverse = TRUE)) / (2*Nlen)
  FF <- FF [1: Nlen ,]
  return((t(mcond) %*% FF) / Nlen)
}

##########################################################
### First attempt: calling fft from RcppArmadillo
##########################################################
sourceCpp("msvec_new.cpp")

##########################################################
### Second attempt: smarter coding using RcppArmadillo
##########################################################
sourceCpp("msvec_new_3.cpp")

##########################################################
### Third attempt: further removing matrices
##########################################################
sourceCpp("msvec_new_4.cpp")

N <- 1e3
p <- 3
mcond <- matrix(rnorm(N*p), ncol = p, nrow = N)
mcond <- scale(mcond, center = TRUE, scale = FALSE)

b <- 100

# testing if equal, may take time, so commented out
all.equal(HAC.new(mcond, method = "bartlett", b = b), msveC_new_1(chain = mcond, b = b), msveC_new_2(chain = mcond, b = b))

benchmark(HAC.new(mcond, method = "bartlett", b = b), msveC_new_1(chain = mcond, b = b), msveC_new_2(chain = mcond, b = b), replications = 100)

## See how the estimators scale - for N = 1e6, more than twice the time.
N <- 1e5
p <- 30
mcond <- matrix(rnorm(N*p), ncol = p, nrow = N)
mcond <- scale(mcond, center = TRUE, scale = FALSE)
b <- 100

benchmark(HAC.new(mcond, method = "bartlett", b = b), msveC_new_1(chain = mcond, b = b), msveC_new_2(chain = mcond, b = b), replications = 10)

