library(mcmc)

HAC1 <- function (mcond , m )
{
  require ( fftwtools )
  dimmcond <- dim( mcond )
  Nlen <- dimmcond [1]
  qlen <- dimmcond [2]
  ww <- weight(m,Nlen)
  ww <- c(1, ww [1:( Nlen -1)], 0, ww [( Nlen -1) :1])
  ww <- Re( fftw (ww))
  FF <- rbind (mcond , matrix (0, Nlen , qlen ))
  FF <- (mvfftw (FF))
  FF <- FF * matrix ( rep(ww , qlen ), ncol = qlen )
  FF <- Re( mvfftw (FF , inverse = TRUE )) / (2* Nlen )
  FF <- FF [1: Nlen ,]
  return ((t( mcond ) %*% FF) / Nlen )
}

weight = function(m, dimN){
  w = numeric(dimN)
  w[1:(2*m+1)] = 1
  return (w)
}


initseq_new = function(foo)
{
  truncs <-  apply(foo, 2, function(t) length(initseq(t)[[2]]) -1)
  print(truncs)
  chosen <- min(truncs)

  sig <- HAC1(foo, chosen)

  # for(m in 0:(n/2)){
  #   sig = HAC1(foo, m)
    
  #   if(!(any(eigen(sig)$values <= 0))){
  #     sn = m
  #     break
  #   }
  #   if(sn > (n/2 - 1)){
  #     stop("Not enough samples")
  #   }
  # }
  
  # Dtm = det(sig)
  # for(m in (sn+1):(n/2)){
  #   sig1 = HAC1(foo, m)
  #   dtm = det(sig1)
    
  #   if(dtm <= Dtm[m-sn])
  #     break
    
  #   Dtm = c(Dtm, dtm)
  #   sig = sig1
    
  # }
  return(list("Sig" = sig, "trunc" = chosen) )
}


library(Rcpp)
library(rbenchmark)
library(mAr)
sourceCpp("inseq.cpp")

p <- 5
n <- 1e4
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
foo <- matrix(rnorm(p^2, sd = 1), nrow = p)
foo <- foo %*% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)

foo <- scale(as.matrix(mAr.sim(rep(0,p), phi, omega, N = n)), scale = F)

a = initseq_new(foo)
b = inseq(foo)
c(a$trunc, b$trunc)
c(det(a$Sig), det(b$Sig))^(1/p)
all.equal(a$Sig,b$Sig)

benchmark(initseq_new(foo),inseq(foo), replications = 2)


library(lineprof)
library(mcmcse)
lineprof(mcse.initseq(foo))
lineprof(initseq_new(foo))
