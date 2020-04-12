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

initseq_new = function(foo){
  
  for(m in 0:(n/2)){
    sig = HAC1(foo, m)
    
    if(!(any(eigen(sig)$values <= 0))){
      sn = m
      break
    }
    if(sn > (n/2 - 1)){
      stop("Not enough samples")
    }
  }
  
  Dtm = det(sig)
  for(m in (sn+1):(n/2)){
    sig1 = HAC1(foo, m)
    dtm = det(sig1)
    
    if(dtm <= Dtm[m-sn])
      break
    
    Dtm = c(Dtm, dtm)
    sig = sig1
    
  }
  
  return (sig)
  
}