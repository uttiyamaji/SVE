library(Rcpp)
library(rbenchmark)
library(mAr)
sourceCpp("inseq_fft_c++.cpp")
sourceCpp("inseq.cpp")

p = 30
n = 1e5

foo <- scale(matrix(rnorm(n*p), nrow = n, ncol = p), scale = F)

a = initseq_new(foo)
b = inseq(foo)

all.equal(a$Sig,b$Sig)

benchmark(initseq_new(foo),inseq(foo), replications = 10)


### Correlated sampling problem ####
p <- 3
n <- 1e3
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
foo <- matrix(rnorm(p^2), nrow = p)
foo <- foo %*% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)

foo <- scale(as.matrix(mAr.sim(rep(0,p), phi, omega, N = n)), scale = F)
a <- initseq_new(foo)
b <- inseq(foo)
all.equal(a$Sig,b$Sig)

benchmark(initseq_new(foo),inseq(foo), replications = 10)



## another setting
### Correlated sampling problem ####
p <- 3
n <- 1e5
omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
foo <- matrix(rnorm(p^2), nrow = p)
foo <- foo %*% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)

foo <- scale(as.matrix(mAr.sim(rep(0,p), phi, omega, N = n)), scale = F)
a <- initseq_new(foo)
b <- inseq(foo)
all.equal(a$Sig,b$Sig)

benchmark(initseq_new(foo),inseq(foo), replications = 10)

