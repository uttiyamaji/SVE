# p = 30
# n = 1e5

# foo <- scale(matrix(rnorm(n*p), nrow = n, ncol = p), scale = F)

# sourceCpp("roncalli_vs_usual.cpp")
# rbenchmark::benchmark(roncalli(foo,20), usual(foo, 20), replications = 10)

library(mAr)
library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)
sourceCpp("inseq.cpp")
sourceCpp("inseq_roncalli.cpp")

p <- 30
n <- 1e5

omega <- 5*diag(1,p)

## Making correlation matrix var(1) model
set.seed(100)
foo <- matrix(rnorm(p^2, sd = 4), nrow = p)
foo <- foo %*% t(foo)
phi <- foo / (max(eigen(foo)$values) + 1)

out <- as.matrix(mAr.sim(rep(0,p), phi, omega, N = n))
out <- scale(out, scale = FALSE)

all.equal(inseq_test(out), inseq(out))
benchmark(inseq_test(out), inseq(out), replications = 2)


