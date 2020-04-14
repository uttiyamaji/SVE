p = 30
n = 1e5

foo <- scale(matrix(rnorm(n*p), nrow = n, ncol = p), scale = F)


sourceCpp("inseq_fft_c++.cpp")
sourceCpp("inseq.cpp")

a = initseq_new(foo)
b = inseq(foo)

all.equal(a$Sig,b$Sig)

library(rbenchmark)
benchmark(initseq_new(foo),inseq(foo), replications = 10)