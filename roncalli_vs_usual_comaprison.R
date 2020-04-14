p = 30
n = 1e5

foo <- scale(matrix(rnorm(n*p), nrow = n, ncol = p), scale = F)

sourceCpp("roncalli_vs_usual.cpp")
rbenchmark::benchmark(roncalli(foo,20), usual(foo, 20), replications = 10)



# 
# sourceCpp("inseq_fft_c++.cpp")
# sourceCpp("inseq.cpp")
# 
# a = initseq_new(foo)
# b = inseq(foo)
# 
# all.equal(a$Sig,b$Sig)
# 
# library(rbenchmark)
# benchmark(initseq_new(foo),inseq(foo), replications = 10)



