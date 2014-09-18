library(microbenchmark)
library(scca)
library(sccaf)

set.seed(42)
n <- 10
p <- 50
q <- 20
X <- matrix(rnorm(n*p),ncol=p)
Y <- matrix(rnorm(n*q),ncol=q)

r_res <- scca::NIPALS(X,Y)
cpp_res <- sccaf::nipals(X, Y)
cpp_res$a1 <- round(cpp_res$a1, 4)
cpp_res$b1 <- round(cpp_res$b1, 4)

bench <- microbenchmark(
    R_time = scca::NIPALS(X, Y),
    Cpp_time = sccaf::nipals(X, Y),
    times = 1000
    )
