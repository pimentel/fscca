library(microbenchmark)
library(scca)
library(sccaf)

set.seed(42)

n <- 1000
p <- 50
q <- 20
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*q), ncol=q)

X <- scale(X, scale = T, center = F)
Y <- scale(Y, scale = T, center = F)

r_res <- scca::NIPALS(X,Y)
cpp_res <- sccaf::nipals(X, Y)
cpp_res$a1 <- round(cpp_res$a1, 4)
cpp_res$b1 <- round(cpp_res$b1, 4)

all.equal(r_res, cpp_res)

bench <- microbenchmark(
    R_time = scca::NIPALS(X, Y),
    Cpp_time = sccaf::nipals(X, Y),
    CCA_time = CCA::cc(X, Y),
    times = 1000
    )

print(bench)

n <- 20
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*q), ncol=q)

bench <- microbenchmark(
    R_time = scca::NIPALS(X, Y),
    Cpp_time = sccaf::nipals(X, Y),
    times = 1000
    )

print(bench)

ggplot(bench, aes(expr, time, fill = expr)) + geom_boxplot() +
    scale_y_log10()

# Unit: milliseconds
#      expr       min        lq    median        uq       max neval
#    R_time 12.124802 12.485295 13.495677 13.879834 139.45854  1000
#  Cpp_time  8.232687  8.380976  8.692669  8.863921  12.58179  1000
#  CCA_time 52.506479 54.135987 54.677352 56.211481 178.46186  1000

computeCov <- function(a, b, X, Y)
{
    crossprod(a, t(X)) %*% (Y %*% b) / (nrow(X) - 1)
}


library(CCA)
cc_res <- cc(X, Y)

computeCov(cc_res$xcoef[,1], cc_res$ycoef[,1], X, Y)
computeCov(cpp_res$a1, cpp_res$b1, X, Y)
computeCov(r_res$a1, r_res$b1, X, Y)

