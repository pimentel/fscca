library(microbenchmark)
library(scca)
library(fscca)

set.seed(42)

n <- 10
n <- 1000
p <- 50
q <- 20
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*q), ncol=q)

X <- scale(X, scale = T, center = F)
Y <- scale(Y, scale = T, center = F)

r_res <- scca::NIPALS(X,Y)
cpp_res <- fscca::nipals(X, Y)
cpp_res$a1 <- round(cpp_res$a1, 4)
cpp_res$b1 <- round(cpp_res$b1, 4)

all.equal(r_res, cpp_res)

bench <- microbenchmark(
    R_time = scca::NIPALS(X, Y),
    Cpp_time = fscca::nipals(X, Y),
    CCA_time = CCA::cc(X, Y),
    times = 1000
    )

print(bench)

n <- 20
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*q), ncol=q)

bench <- microbenchmark(
    R_time = scca::NIPALS(X, Y),
    Cpp_time = fscca::nipals(X, Y),
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

################################################################################

hi <- sparse_nipals(X[1:10,], Y[1:10,], 2, 3)
bye <- NIPALS.sparse(X[1:10, ], Y[1:10,], 2, 3, "LASSO")

bench_sparse <- microbenchmark(
    R_time = scca::NIPALS.sparse(X[1:10,], Y[1:10,], 2, 3, "LASSO"),
    Cpp_time = fscca::sparse_nipals(X[1:10,], Y[1:10,], 2, 3),
    times = 1000
    )

sum(abs(hi$a))
sum(abs(bye$a))

sum(abs(hi$b))
sum(abs(bye$b))

s_res <- scca(X,Y, "LASSO")
all.equal(s_res$U[,1], as.numeric(X %*% s_res$A[,1]))
s_res_1 <- NIPALS.sparse(X, Y, 1, 3, "LASSO")
s_res_3 <- NIPALS.sparse(X, Y, 1, 3, "LASSO")
s_res_05 <- NIPALS.sparse(X, Y, 0.000005, 0.005, "LASSO")
c(sum(abs(s_res_1$a1)), sum(abs(s_res_1$b1)))
c(sum(abs(s_res_3$a1)), sum(abs(s_res_3$b1)))
c(sum(abs(s_res_05$a1)), sum(abs(s_res_05$b1)))

c_res_1 <- sparse_nipals(X, Y, 1, 3)
c_res_3 <- sparse_nipals(X, Y, 3, 3)
c_res_10 <- sparse_nipals(X, Y, 10, 3)
c_res_40 <- sparse_nipals(X, Y, 10, 40)
c_res_05 <- sparse_nipals(X, Y, 0.5, 0.5)
c_res_05 <- sparse_nipals(X, Y, 0.000005, 0.005)

c(sum(abs(c_res_1$a)), sum(abs(c_res_1$b)))
c(sum(abs(c_res_3$a)), sum(abs(c_res_3$b)))
c(sum(abs(c_res_10$a)), sum(abs(c_res_10$b)))
c(sum(abs(c_res_40$a)), sum(abs(c_res_40$b)))
c(sum(abs(c_res_05$a)), sum(abs(c_res_05$b)))

p_res_05 <- PMA::CCA(X, Y, penaltyx = 0.5, penaltyz = 0.5)
p_res_005 <- PMA::CCA(X, Y, penaltyx = 0.05, penaltyz = 0.05)
p_res_01 <- PMA::CCA(X, Y, penaltyx = 0.01, penaltyz = 0.01)
p_res_04 <- PMA::CCA(X, Y, penaltyx = 0.04, penaltyz = 0.04)
sqrt(sum(c_res$a^2))
