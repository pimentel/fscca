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

hi <- sparse_nipals(X[1:10,], Y[1:10,], "lasso", "lasso", 2, 3)
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
s_res_3 <- NIPALS.sparse(X, Y, 3, 3, "LASSO")
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

# PMA simulations
u <- matrix(c(rep(1,25),rep(0,75)),ncol=1)
v1 <- matrix(c(rep(1,50),rep(0,450)),ncol=1)
v2 <- matrix(c(rep(0,50),rep(1,50),rep(0,900)),ncol=1)
x <- u%*%t(v1) + matrix(rnorm(100*500),ncol=500)
z <- u%*%t(v2) + matrix(rnorm(100*1000),ncol=1000)

meow <- sparse_nipals(x, z, 1, 1)
meow <- sparse_nipals(x, z, 10, 20)
c(sum(abs(meow$a)), sum(abs(meow$b)))
c(sum(abs(out$u[,1])), sum(abs(out$v[,1]))
hist(abs(out$u[,1]))
hist(abs(meow$a))


hi <- NIPALS.sparse(X, Y, 40000000, 40000000, "LASSO")
hi <- NIPALS.sparse(X, Y, 2.5, 2.5, "LASSO")
c(sum(abs(hi$a)), sum(abs(hi$b)))
bye <- sparse_nipals(X, Y, 2.5, 2.5)
bye <- sparse_nipals(X, Y, 40000000, 40000000)
c(sum(abs(bye$a)), sum(abs(bye$b)))

t(X %*% bye$a) %*% (Y %*% bye$b)
t(X %*% hi$a) %*% (Y %*% hi$b)


set.seed(42)
n <- 1000

p <- 50
q <- 20
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*q), ncol=q)

X <- scale(X, scale = T, center = F)
Y <- scale(Y, scale = T, center = F)

bench <- microbenchmark(
    R_time = scca::NIPALS.sparse(X, Y, 2.0, 3.0, "LASSO"),
    Cpp_time = fscca::sparse_nipals(X, Y, "lasso", "lasso", 2.0, 3.0),
    times = 1000
    )

# Unit: milliseconds
#      expr      min       lq   median       uq       max neval
#    R_time 45.70580 47.67376 57.70300 64.41011 259.85699  1000
#  Cpp_time 13.15465 13.30561 13.44036 14.95119  33.17595  1000
# Unit: milliseconds
#      expr      min       lq   median       uq       max neval
#    R_time 45.59999 46.62820 47.82320 59.41712 232.38032  1000
#  Cpp_time 13.12863 13.23349 13.29242 13.35226  17.89528  1000



set.seed(42)
n <- 40

p <- 2000
q <- 2200
X <- matrix(rnorm(n*p), ncol=p)
Y <- matrix(rnorm(n*q), ncol=q)

X <- scale(X, scale = T, center = F)
Y <- scale(Y, scale = T, center = F)


bench <- microbenchmark(
    R_time = scca::NIPALS.sparse(X, Y, 3.0, 3.0, "LASSO"),
    Cpp_time = fscca::sparse_nipals(X, Y, "lasso", "lasso", 3.0, 3.0),
    times = 50
    )

z <- sample.int(40, 32)

bench <- microbenchmark(
    R = X[z,],
    Cpp = get_submatrix(X, z),
    times = 10000)

ggplot(bench, aes(expr, time)) + geom_boxplot()
    + ylim(0, 600000)

m <- rnorm(ncol(X))

bench <- microbenchmark(
    R = X[z,] %*% m,
    Cpp = get_submatrix_mult(X, z, m),
    Cpp_ptr = get_submatrix_mult_ptr(X, z, m),
    times = 10000)

ggplot(bench, aes(expr, time)) + geom_boxplot() + ylim(0, 800000)


# testing shuffle
x <- 1:4
shuf_dist <- lapply(1:1000, function(i)
    {
        shuf_gen <- lapply(1:100000, function(it)
            {
                paste(shuffle(x), collapse = "")
            })
        table(unlist(shuf_gen))
    })

s_d <- rbindlist(lapply(shuf_dist, data.frame))
ggplot(s_d, aes(factor(Var1), Freq)) + geom_boxplot()


# split into k-groups

kfold <- groups_to_rows(split_in_groups(10, 3), 3)
