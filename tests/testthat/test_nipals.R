context("NIPALS CCA")

n <- 100
p <- 15
q <- 10
X <- matrix(rnorm(n*p), ncol = p)
Y <- matrix(rnorm(n*q), ncol = q)

test_that("input dimensions are ok",
    {
        expect_that(nipals(X[1:10, ], Y[1:3,]), throws_error())
    })

test_that("out dimensions are correct",
    {
        res <- nipals(X, Y)
        expect_equal(nrow(res$a1), p)
        expect_equal(nrow(res$b1), q)
    })
