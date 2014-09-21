context("NIPALS CCA")

test_that("dimensions are correct",
    {
        n <- 100
        p <- 15
        q <- 10
        X <- matrix(rnorm(n*p), ncol = p)
        Y <- matrix(rnorm(n*q), ncol = q)
        res <- nipals(X, Y)
        expect_equal(nrow(res$a1), p)
        expect_equal(nrow(res$b1), q)
    })
