context("integration tests")

## use similar parameters as in paper
tr <- 0.2
sr <- 0.05
to <- 0.21
so <- 0.05
x <- 1
y <- 1

## integration tests
test_that("joint posterior of theta and alpha integrates to one", {
    ## define nested integration function
    .intFunAlpha <- function(alpha, tr, sr, to, so, x, y) {
        intFunTheta <- function(theta) {
            postPP(theta = theta, alpha = alpha, tr = tr, sr = sr, to = to,
                   so = so, x = x, y = y)
        }
        res <- stats::integrate(f = intFunTheta, lower = -Inf, upper = Inf)$value
        return(res)
    }
    ## vectorize function (expected by integrate)
    intFunAlpha <- Vectorize(FUN = .intFunAlpha)
    ## integrate nested integration function
    res <- stats::integrate(f = intFunAlpha, lower = 0, upper = 1, tr = tr, sr = sr, to = to,
                            so = so, x = x, y = y)
    expect_equal(object = res$value, expected = 1)
})


test_that("marginal posterior of alpha integrates to one", {
    ## integrate marginal posterior density of alpha
    res <- stats::integrate(f = postPPalpha, lower = 0, upper = 1, tr = tr,
                            sr = sr, to = to, so = so, x = x, y = y)
    expect_equal(object = res$value, expected = 1)
})


## ## TODO this doesn't work yet.... some parts where integral divergent
## thetaseq <- seq(-3, 3, length.out = 1000)
## postdens <- postPPtheta(theta = thetaseq, tr = tr, sr = sr, to = to, so = so,
##                         x = x, y = y)
## plot(thetaseq, postdens, type = "l")

## test_that("marginal posterior of theta integrates to one", {
##     ## integrate marginal posterior density of theta
##     res <- stats::integrate(f = postPPtheta, lower = -Inf, upper = Inf, tr = tr,
##                             sr = sr, to = to, so = so, x = x, y = y)
##     expect_equal(object = res$value, expected = 1)
## })
