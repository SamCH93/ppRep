library(tinytest)
library(ppRep)
library(cubature)

## use similar parameters as in paper
tr <- 0.2
sr <- 0.05
to <- 0.21
so <- 0.05
x <- 1
y <- 1

## test whether posterior integrates to 1
## -----------------------------------------------------------------------------

## joint posterior
.intFunAlpha <- function(alpha, tr, sr, to, so, x, y) {
    intFunTheta <- function(theta) {
        postPP(theta = theta, alpha = alpha, tr = tr, sr = sr, to = to,
               so = so, x = x, y = y)
    }
    res <- stats::integrate(f = intFunTheta, lower = -Inf, upper = Inf)$value
    return(res)
}
intFunAlpha <- Vectorize(FUN = .intFunAlpha)
res <- stats::integrate(f = intFunAlpha, lower = 0, upper = 1, tr = tr, sr = sr, to = to,
                        so = so, x = x, y = y)
expect_equal(res$value, 1,
             info = "joint posterior of theta and alpha integrates to one")

## marginal posterior of alpha
res2 <- stats::integrate(f = postPPalpha, lower = 0, upper = 1, tr = tr,
                         sr = sr, to = to, so = so, x = x, y = y)
expect_equal(res2$value, 1,
             info = "marginal posterior of alpha integrates to one")

## marginal posterior density of theta
intFun <- function(theta) {
    postPPtheta(theta = theta, tr = tr, sr = sr, to = to, so = so, x = x,
                y = y, rel.tol = .Machine$double.eps^0.75)
}
## need to use cubature library because stats::integrate can't handle this integral
res3 <- hcubature(f = intFun, lowerLimit = -Inf, upperLimit = Inf, tol = 1e-12)$integral
expect_equal(res3, 1,
             info = "marginal posterior of theta integrates to one")
## thetaseq <- seq(-3, 3, length.out = 1000)
## postdens <- postPPtheta(theta = thetaseq, tr = tr, sr = sr, to = to, so = so,
##                         x = x, y = y, rel.tol = .Machine$double.eps^0.5)
## plot(thetaseq, postdens, type = "l")
