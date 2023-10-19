#' @title Variance of normalized power prior
#' @description This function computes the variance of a normalized power prior
#'     conditional on a fixed power parameter and an initial normal prior for
#'     the effect size.
#' @param vardata Variance of the data.
#' @param priorvar Variance parameter of initial normal prior. Defaults to
#'     \code{Inf} (uniform prior).
#' @param alpha Power parameter. Indicates to which power the likelihood of the
#'     data is raised. Can be set to a number in [0, 1]. Defaults to \code{1}.
#'
#' @return Prior variance
#'
#' @author Samuel Pawel
#'
#' @keywords internal
postNormVar <- function(vardata, priorvar, alpha = 1) {
    ## separate in two cases to avoid some numerical problems in special
    ## situation when flat initial prior (v -> Inf) and alpha = 0
    if (!is.finite(priorvar)) {
        postvar <- vardata/alpha
    } else {
        postvar <- 1/(alpha/vardata + 1/priorvar)
    }
    return(postvar)
}

#' @title Mean of normalized power prior
#' @description This function computes the mean of a normalized power prior
#'     conditional on a fixed power parameter and an initial normal prior for
#'     the effect size.
#' @param dat Data.
#' @param vardata Variance of the data.
#' @param priormean Mean parameter of initial normal prior. Defaults to
#'     \code{0}.
#' @param priorvar Variance parameter of initial normal prior. Defaults to
#'     \code{Inf} (uniform prior).
#' @param alpha Power parameter. Indicates to which power the likelihood of the
#'     data is raised. Can be set to a number in [0, 1]. Defaults to \code{1}.
#'
#' @return Prior mean
#'
#' @author Samuel Pawel
#'
#' @keywords internal
postNormMean <- function(dat, vardata, priormean, priorvar, alpha = 1) {
    ## separate in two cases to avoid some numerical problems in special
    ## situation when flat initial prior (v -> Inf) and alpha = 0
    if (!is.finite(priorvar)) {
        postmean <- dat
    } else {
        postvar <- postNormVar(vardata = vardata, priorvar = priorvar,
                               alpha = alpha)
        postmean <- postvar*(dat*alpha/vardata + priormean/priorvar)
    }
    return(postmean)
}


#' @title Posterior density of effect size and power parameter
#'
#' @description This function computes the posterior density of effect size
#'     \eqn{\theta}{theta} and power parameter \eqn{\alpha}{alpha} assuming a
#'     normal likelihood for original and replication effect estimate. A power
#'     prior for \eqn{\theta}{theta} is constructed by updating an initial
#'     normal prior \eqn{\theta \sim \mathrm{N}(\code{m}, \code{v})}{theta ~
#'     N(m, v)} with the likelihood of the original data raised to the power of
#'     \eqn{\alpha}{alpha}. A marginal beta prior \eqn{\alpha \sim
#'     \mbox{Beta}(\code{x},\code{y})}{alpha ~ Beta(x, y)} is assumed.
#'
#' @param theta Effect size. Has to be of length one or the same length as
#'     \code{alpha}.
#' @param alpha Power parameter. Has to be of length one or the same length as
#'     \code{theta}.
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to \code{1}.
#' @param y Number of failures parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to \code{1}.
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{0}.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{Inf} (uniform prior).
#' @param ... Additional arguments passed to \code{stats::integrate}.
#'
#' @return Posterior density
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{postPPalpha}}, \code{\link{postPPtheta}}, \code{\link{plotPP}}
#'
#' @examples
#' alpha <- seq(0, 1, length.out = 200)
#' theta <- seq(0, 0.3, length.out = 200)
#' parGrid <- expand.grid(alpha = alpha, theta = theta)
#' postdens <- postPP(theta = parGrid$theta, alpha = parGrid$alpha, tr = 0.1,
#'                    sr = 0.05, to = 0.2, so = 0.05)
#' postdensMat <- matrix(data = postdens, ncol = 200, byrow = TRUE)
#' filled.contour(x = theta, y = alpha, z = postdensMat,
#'                xlab = bquote("Effect size" ~ theta),
#'                ylab = bquote("Power parameter" ~ alpha), nlevels = 15,
#'                color.palette = function(n) hcl.colors(n = n, palette = "viridis"))
#' @export
postPP <- function(theta, alpha, tr, sr, to, so, x = 1, y = 1, m = 0, v = Inf,
                   ...) {
    ## input checks
    stopifnot(
        1 <= length(theta),
        1 <= length(alpha),
        ((length(theta) == length(alpha) |
          (length(theta) == 1 & length(alpha) > 1) |
          (length(theta) > 1 & length(alpha) == 1))),

        any(!is.numeric(theta)) == FALSE,
        any(!is.finite(theta)) == FALSE,

        any(!is.numeric(alpha)) == FALSE,
        any(!is.finite(alpha)) == FALSE,
        any(!(0 <= alpha)) == FALSE,
        any(!(alpha <= 1)) == FALSE,

        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 <= x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 <= y,

        length(m) == 1,
        is.numeric(m),
        is.finite(m),

        length(v) == 1,
        is.numeric(v),
        0 < v
    )

    ## compute normalizing constant just once
    nC <- margLik(tr = tr, sr = sr, to = to, so = so, x = x, y = y, m = m,
                  v = v, ... = ...)
    if (is.nan(nC)) {
        out <- rep(x = NaN, times = pmax(length(alpha), length(theta)))
        return(out)
    }

    ## compute mean and variance of normalized power prior conditional on alpha
    pvar <- postNormVar(vardata = so^2, priorvar = v, alpha = alpha)
    pmean <- postNormMean(dat = to, vardata = so^2, priormean = m, priorvar = v,
                          alpha = alpha)

    ## compute posterior density
    densProp <- stats::dnorm(x = tr, mean = theta, sd = sr) *
        stats::dnorm(x = theta, mean = pmean, sd = sqrt(pvar)) *
        stats::dbeta(x = alpha, shape1 = x, shape2 = y)
    dens <- densProp / nC
    return(dens)
}

#' @rdname posterioralpha
#' @title Marginal posterior distribution of power parameter
#'
#' @description These functions compute the marginal posterior of the power
#'     parameter \eqn{\alpha}{alpha}. A power prior for \eqn{\theta}{theta} is
#'     constructed by updating an initial normal prior \eqn{\theta \sim
#'     \mathrm{N}(\code{m}, \code{v})}{theta ~ N(m, v)} with the likelihood of
#'     the original data raised to the power of \eqn{\alpha}{alpha}. A marginal
#'     beta prior \eqn{\alpha \sim \mbox{Beta}(\code{x},\code{y})}{alpha ~
#'     Beta(x, y)} is assumed.
#'
#' @param alpha Power parameter. Can be a vector.
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter of beta prior \eqn{\alpha}{alpha}.
#'     Defaults to \code{1}.
#' @param y Number of failures parameter of beta prior \eqn{\alpha}{alpha}.
#'     Defaults to \code{1}.
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{0}.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{Inf} (uniform prior).
#' @param ... Additional arguments passed to \code{stats::integrate}.
#'
#' @return \code{postPPalpha} returns the marginal posterior density of the power
#'     parameter.
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{postPP}}, \code{\link{postPPtheta}}, \code{\link{plotPP}}
#'
#' @examples
#' alpha <- seq(0, 1, 0.001)
#' margpostdens <- postPPalpha(alpha = alpha, tr = 0.1, to = 0.2, sr = 0.05, so = 0.05)
#' plot(alpha, margpostdens, type = "l", xlab = bquote("Power parameter" ~ alpha),
#'      ylab = "Marginal posterior density", las = 1)
#' @export
postPPalpha <- function(alpha, tr, sr, to, so, x = 1, y = 1, m = 0, v = Inf,
                        ...) {
    ## input checks
    stopifnot(
        1 <= length(alpha),
        any(!is.numeric(alpha)) == FALSE,
        any(!is.finite(alpha)) == FALSE,
        any(!(0 <= alpha)) == FALSE,
        any(!(alpha <= 1)) == FALSE,

        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 <= x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 <= y,

        length(m) == 1,
        is.numeric(m),
        is.finite(m),

        length(v) == 1,
        is.numeric(v),
        0 < v
    )

    ## compute normalizing constant just once
    nC <- margLik(tr = tr, sr = sr, to = to, so = so, x = x, y = y, m = m,
                  v = v, ... = ...)
    if (is.nan(nC)) {
        out <- rep(x = NaN, times = length(alpha))
        return(out)
    }

    ## compute mean and variance of normalized power prior conditional on alpha
    pvar <- postNormVar(vardata = so^2, priorvar = v, alpha = alpha)
    pmean <- postNormMean(dat = to, vardata = so^2, priormean = m, priorvar = v,
                          alpha = alpha)

    ## compute marginal posterior density
    margdensProp <- stats::dnorm(x = tr, mean = pmean, sd = sqrt(sr^2 + pvar)) *
        stats::dbeta(x = alpha, shape1 = x, shape2 = y)
    margdens <- margdensProp / nC
    return(margdens)
}

#' @rdname posterioralpha
#'
#' @param level Credibility level of the highest posterior density interval.
#'     Defaults to \code{0.95}.
#'
#' @return \code{postPPalphaHPD} returns the highest marginal posterior density
#'     interval of the power parameter.
#' @export
postPPalphaHPD <- function(level = 0.95, tr, sr, to, so, x = 1, y = 1, m = 0,
                           v = Inf, ...) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1
    )
    ## posterior quantile function
    quantileFun <- function(q, ...) {
        if (q == 0) {
            res <- 0
        } else if (q == 1) {
            res <- 1
        } else {
            densFun <- function(alpha) {
                postPPalpha(alpha = alpha, tr = tr, sr = sr, to = to, so = so,
                            x = x, y = y, m = m, v = v, ... = ...)
            }
            rootFun <- function(x) {
                stats::integrate(f = densFun, lower = 0, upper = x,
                                 ... = ...)$value - q
            }
            res <- stats::uniroot(f = rootFun, interval = c(0, 1))$root
        }
        return(res)
    }

    ## find narrowest interval
    optFun. <- function(qLow) {
        width <- quantileFun(q = qLow + level) - quantileFun(q = qLow)
        return(width)
    }
    optFun <- Vectorize(FUN = optFun.)
    minLower <- try(stats::optim(par = (1 - level)/2, fn = optFun,
                                 method = "L-BFGS-B", lower = 0,
                                 upper = 1 - level)$par)
    if (inherits(minLower, "try-error")) {
        CI <- c("lower" = NaN, "upper" = NaN)
    } else {
        CI <- c("lower" = quantileFun(q = minLower),
                "upper" = quantileFun(q = minLower + level))
    }
    return(CI)
}

## ## example
## tr <- 0.21
## sr <- 0.05
## to <- 0.3
## so <- 0.04
## x <- 1
## y <- 1
## m <- 0
## v <- Inf
## ci <- postPPalphaHPD(level = 0.95, tr = tr, sr = sr, to = to, so = so, x = x,
##                      y = y, m = m, v = v)
## aseq <- seq(0, 1, 0.001)
## aFun <- function(a) {
##     postPPalpha(alpha = a, tr = tr, sr = sr, to = to, so = so, x = x, y = y,
##                 m = m, v = v)
## }
## integrate(aFun, lower = ci[1], upper = ci[2])
## amax <- optim(par = 0, fn = function(x) -log(aFun(x)),
##               lower = 0, upper = 1, method = "Brent")$par
## plot(aseq, aFun(aseq), type = "l", xlab = bquote(alpha), ylab = "Density",
##      ylim = c(0, aFun(amax)*1.05))
## arrows(x0 = ci[1], x1 = ci[2], y0 = aFun(amax)*1.05, angle = 90, code = 3, length = 0.05)


#' @rdname posteriortheta

#' @title Marginal posterior distribution of effect size
#'
#' @description These functions compute the marginal posterior of the effect
#'     size \eqn{\theta}{theta}. A power prior for \eqn{\theta}{theta} is
#'     constructed by updating an initial normal prior \eqn{\theta \sim
#'     \mathrm{N}(\code{m}, \code{v})}{theta ~ N(m, v)} with likelihood of the
#'     original data raised to the power of \eqn{\alpha}{alpha}. The power
#'     parameter \eqn{\alpha}{alpha} can either be fixed to some value between 0
#'     and 1 or it can have a beta prior distribution \eqn{\alpha \sim
#'     \mbox{Beta}(\code{x}, \code{y})}{alpha ~ Beta(x,y)}.
#'
#' @param theta Effect size. Can be a vector.
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter for beta prior of power parameter
#'     \eqn{\alpha}{alpha}. Defaults to \code{1}. Is only taken into account
#'     when \code{alpha = NA}.
#' @param y Number of failures parameter for beta prior of power parameter
#'     \eqn{\alpha}{alpha}. Defaults to \code{1}. Is only taken into account
#'     when \code{alpha = NA}.
#' @param alpha Power parameter. Can be set to a number between 0 and 1.
#'     Defaults to \code{NA} (a beta prior on the power parameter).
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{0}.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{Inf} (uniform prior).
#' @param hypergeo Logical indicating whether for uniform priors, the marginal
#'     posterior should be computed with the hypergeometric function. Defaults
#'     to \code{FALSE} (using numerical integration instead).
#' @param ... Additional arguments passed to \code{stats::integrate} or
#'     \code{hypergeo::genhypergeo} (depending on the \code{hypergeo} argument).
#'
#' @return \code{postPPtheta} returns the marginal posterior density of the
#'     effect size.
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{postPP}}, \code{\link{postPPalpha}}, \code{\link{plotPP}}
#'
#' @examples
#' theta <- seq(0, 0.6, 0.001)
#' margpostdens <- postPPtheta(theta = theta, tr = 0.1, to = 0.2, sr = 0.05, so = 0.05)
#' plot(theta, margpostdens, type = "l", xlab = bquote("Effect size" ~ theta),
#'      ylab = "Marginal posterior density", las = 1)
#' @export
postPPtheta <- function(theta, tr, sr, to, so, x = 1, y = 1, alpha = NA, m = 0,
                        v = Inf, hypergeo = FALSE, ...) {
    ## input checks
    stopifnot(
        1 <= length(theta),
        any(!is.numeric(theta)) == FALSE,
        ## any(!is.finite(theta)) == FALSE,

        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 <= x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 <= y,

        length(alpha) == 1,
        (is.na(alpha) |
         ((is.numeric(alpha)) &
          (is.finite(alpha)) &
          (0 <= alpha) &
          (alpha <= 1))),

        length(m) == 1,
        is.numeric(m),
        is.finite(m),

        length(v) == 1,
        is.numeric(v),
        0 < v,

        length(hypergeo) == 1,
        is.logical(hypergeo),
        !is.na(hypergeo)
    )

    ## alpha fixed
    if (!is.na(alpha)) {
        ## compute mean and variance of normalized power prior conditional on
        ## alpha
        pvar <- postNormVar(vardata = so^2, priorvar = v, alpha = alpha)
        pmean <- postNormMean(dat = to, vardata = so^2, priormean = m,
                              priorvar = v, alpha = alpha)
        ## mean and variance of posterior
        postVar <- 1/(1/sr^2 + 1/pvar)
        postMean <- (tr/sr^2 + pmean/pvar)*postVar
        margdens <- stats::dnorm(x = theta, mean = postMean, sd = sqrt(postVar))
    } else { ## alpha random

        ## compute normalizing constant just once
        nC <- margLik(tr = tr, sr = sr, to = to, so = so, x = x, y = y, m = m,
                      v = v, ... = ...)
        if (is.nan(nC)) {
            out <- rep(x = NaN, times = length(theta))
            return(out)
        }

        ## compute marginal posterior
        ## 1) when flat prior for effect size (v = Inf), can compute using
        ## confluent hypergeometric function
        if (!is.finite(v) & hypergeo == TRUE) {
            margdens <- stats::dnorm(x = tr, mean = theta, sd = sr) *
                abs(hypergeo::genhypergeo(U = x + 0.5, L = x + y + 0.5,
                                          z = -(to - theta)^2/(2*so^2),
                                          ... = ...)) *
                beta(a = x + 0.5, b = y) / nC / sqrt(2*pi*so^2) /
                beta(a = x, b = y)
        } else {
            ## 2) otherwise use numerical integration
            margdens <- vapply(X = theta, FUN = function(thetai) {
                ## integrate out power parameter alpha
                intFun <- function(alpha) {
                    ## compute mean and variance of normalized power prior
                    ## conditional on alpha
                    pvar <- postNormVar(vardata = so^2, priorvar = v, alpha = alpha)
                    pmean <- postNormMean(dat = to, vardata = so^2,
                                          priormean = m, priorvar = v,
                                          alpha = alpha)
                    stats::dnorm(x = thetai, mean = pmean, sd = sqrt(pvar)) *
                        stats::dbeta(x = alpha, shape1 = x, shape2 = y)
                }
                int <- try(stats::integrate(f = intFun, lower = 0, upper = 1,
                                            ... = ...)$value)
                if (inherits(int, "try-error")) {
                    margdens_i <- NaN
                    warnString <- paste("Numerical problems integrating out power parameter",
                                        "from posterior. \nTry adjusting integration options",
                                        "with ... argument. \nSee ?stats::integrate for",
                                        "available options.")
                    warning(warnString)
                }
                else {
                    margdens_i <- stats::dnorm(x = tr, mean = thetai, sd = sr) * int / nC
                }
                return(margdens_i)
            }, FUN.VALUE = 1)
        }
    }
    return(margdens)
}



#' @rdname posteriortheta
#'
#' @param level Credibility level of the highest posterior density interval.
#'     Defaults to \code{0.95}.
#' @param thetaRange The numerical search range for the effect size. Defaults to
#'     the \code{level*100}\% confidence inteval range inflated by a factor of
#'     three. We recommend changing this argument only if there are numerical
#'     problems in calculating the HPD interval.
#' @param quantileRange The numerical search range for the lower posterior
#'     quantile of the HPD interval. Defaults to the range between \code{(1 -
#'     level)*0.2} and \code{(1 - level)*0.8}. We recommend changing this
#'     argument only if there are numerical problems in calculating the HPD
#'     interval.
#'
#' @return \code{postPPthetaHPD} returns the highest marginal posterior density
#'     interval of the effect size (this may take a while).
#' @export
postPPthetaHPD <- function(level, tr, sr, to, so, x = 1, y = 1, alpha = NA,
                           m = 0, v = Inf,
                           thetaRange = tr + c(-1, 1)*stats::qnorm(p = (1 + level)/2)*sr*3,
                           quantileRange = c((1 - level)*0.2, (1 - level)*0.8),
                           ...) {
    ## input checks
    stopifnot(
        length(level) == 1,
        is.numeric(level),
        is.finite(level),
        0 < level, level < 1,

        length(thetaRange) == 2,
        all(is.numeric(thetaRange)),
        all(is.finite(thetaRange)),
        thetaRange[1] < thetaRange[2],

        length(quantileRange) == 2,
        all(is.numeric(quantileRange)),
        all(is.finite(quantileRange)),
        all(quantileRange <= 1 - level),
        all(quantileRange >= 0),
        quantileRange[1] < quantileRange[2]

    )
    ## posterior quantile function
    quantileFun <- function(q, ...) {
        if (q == 0) {
            res <- -Inf
        } else if (q == 1) {
            res <- Inf
        } else {
            densFun <- function(theta) {
                postPPtheta(theta = theta, tr = tr, sr = sr, to = to, so = so,
                            x = x, y = y, alpha = alpha, m = m, v = v,
                            ## this integral requires typically higher precision
                            rel.tol = .Machine$double.eps^0.75)
            }
            rootFun <- function(x) {
                stats::integrate(f = densFun, lower = -Inf, upper = x,
                                 ... = ...)$value - q
            }
            res <- stats::uniroot(f = rootFun, interval = thetaRange)$root
        }
        return(res)
    }

    ## find narrowest interval
    optFun. <- function(qLow) {
        width <- quantileFun(q = qLow + level) - quantileFun(q = qLow)
        return(width)
    }
    optFun <- Vectorize(FUN = optFun.)
    minLower <- try(stats::optim(par = (1 - level)/2, fn = optFun,
                                 method = "L-BFGS-B", lower = quantileRange[1],
                                 upper = quantileRange[2])$par)
    if (inherits(minLower, "try-error")) {
        CI <- c("lower" = NaN, "upper" = NaN)
    } else {
        CI <- c("lower" = quantileFun(q = minLower),
                "upper" = quantileFun(q = minLower + level))
    }
    return(CI)
}

## ## example for theta
## tr <- 0.21
## sr <- 0.05
## to <- 0.3
## so <- 0.04
## x <- 1
## y <- 1
## m <- 0
## v <- Inf
## ci <- postPPthetaHPD(level = 0.95, tr = tr, sr = sr, to = to, so = so, x = x,
##                      y = y, m = m, v = v)
## tseq <- seq(tr - 5*sr, tr + 5*sr, length.out = 200)
## tFun <- function(t) {
##     postPPtheta(theta = t, tr = tr, sr = sr, to = to, so = so, x = x, y = y,
##                 m = m, v = v,
##                 rel.tol = .Machine$double.eps^0.75)
## }
## integrate(tFun, lower = ci[1], upper = ci[2])
## tmax <- optim(par = mean(tseq), fn = function(x) -log(tFun(x)),
##               lower = min(tseq), upper = max(tseq), method = "Brent")$par
## plot(tseq, tFun(tseq), type = "l", xlab = bquote(theta), ylab = "Density")
## arrows(x0 = ci[1], x1 = ci[2], y0 = tFun(tmax)*1.01, angle = 90, code = 3, length = 0.05)



## tr <- 0.1
## to <- 0.2
## sr <- so <- 0.05
## x <- y <- 1
## theta <- seq(-10, 10, length.out = 5000)
## margpostdens <- postPPtheta(theta = theta, tr = tr, to = to, sr = sr, so = so,
##                             hypergeo = FALSE,
##                             rel.tol = .Machine$double.eps^0.75)
## plot(theta, margpostdens, type = "l", xlab = bquote("Effect size" ~ theta),
##      ylab = "Marginal posterior density", las = 1)


#' @title Plot joint and marginal posterior distributions
#'
#' @description This convenience function computes and, if desired, visualizes
#'     the joint posterior density of effect size \eqn{\theta}{theta} and power
#'     parameter \eqn{\alpha}{alpha}, as well as the marginal posterior
#'     densities of effect size \eqn{\theta}{theta} and power parameter
#'     \eqn{\alpha}{alpha} individually. See the functions \code{\link{postPP}},
#'     \code{\link{postPPalpha}}, and \code{\link{postPPtheta}} for more details
#'     on their computation.
#'
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to \code{1}.
#' @param y Number of failures parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to \code{1}.
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{0}.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to \code{Inf} (uniform prior).
#' @param thetaRange Range of effect sizes. Defaults to three standard errors
#'     around the replication effect estimate.
#' @param alphaRange Range of power parameters. Defaults to the range between
#'     zero and one.
#' @param nGrid Number of grid points. Defaults to \code{100}.
#' @param plot Logical indicating whether data should be plotted. If
#'     \code{FALSE} only the data used for plotting are returned.
#' @param CI Logical indicating whether 95\% highest posterior credible
#'     intervals should be plotted. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to \code{stats::integrate} for
#'     computation of posterior densities and highest posterior density credible
#'     intervals.
#'
#' @return Plots joint and marginal posterior densities, invisibly returns a
#'     list with the data for the plots.
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{postPP}}, \code{\link{postPPalpha}}, \code{\link{postPPtheta}}
#'
#' @examples
#' plotPP(tr = 0.2, sr = 0.05, to = 0.15, so = 0.05)
#' @export

plotPP <- function(tr, sr, to, so, x = 1, y = 1, m = 0, v = Inf,
                   thetaRange = c(tr - 3*sr, tr + 3*sr),
                   alphaRange = c(0, 1), nGrid = 100, plot = TRUE,
                   CI = FALSE, ...) {
    ## input checks
    stopifnot(
        length(tr) == 1,
        is.numeric(tr),
        is.finite(tr),

        length(to) == 1,
        is.numeric(to),
        is.finite(to),

        length(sr) == 1,
        is.numeric(sr),
        is.finite(sr),
        0 < sr,

        length(so) == 1,
        is.numeric(so),
        is.finite(so),
        0 < so,

        length(x) == 1,
        is.numeric(x),
        is.finite(x),
        0 <= x,

        length(y) == 1,
        is.numeric(y),
        is.finite(y),
        0 <= y,

        length(m) == 1,
        is.numeric(m),
        is.finite(m),

        length(v) == 1,
        is.numeric(v),
        0 < v,

        length(alphaRange) == 2,
        all(is.numeric(alphaRange)),
        all(is.finite(alphaRange)),
        alphaRange[2] > alphaRange[1],
        0 <= alphaRange[1],
        alphaRange[2] <= 1,

        length(thetaRange) == 2,
        all(is.numeric(thetaRange)),
        all(is.finite(thetaRange)),
        alphaRange[2] > alphaRange[1],

        length(nGrid) == 1,
        is.numeric(nGrid),
        is.finite(nGrid),
        0 < nGrid,

        length(plot) == 1,
        is.logical(plot),
        !is.na(plot),

        length(CI) == 1,
        is.logical(CI),
        !is.na(CI)
    )

    ## grids for computing the densities
    alphaGrid <- seq(alphaRange[1], alphaRange[2], length.out = nGrid)
    thetaGrid <- seq(thetaRange[1], thetaRange[2], length.out = nGrid)
    jointGrid <- expand.grid(alpha = alphaGrid, theta = thetaGrid)

    ## compute densities
    jointdens <- postPP(theta = jointGrid$theta, alpha = jointGrid$alpha, tr = tr,
                        sr = sr, to = to, so = so, x = x, y = y, m = m, v = v, ... = ...)
    alphadens <- postPPalpha(alpha = alphaGrid, tr = tr, sr = sr, to = to,
                             so = so, x = x, y = y, m = m, v = v, ... = ...)
    thetadens <- postPPtheta(theta = thetaGrid, tr = tr, sr = sr, to = to,
                             so = so, x = x, y = y, m = m, v = v, ... = ...)

    ## compute HPD intervals
    if (CI == TRUE) {
        alphaCI <- postPPalphaHPD(level = 0.95, tr = tr, sr = sr, to = to,
                                  so = so, x = x, y = y, m = m, v = v, ... = ...)
        thetaCI <- postPPthetaHPD(level = 0.95, tr = tr, sr = sr, to = to,
                                  so = so, x = x, y = y, m = m, v = v, ... = ...)
    } else {
        alphaCI <- c(NA, NA)
        thetaCI <- c(NA, NA)
    }

    ## plot posterior distributions
    if (plot) {
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))
        graphics::layout(mat = matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE))
        graphics::par(mar = c(1, oldpar$mar[2:4]))
        ## joint distribution
        jointdensMat <- matrix(data = jointdens, ncol = nGrid, byrow = TRUE)
        graphics::image(x = thetaGrid, y = alphaGrid, z = jointdensMat,
                        xlab = "",
                        ylab = bquote("Power parameter" ~ alpha),
                        main = bquote(plain("Joint posterior density")),
                        col = grDevices::hcl.colors(n = 100, palette = "Blues 3",
                                                    rev = TRUE),
                        las = 1)
        graphics::mtext(text = bquote("Effect size" ~ theta), side = 1, line = 2.5, cex = 0.9)
        graphics::contour(x = thetaGrid, y = alphaGrid, z = jointdensMat, add = TRUE,
                          drawlabels = FALSE, nlevels = 5, col = "#00000080")
        graphics::par(mar = oldpar$mar)
        ## power parameter
        plot(alphaGrid, alphadens, xlab = bquote("Power parameter" ~ alpha),
             ylab = "Marginal posterior density", type = "l", las = 1,
             ylim = c(0, max(alphadens)*1.1))
        if (CI == TRUE) {
            graphics::arrows(x0 = alphaCI[1], x1 = alphaCI[2],
                             y0 = max(alphadens)*1.1, length = 0.05, angle = 90,
                             code = 3)
        }
        ## effect size
        plot(thetaGrid, thetadens, xlab = bquote("Effect size" ~ theta),
             ylab = "Marginal posterior density", type = "l", las = 1,
             ylim = c(0, max(thetadens)*1.1))
        if (CI == TRUE) {
            graphics::arrows(x0 = thetaCI[1], x1 = thetaCI[2],
                             y0 = max(thetadens)*1.1, length = 0.05, angle = 90,
                             code = 3)
        }
    }

    ## return plot data
    out <- list(jointDF = data.frame(jointGrid, density = jointdens),
                alphaDF = data.frame(alpha = alphaGrid,
                                     density = alphadens),
                thetaDF = data.frame(theta = thetaGrid,
                                     density = thetadens),
                ciDF = data.frame(lower = c(alphaCI[1], thetaCI[1]),
                                  upper = c(alphaCI[2], thetaCI[2]),
                                  parameter = c("power parameter", "effect size")))
    invisible(out)

}

## ## generic (retired) function
## ## function to compute specified quantile based on density function
## quantileFun <- function(densFun, q, support = c(-Inf, Inf),
##                         searchInt = c(support[1], support[2]),
##                         cubature = TRUE, ...) {
##     if (q == 0) {
##         res <- support[1]
##     } else if (q == 1) {
##         res <- support[2]
##     } else {
##         if (cubature == TRUE) {
##             cubFun <- function(x) matrix(densFun(x), ncol = 1)
##             rootFun <- function(x) {
##                 cubature::hcubature(f = cubFun, lowerLimit = support[1],
##                                     upperLimit = x,
##                                     vectorInterface = TRUE,
##                                     ... = ...)$integral - q
##             }
##         } else {
##             rootFun <- function(x) {
##                 stats::integrate(f = densFun, lower = support[1],
##                                  upper = x,
##                                  ... = ...)$value - q
##             }
##         }
##         res <- stats::uniroot(f = rootFun, interval = searchInt)$root
##     }
##     return(res)
## }


## ## function to compute highest posterior density interval
## hpdCI <- function(densFun, level = 0.95, support = c(-Inf, Inf),
##                   searchInt = c(support[1], support[2]), levelInt = c(0, 1),
##                   cubature = TRUE) {
##     ## determine the credible interval with the smallest width
##     optFun. <- function(alpha) {
##         width <- quantileFun(densFun = densFun, q = alpha + level,
##                              support = support, searchInt = searchInt,
##                              cubature = cubature) -
##             quantileFun(densFun = densFun, q = alpha, support = support,
##                         searchInt = searchInt, cubature = cubature)

##         return(width)
##     }
##     optFun <- Vectorize(FUN = optFun.)
##     minLower <- stats::optim(par = (1 - level)/2, fn = optFun,
##                              method = "L-BFGS-B", lower = levelInt[1],
##                              upper = levelInt[2] - level)$par
##     CI <- c(quantileFun(densFun = densFun, q = minLower, support = support,
##                         searchInt = searchInt, cubature = cubature),
##             quantileFun(densFun = densFun, q = minLower + level,
##                         support = support, searchInt = searchInt, cubature = cubature))
##     return(CI)
## }
