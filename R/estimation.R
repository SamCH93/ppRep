#' @title Variance of normalized power prior
#' @description This function computes the variance of a normalized power prior
#'     conditional on a fixed power parameter and an initial normal prior for
#'     the effect size.
#' @param vardata Variance of the data.
#' @param priorvar Variance parameter of initial normal prior. Defaults to Inf
#'     (uniform prior).
#' @param alpha Power parameter. Indicates to which power the likelihood of the
#'     data is raised. Can be set to a number in [0, 1]. Defaults to 1.
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
#' @param priormean Mean parameter of initial normal prior. Defaults to 0.
#' @param priorvar Variance parameter of initial normal prior. Defaults to Inf
#'     (uniform prior).
#' @param alpha Power parameter. Indicates to which power the likelihood of the
#'     data is raised. Can be set to a number in [0, 1]. Defaults to 1.
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
#'     Defaults to 1.
#' @param y Number of failures parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to 1.
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to 0.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to Inf (uniform prior).
#' @param ... Additional arguments for integration function.
#'
#' @return Posterior density
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{postPPalpha}}, \code{\link{postPPtheta}}
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


#' @title Marginal posterior density of power parameter
#'
#' @description This function computes the marginal posterior density of the
#'     power parameter \eqn{\alpha}{alpha}. A power prior for
#'     \eqn{\theta}{theta} is constructed by updating an initial normal prior
#'     \eqn{\theta \sim \mathrm{N}(\code{m}, \code{v})}{theta ~ N(m, v)} with
#'     the likelihood of the original data raised to the power of
#'     \eqn{\alpha}{alpha}. A marginal beta prior \eqn{\alpha \sim
#'     \mbox{Beta}(\code{x},\code{y})}{alpha ~ Beta(x, y)} is assumed.
#'
#' @param alpha Power parameter. Can be a vector.
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter of beta prior \eqn{\alpha}{alpha}.
#'     Defaults to 1.
#' @param y Number of failures parameter of beta prior \eqn{\alpha}{alpha}.
#'     Defaults to 1.
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to 0.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to Inf (uniform prior).
#' @param ... Additional arguments for integration function.
#'
#' @return Marginal posterior density of power parameter
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{postPP}}, \code{\link{postPPtheta}}
#'
#' @examples
#' alpha <- seq(0, 1, 0.001)
#' margpostdens <- postPPalpha(alpha = alpha, tr = 0.1, to = 0.2, sr = 0.05, so = 0.05)
#' plot(alpha, margpostdens, type = "l", xlab = bquote("Power paramter" ~ alpha),
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

#' @title Marginal posterior density of effect size
#'
#' @description This function computes the marginal posterior density of the
#'     effect size \eqn{\theta}{theta}. A power prior for \eqn{\theta}{theta} is
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
#' @param x Number of successes parameter for beta prior of power
#'     parameter \eqn{\alpha}{alpha}. Defaults to 1. Is only taken into account
#'     when \code{alpha = NA}.
#' @param y Number of failures parameter for beta prior of power
#'     parameter \eqn{\alpha}{alpha}. Defaults to 1. Is only taken into account
#'     when \code{alpha = NA}.
#' @param alpha Power parameter. Can be set to a number between 0 and 1.
#'     Defaults to \code{NA}.
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to 0.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to Inf (uniform prior).
#' @param ... Additional arguments for integration function.
#'
#' @return Marginal posterior density of effect size
#'
#' @author Samuel Pawel
#'
#' @seealso \code{\link{postPP}}, \code{\link{postPPalpha}}
#'
#' @examples
#' theta <- seq(0, 0.6, 0.001)
#' margpostdens <- postPPtheta(theta = theta, tr = 0.1, to = 0.2, sr = 0.05, so = 0.05)
#' plot(theta, margpostdens, type = "l", xlab = bquote("Effect size" ~ theta),
#'      ylab = "Marginal posterior density", las = 1)
#' @export
postPPtheta <- function(theta, tr, sr, to, so, x = 1, y = 1, alpha = NA, m = 0,
                        v = Inf, ...) {
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
        0 < v
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
        if (!is.finite(v)) {
            margdens <- stats::dnorm(x = tr, mean = theta, sd = sr) *
                abs(hypergeo::genhypergeo(U = x + 0.5, L = x + y + 0.5,
                                          z = -(to - theta)^2/(2*so^2))) *
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


#' @title Joint and marginal posterior density plots
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
#'     Defaults to 1.
#' @param y Number of failures parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to 1.
#' @param m Mean parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to 0.
#' @param v Variance parameter of initial normal prior for \eqn{\theta}{theta}.
#'     Defaults to Inf (uniform prior).
#' @param thetaRange Range of effect sizes. Defaults to three standard errors
#'     around the replication effect estimate.
#' @param alphaRange Range of power parameters. Defaults to the range between
#'     zero and one.
#' @param nGrid Number of grid points. Defaults to 100.
#' @param plot Logical indicating whether data should be plotted. If
#'     \code{FALSE} only the data used for plotting are returned.
#' @param ... Additional arguments passed to \code{plot} function.
#'
#' @return Plots joint and marginal posterior densities, invisibly returns the
#'     data for the plots.
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
                   alphaRange = c(0, 1), nGrid = 100, plot = TRUE, ...) {
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
        !is.na(plot)
    )

    ## grids for computing the densities
    alphaGrid <- seq(alphaRange[1], alphaRange[2], length.out = nGrid)
    thetaGrid <- seq(thetaRange[1], thetaRange[2], length.out = nGrid)
    jointGrid <- expand.grid(alpha = alphaGrid, theta = thetaGrid)

    ## compute densities
    jointdens <- postPP(theta = jointGrid$theta, alpha = jointGrid$alpha, tr = tr,
                        sr = sr, to = to, so = so, x = x, y = y, m = m, v = v)
    alphadens <- postPPalpha(alpha = alphaGrid, tr = tr, sr = sr, to = to,
                             so = so, x = x, y = y, m = m, v = v)
    thetadens <- postPPtheta(theta = thetaGrid, tr = tr, sr = sr, to = to,
                             so = so, x = x, y = y, m = m, v = v)

    if (plot) {
        oldpar <- graphics::par("mfrow")
        graphics::layout(mat = matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE))
        ## joint density
        jointdensMat <- matrix(data = jointdens, ncol = nGrid, byrow = TRUE)
        graphics::image(x = thetaGrid, y = alphaGrid, z = jointdensMat,
                        xlab = bquote("Effect size" ~ theta),
                        ylab = bquote("Power parameter" ~ alpha),
                        main = bquote(plain("Joint posterior density")),
                        col = grDevices::hcl.colors(n = 100, palette = "Blues 3",
                                                    rev = TRUE),
                        las = 1)
        graphics::contour(x = thetaGrid, y = alphaGrid, z = jointdensMat, add = TRUE,
                          drawlabels = FALSE, nlevels = 5, col = "#00000080")
        ## power parameter
        plot(alphaGrid, alphadens, xlab = bquote("Power parameter" ~ alpha),
             ylab = "Marginal posterior density", type = "l", las = 1)
        ## effect size
        plot(thetaGrid, thetadens, xlab = bquote("Effect size" ~ theta),
             ylab = "Marginal posterior density", type = "l", las = 1)
        graphics::par(mfrow = oldpar)
    }

    ## return plot data
    out <- list(jointDF = data.frame(jointGrid, density = jointdens),
                alphaDF = data.frame(alpha = alphaGrid,
                                     density = alphadens),
                thetaDF = data.frame(theta = thetaGrid,
                                     density = thetadens))
    invisible(out)

}
