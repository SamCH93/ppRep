#' @title Posterior density of effect size and power parameter
#'
#' @description This function computes the posterior density of effect size
#'     \eqn{\theta}{theta} and power parameter \eqn{\alpha}{alpha} assuming a
#'     normal likelihood for original and replication effect estimate. A power
#'     prior for \eqn{\theta}{theta} is constructed by updating an initial flat
#'     prior with the likelihood of the original data raised to the power of
#'     \eqn{\alpha}{alpha}. A marginal beta prior \eqn{\alpha \sim
#'     \mbox{Beta}(x, y)}{alpha ~ Beta(x, y)} is assumed.
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
postPP <- function(theta, alpha, tr, sr, to, so, x = 1, y = 1, ...) {
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
        0 <= y
    )

    ## compute normalizing constant just once
    nC <- margLik(tr = tr, sr = sr, to = to, so = so, x = x, y = y,
                       ... = ...)
    if (is.nan(nC)) {
        out <- rep(x = NaN, times = pmax(length(alpha), length(theta)))
        return(out)
    }

    ## compute posterior density
    densProp <- stats::dnorm(x = tr, mean = theta, sd = sr) *
        stats::dnorm(x = theta, mean = to, sd = so/sqrt(alpha)) *
        stats::dbeta(x = alpha, shape1 = x, shape2 = y)
    dens <- densProp / nC
    return(dens)
}


#' @title Marginal posterior density of power parameter
#'
#' @description This function computes the marginal posterior density of the
#'     power parameter \eqn{\alpha}{alpha}. A power prior for
#'     \eqn{\theta}{theta} is constructed by updating an initial flat prior with
#'     the likelihood of the original data raised to the power of
#'     \eqn{\alpha}{alpha}. A marginal beta prior \eqn{\alpha \sim
#'     \mbox{Beta}(x, y)}{alpha ~ Beta(x, y)} is assumed.
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
postPPalpha <- function(alpha, tr, sr, to, so, x = 1, y = 1, ...) {
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
        0 <= y
    )

    ## compute normalizing constant just once
    nC <- margLik(tr = tr, sr = sr, to = to, so = so, x = x, y = y,
                       ... = ...)
    if (is.nan(nC)) {
        out <- rep(x = NaN, times = length(alpha))
        return(out)
    }

    ## compute marginal posterior density
    margdensProp <- stats::dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2/alpha)) *
        stats::dbeta(x = alpha, shape1 = x, shape2 = y)
    margdens <- margdensProp / nC
    return(margdens)
}

#' @title Marginal posterior density of effect size
#'
#' @description This function computes the marginal posterior density of the
#'     effect size \eqn{\theta}{theta}. A power prior for \eqn{\theta}{theta} is
#'     constructed by updating an initial flat prior with likelihood of the
#'     original data raised to the power of \eqn{\alpha}{alpha}. The power
#'     parameter \eqn{\alpha}{alpha} can either be fixed to some value between 0
#'     and 1, or it can have a Beta distribution \eqn{\alpha \sim
#'     \mbox{Beta}(\code{x}, \code{y})}{alpha ~ Beta(\code{x}, \code{y})}.
#'
#' @param theta Effect size. Can be a vector.
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter \eqn{x}{x} for beta prior of power
#'     parameter \eqn{\alpha}{alpha}. Defaults to 1. Is only taken into account
#'     when \code{alpha = NA}.
#' @param y Number of failures parameter \eqn{y}{y} for beta prior of power
#'     parameter \eqn{\alpha}{alpha}. Defaults to 1. Is only taken into account
#'     when \code{alpha = NA}.
#' @param alpha Power parameter. Can be set to a number between 0 and 1.
#'     Defaults to \code{NA}.
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
postPPtheta <- function(theta, tr, sr, to, so, x = 1, y = 1, alpha = NA, ...) {
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
          (alpha <= 1)))
    )

    ## alpha fixed
    if (!is.na(alpha)) {
        postVar <- 1/(1/sr^2 + alpha/so^2)
        postMean <- (tr/sr^2 + to*alpha/so^2)*postVar
        margdens <- stats::dnorm(x = theta, mean = postMean, sd = sqrt(postVar))
    } else { ## alpha random

        ## compute normalizing constant just once
        nC <- margLik(tr = tr, sr = sr, to = to, so = so, x = x, y = y,
                      ... = ...)
        if (is.nan(nC)) {
            out <- rep(x = NaN, times = length(theta))
            return(out)
        }

        ## compute marginal posterior
        margdens <- vapply(X = theta, FUN = function(thetai) {
            ## integrate out power parameter alpha
            intFun <- function(alpha) {
                stats::dnorm(x = thetai, mean = to, sd = so/sqrt(alpha)) *
                    stats::dbeta(x = alpha, shape1 = x, shape2 = y)
            }
            int <- try(stats::integrate(f = intFun, lower = 0, upper = 1,
                                        ... = ...)$value)
            if (class(int) == "try-error") {
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
    return(margdens)
}
