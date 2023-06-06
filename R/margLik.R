#' @title Marginal likelihood of replication effect estimate
#'
#' @description This function computes the marginal likelihood of the
#'     replication effect estimate \code{tr} under the power prior model
#'     \deqn{f(\code{tr}|\code{to}, \code{so}, \code{sr}, \code{x}, \code{y}) =
#'     \int_0^1 \int_{-\infty}^{\infty} \mathrm{N}(\code{tr}; \theta,
#'     \code{sr}^2) \times \mathrm{N}(\theta; \mu, \phi)
#'     \times \mbox{Beta}(\alpha; \code{x}, \code{y}) ~\mbox{d}\theta~
#'     \mbox{d}\alpha}{int int N(tr;theta,sr^2) N(theta; mu, phi)
#'     Beta(alpha; x, y) dtheta dalpha} with \eqn{\phi = 1/(1/\code{v} +
#'     \alpha/\code{so}^2)}{phi = 1/(1/v + alpha/so^2)} and \eqn{\mu =
#'     \phi\{(\alpha\times\code{to})/\code{so}^2 + \code{m}/\code{v}\}}{mu =
#'     phi*(alpha*to/so^2 + m/v)} using numerical integration.
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
#' @param ... Additional arguments passed to \code{stats::integrate}.
#'
#' @return Marginal likelihood
#'
#' @author Samuel Pawel
#' @export
margLik <- function(tr, to, sr, so, x = 1, y = 1, m = 0, v = Inf, ...) {
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
        0 < v
    )

    ## integrate out alpha with numerical integration
    intFun <- function(alpha) {
        ## mean and variance of normalized power prior conditional on alpha
        pvar <- postNormVar(vardata = so^2, priorvar = v, alpha = alpha)
        pmean <- postNormMean(dat = to, vardata = so^2, priormean = m,
                              priorvar = v, alpha = alpha)

        ## effect size theta can be integrated analytically
        stats::dnorm(x = tr, mean = pmean, sd = sqrt(sr^2 + pvar)) *
            stats::dbeta(x = alpha, shape1 = x, shape2 = y)
    }
    res <- try(stats::integrate(f = intFun, lower = 0, upper = 1,
                                ... = ...)$value)
    if (inherits(res, "try-error")) {
        warnString <- paste("Numerical problems computing normalizing constant.",
                            "Try adjusting integration options \nwith ... argument.",
                            "See ?stats::integrate for available options.")
        warning(warnString)
        const <- NaN
    } else {
        const <- res
    }
    return(res)
}
