#' @title Marginal likelihood of replication effect estimate under power prior model
#'
#' @description This function computes the marginal likelihood of the
#'     replication effect estimate \code{tr} under the power prior model
#'     \deqn{f(\code{tr}|\code{to}, \code{so}, \code{sr}, \code{x}, \code{y}) =
#'     \int_0^1 \int_{-\infty}^{\infty} \mathrm{N}(\code{tr}; \theta,
#'     \code{sr}^2) \times \mathrm{N}(\theta; \code{to}, \code{so}^2/\alpha)
#'     \times \mbox{Beta}(\alpha; \code{x}, \code{y}) ~\mbox{d}\theta~
#'     \mbox{d}\alpha}{int int N(tr;theta,sr^2) N(theta; to, so^2/alpha)
#'     Beta(alpha; x, y) dtheta dalpha} using numerical integration.
#'
#' @param tr Effect estimate of the replication study.
#' @param to Effect estimate of the original study.
#' @param sr Standard error of the replication effect estimate.
#' @param so Standard error of the replication effect estimate.
#' @param x Number of successes parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to 1.
#' @param y Number of failures parameter of beta prior for \eqn{\alpha}{alpha}.
#'     Defaults to 1.
#' @param ... Additional arguments for integration function
#'
#' @return Marginal likelihood
#'
#' @author Samuel Pawel
#' @export
margLik <- function(tr, to, sr, so, x = 1, y = 1, ...) {
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
        0 <= y
    )

    ## integrate out alpha with numerical integration
    intFun <- function(alpha) {
        ## effect size theta can be integrated analytically
        stats::dnorm(x = tr, mean = to, sd = sqrt(sr^2 + so^2/alpha)) *
            stats::dbeta(x = alpha, shape1 = x, shape2 = y)
    }
    res <- try(stats::integrate(f = intFun, lower = 0, upper = 1,
                                ... = ...)$value)
    if (class(res) == "try-error") {
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
