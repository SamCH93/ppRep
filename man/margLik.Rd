% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/margLik.R
\name{margLik}
\alias{margLik}
\title{Marginal likelihood of replication effect estimate}
\usage{
margLik(tr, to, sr, so, x = 1, y = 1, m = 0, v = Inf, ...)
}
\arguments{
\item{tr}{Effect estimate of the replication study.}

\item{to}{Effect estimate of the original study.}

\item{sr}{Standard error of the replication effect estimate.}

\item{so}{Standard error of the original effect estimate.}

\item{x}{Number of successes parameter of beta prior for \eqn{\alpha}{alpha}.
Defaults to \code{1}.}

\item{y}{Number of failures parameter of beta prior for \eqn{\alpha}{alpha}.
Defaults to \code{1}.}

\item{m}{Mean parameter of initial normal prior for \eqn{\theta}{theta}.
Defaults to \code{0}.}

\item{v}{Variance parameter of initial normal prior for \eqn{\theta}{theta}.
Defaults to \code{Inf} (uniform prior).}

\item{...}{Additional arguments passed to \code{stats::integrate}.}
}
\value{
Marginal likelihood
}
\description{
This function computes the marginal likelihood of the
    replication effect estimate \code{tr} under the power prior model
    \deqn{f(\code{tr}|\code{to}, \code{so}, \code{sr}, \code{x}, \code{y}) =
    \int_0^1 \int_{-\infty}^{\infty} \mathrm{N}(\code{tr}; \theta,
    \code{sr}^2) \times \mathrm{N}(\theta; \mu, \phi)
    \times \mbox{Beta}(\alpha; \code{x}, \code{y}) ~\mbox{d}\theta~
    \mbox{d}\alpha}{int int N(tr;theta,sr^2) N(theta; mu, phi)
    Beta(alpha; x, y) dtheta dalpha} with \eqn{\phi = 1/(1/\code{v} +
    \alpha/\code{so}^2)}{phi = 1/(1/v + alpha/so^2)} and \eqn{\mu =
    \phi\{(\alpha\times\code{to})/\code{so}^2 + \code{m}/\code{v}\}}{mu =
    phi*(alpha*to/so^2 + m/v)} using numerical integration.
}
\author{
Samuel Pawel
}
