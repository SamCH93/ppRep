# ppRep

**ppRep** is an R package for Bayesian estimation and testing of effect sizes
based on original and replication study using the power prior framework.

## Installation

```r
remotes::install_github(repo = "SamCH93/ciCalibrate")
```

## Usage

``` r
library("ppRep")

## data from Protzko et al. (2020)
to <- 0.2 # original effect estimate
so <- 0.05 # original standard error
tr <- 0.1 # replication effect estimate
sr <- 0.05 # replication standard error

## grid of power parameter and effect size for joint posterior
alpha <- seq(0, 1, length.out = 200)
theta <- seq(0, 0.3, length.out = 200)
parGrid <- expand.grid(alpha = alpha, theta = theta)

## compute and plot posterior density
postdens <- postPP(theta = parGrid$theta, alpha = parGrid$alpha, tr = tr,
                   sr = sr, to = to, so = so)
postdensMat <- matrix(data = postdens, ncol = 200, byrow = TRUE)
filled.contour(x = theta, y = alpha, z = postdensMat,
               xlab = bquote("Effect size" ~ theta),
               ylab = bquote("Power parameter" ~ alpha), nlevels = 15,
               color.palette = function(n) hcl.colors(n = n, palette = "viridis"))

```
![Plot of joint posterior distribution of power parameter and effect size.](posterior.png)

## References

* Pawel, S., Aust, F., Held, L., and Wagenmakers, E.-J. (2022). Power Priors for
  Replication Studies. To appear on arXiv soon.
  [10.48550/arXiv.XXXX.XXXXX](https://doi.org/10.48550/arXiv.XXXX.XXXXX)

