# ppRep

**ppRep** is an R package for Bayesian estimation and testing of effect sizes
based on original and replication study using a power prior framework for
dynamic discounting of the original data. For more information, see the preprint

Pawel, S., Aust, F., Held, L., and Wagenmakers, E.-J. (2022). Power Priors for
  Replication Studies.
  doi:[10.48550/arXiv.2207.14720](https://doi.org/10.48550/arXiv.2207.14720)

## Installation

```r
remotes::install_github(repo = "SamCH93/ppRep")
```

## Usage

``` r
library("ppRep")

## data from Protzko et al. (2020)
to <- 0.2 # original effect estimate
so <- 0.05 # original standard error
tr <- 0.1 # replication effect estimate
sr <- 0.05 # replication standard error

## compute and plot posterior density
plotPP(tr = tr, sr = sr, to = to, so = so)
```
![Plot of joint posterior distribution of power parameter and effect size.](posterior.png)

