
<!-- README.md is generated from README.Rmd. Please edit that file -->

# msaenet <img src="man/figures/logo.png" align="right" width="120" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/nanxstats/msaenet/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nanxstats/msaenet/actions/workflows/R-CMD-check.yaml)
[![CRAN
Version](https://www.r-pkg.org/badges/version/msaenet)](https://cran.r-project.org/package=msaenet)
[![Downloads from the RStudio CRAN
mirror](https://cranlogs.r-pkg.org/badges/msaenet)](https://cran.r-project.org/package=msaenet)
<!-- badges: end -->

`msaenet` implements the multi-step adaptive elastic-net (MSAENet)
algorithm for feature selection in high-dimensional regressions proposed
in Xiao and Xu (2015)
\<[DOI:10.1080/00949655.2015.1016944](https://www.tandfonline.com/doi/full/10.1080/00949655.2015.1016944)\>
([PDF](https://nanx.me/papers/msaenet.pdf)).

Nonconvex multi-step adaptive estimations based on MCP-net or SCAD-net
are also supported.

Check `vignette("msaenet")` for a quick-start.

## Installation

You can install `msaenet` from CRAN:

``` r
install.packages("msaenet")
```

Or try the development version on GitHub:

``` r
remotes::install_github("nanxstats/msaenet")
```

## Citation

> Nan Xiao and Qing-Song Xu. (2015). Multi-step adaptive elastic-net:
> reducing false positives in high-dimensional variable selection.
> *Journal of Statistical Computation and Simulation* 85(18), 3755–3765.

## Gallery

### Adaptive Elastic-Net / Multi-Step Adaptive Elastic-Net

<img src="man/figures/README-msaenet-1.png" width="100%" />

### Adaptive MCP-Net / Multi-Step Adaptive MCP-Net

<img src="man/figures/README-msamnet-1.png" width="100%" />

### Adaptive SCAD-Net / Multi-Step Adaptive SCAD-Net

<img src="man/figures/README-msasnet-1.png" width="100%" />

## Contribute

To contribute to this project, please take a look at the [Contributing
Guidelines](CONTRIBUTING.md) first. Please note that this project is
released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.

## License

msaenet is free and open source software, licensed under GPL-3.
