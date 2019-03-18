# msaenet  <a href="https://nanx.me/msaenet/"><img src="https://i.imgur.com/IUrYMK6.png" align="right" alt="logo" height="180" width="180" /></a>

[![Build Status](https://travis-ci.org/nanxstats/msaenet.svg?branch=master)](https://travis-ci.org/nanxstats/msaenet)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/nanxstats/msaenet?branch=master&svg=true)](https://ci.appveyor.com/project/nanxstats/msaenet)
[![CRAN Version](https://www.r-pkg.org/badges/version/msaenet)](https://cran.r-project.org/package=msaenet)
[![Downloads from the RStudio CRAN mirror](https://cranlogs.r-pkg.org/badges/msaenet)](https://cran.r-project.org/package=msaenet)

`msaenet` implements the multi-step adaptive elastic-net (MSAENet) algorithm for feature selection in high-dimensional regressions proposed in Xiao and Xu (2015) <[DOI:10.1080/00949655.2015.1016944](https://www.tandfonline.com/doi/full/10.1080/00949655.2015.1016944)> ([PDF](https://nanx.me/papers/msaenet.pdf)).

Nonconvex multi-step adaptive estimations based on MCP-net or SCAD-net are also supported.

## Paper Citation

Formatted citation:

Nan Xiao and Qing-Song Xu. (2015). Multi-step adaptive elastic-net: reducing false positives in high-dimensional variable selection. _Journal of Statistical Computation and Simulation_ 85(18), 3755-3765.

BibTeX entry:

```
@article{xiao2015msaenet,
  title={Multi-step adaptive elastic-net: reducing false positives in high-dimensional variable selection},
  author={Xiao, Nan and Xu, Qing-Song},
  journal={Journal of Statistical Computation and Simulation},
  volume={85},
  number={18},
  pages={3755--3765},
  year={2015},
  publisher={Taylor \& Francis}
}
```

## Gallery

### Adaptive Elastic-Net / Multi-Step Adaptive Elastic-Net

<img src="https://nanx.me/msaenet/img/enet.png" width="100%" alt="enet">

### Adaptive MCP-Net / Multi-Step Adaptive MCP-Net

<img src="https://nanx.me/msaenet/img/mnet.png" width="100%" alt="mnet">

### Adaptive SCAD-Net / Multi-Step Adaptive SCAD-Net

<img src="https://nanx.me/msaenet/img/snet.png" width="100%" alt="snet">

## Installation

To download and install `msaenet` from CRAN:

```r
install.packages("msaenet")
```

Or try the development version on GitHub:

```r
# install.packages("devtools")
devtools::install_github("nanxstats/msaenet")
```

[Browse the vignette](https://nanx.me/msaenet/articles/msaenet.html) (can be opened with `vignette("msaenet")` in R) for a quick-start.

## Contribute

To contribute to this project, please take a look at the [Contributing Guidelines](CONTRIBUTING.md) first. Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

## License

msaenet is free and open source software, licensed under GPL-3.
