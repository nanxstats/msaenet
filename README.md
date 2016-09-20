# msaenet

[![Build Status](https://travis-ci.org/road2stat/msaenet.svg?branch=master)](https://travis-ci.org/road2stat/msaenet)
[![CRAN Version](http://www.r-pkg.org/badges/version/msaenet)](http://www.r-pkg.org/badges/version/msaenet)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/msaenet)](http://cranlogs.r-pkg.org/badges/msaenet)

`msaenet` implements the multi-step adaptive elastic-net (MSAENet) algorithm for feature selection in high-dimensional regressions.

## Installation

To download and install `msaenet` from CRAN:

```r
install.packages("msaenet")
```

Or try the development version on GitHub:

```r
# install.packages("devtools")
devtools::install_github("road2stat/msaenet")
```

To load the package in R, simply use

```r
library("msaenet")
```

and you are all set. See [the vignette](http://msaenet.com/doc/) (can also be opened with `vignette("msaenet")` in R) for a quick-start guide.

## Paper Citation

Formatted citation:

Nan Xiao and Qing-Song Xu. (2015). Multi-step adaptive elastic-net: reducing false positives in high-dimensional variable selection. Journal of Statistical Computation and Simulation 85(18), 3755-3765.

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
