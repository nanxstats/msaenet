# msaenet  <a href="http://msaenet.com"><img src="http://nanx.me/images/project-msaenet.png" align="right" alt="logo" height="180" width="180" /></a>

[![Build Status](https://travis-ci.org/road2stat/msaenet.svg?branch=master)](https://travis-ci.org/road2stat/msaenet)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/road2stat/msaenet?branch=master&svg=true)](https://ci.appveyor.com/project/road2stat/msaenet)
[![CRAN Version](http://www.r-pkg.org/badges/version/msaenet)](https://cran.r-project.org/package=msaenet)
[![Downloads from the RStudio CRAN mirror](http://cranlogs.r-pkg.org/badges/msaenet)](https://cran.r-project.org/package=msaenet)

`msaenet` implements the multi-step adaptive elastic-net (MSAENet) algorithm for feature selection in high-dimensional regressions proposed in Xiao and Xu (2015) <[DOI:10.1080/00949655.2015.1016944](http://www.tandfonline.com/doi/full/10.1080/00949655.2015.1016944)> ([PDF](https://drive.google.com/file/d/0B1YdO4YnMkAxeFUtZ3FLY1dLN2s/view)).

Multi-step adaptive estimation based on MCP-net or SCAD-net is also supported.

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

<img src="http://msaenet.com/img/aenet.png" width="49%" alt="aenet">
<img src="http://msaenet.com/img/msaenet.png" width="49%" alt="msaenet">

### Adaptive MCP-Net / Multi-Step Adaptive MCP-Net

<img src="http://msaenet.com/img/amnet.png" width="49%" alt="amnet">
<img src="http://msaenet.com/img/msamnet.png" width="49%" alt="msamnet">

### Adaptive SCAD-Net / Multi-Step Adaptive SCAD-Net

<img src="http://msaenet.com/img/asnet.png" width="49%" alt="asnet">
<img src="http://msaenet.com/img/msasnet.png" width="49%" alt="msasnet">

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

See [the vignette](http://msaenet.com/doc/) (can be opened with `vignette("msaenet")` in R) for a quick-start guide.
