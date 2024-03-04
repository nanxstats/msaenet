# msaenet 3.1.1

## Improvements

- Use a proper, three-component version number following Semantic Versioning.
- Fix warnings about single lambda (#11).
- Fix "lost braces" check notes on r-devel and check notes about `LazyData`.
- Fix code linting issues.
- Use GitHub Actions to build the pkgdown site.

# msaenet 3.1

## Improvements

- Added detailed signal-to-noise ratio (SNR) definition in `msaenet.sim.gaussian()`.
- Updated the example code in the vignette to make it work better with the most recent version of glmnet (2.0-16).
- Updated GitHub repository links due to the handle change.
- Updated the vignette style.

# msaenet 3.0

## New features

- Added a new argument `penalty.factor.init` to support customized penalty factor applied to each coefficient in the initial estimation step. This is useful for incorporating prior information about variable weights, for example, emphasizing specific clinical variables. We thank Xin Wang from University of Michigan for this feedback [[#4](https://github.com/nanxstats/msaenet/issues/4)].

# msaenet 2.9

## Improvements

- New URL for the documentation website: https://nanx.me/msaenet/.

# msaenet 2.8

## New features

- Added a Cleveland dot plot option `type = "dotplot"` in `plot.msaenet()`. This plot offers a direct visualization of the model coefficients at the optimal step.

# msaenet 2.7

## Bug fixes

- Fixed the missing arguments issue when `init = "ridge"`.

# msaenet 2.6

## Improvements

- Added two arguments `lower.limits` and `upper.limits` to support coefficient constraints in `aenet()` and `msaenet()` [[#1](https://github.com/nanxstats/msaenet/issues/1)].

# msaenet 2.5

## Improvements

- Better code indentation style.
- Update gallery images in `README.md`.

# msaenet 2.4

## Improvements

- Improved graphical details for coefficient path plots, following the general graphic style in the ESL (*The Elements of Statistical Learning*) book.
- More options available in `plot.msaenet()` for extra flexibility: it is now possible to set important properties of the label appearance such as position, offset, font size, and axis titles via the new arguments `label.pos`, `label.offset`, `label.cex`, `xlab`, and `ylab`.

# msaenet 2.3

## Improvements

- Reduced model saturation cases and improved speed at the initialization step for MCP-net and SCAD-net based models when `init = "ridge"`, by using the ridge estimation implementation from `glmnet`. As a benefit, we now have a more aligned baseline for the comparison between elastic-net based models and MCP-net/SCAD-net based models when `init = "ridge"`.
- Style improvements in code and examples: reduced whitespace with a new formatting scheme.

# msaenet 2.2

## New features

- Added BIC, EBIC, and AIC in addition to k-fold cross-validation for model selection.
- Added new arguments `tune` and `tune.nsteps` to controls this for selecting the optimal model for each step, and the optimal model among all steps (i.e. the optimal step).
- Added arguments `ebic.gamma` and `ebic.gamma.nsteps` to control the EBIC tuning parameter, if `ebic` is specified by `tune` or `tune.nsteps`.
- Redesigned plot function: now supports two types of plots (coefficient path, screeplot of the optimal step selection criterion), optimal step highlighting, variable labeling, and color palette customization. See `?plot.msaenet` for details.

## Improvements

- Renamed previous argument `gamma` (scaling factor for adaptive weights) to `scale` to avoid possible confusion.
- Reset the default values of candidate concavity parameter `gammas` to be 3.7 for SCAD-net and 3 for MCP-net.
- Unified the supported model `family` in all model types to be `"gaussian"`, `"binomial"`, `"poisson"`, and `"cox"`.

# msaenet 2.1

## New features

- Added functions `msaenet.sim.binomial()`, `msaenet.sim.poisson()`, `msaenet.sim.cox()` to generate simulation data for logistic, Poisson, and Cox regression models.
- Added function `msaenet.fn()` for computing the number of false negative selections in msaenet models.
- Added function `msaenet.mse()` for computing mean squared error (MSE).

## Improvements

- Speed improvements in `msaenet.sim.gaussian()` by more vectorization when generating correlation matrices.
- Added parameters `max.iter` and `epsilon` for MCP-net and SCAD-net related functions to have finer control over convergence criterion. By default, `max.iter = 10000` and `epsilon = 1e-4`.

# msaenet 2.0

## New features

- Added support for adaptive MCP-net. See `?amnet` for details.
- Added support for adaptive SCAD-net. See `?asnet` for details.
- Added support for multi-step adaptive MCP-net (MSAMNet). See `?msamnet` for details.
- Added support for multi-step adaptive SCAD-net (MSASNet). See `?msasnet` for details.
- Added `msaenet.nzv.all()` for displaying the indices of non-zero variables in all adaptive estimation steps.

## Improvements

- More flexible `predict.msaenet` method allowing users to specify prediction type.

# msaenet 1.1

## New features

- Added method `coef` for extracting model coefficients.
  See `?coef.msaenet` for details.

## Improvements

- New documentation website generated by pkgdown,
  with a full set of function documentation and vignettes available.
- Added Windows continuous integration support using AppVeyor.

# msaenet 1.0

## New features

- Initial version of the msaenet package
