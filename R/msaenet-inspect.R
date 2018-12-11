#' Get Indices of Non-Zero Variables
#'
#' Get the indices of non-zero variables from msaenet model objects.
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#'
#' @return Indices vector of non-zero variables in the model.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.nzv
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001
#' )
#' 
#' msaenet.fit <- msaenet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2),
#'   nsteps = 3L, seed = 1003
#' )
#' 
#' msaenet.nzv(msaenet.fit)
#' 
#' # coefficients of non-zero variables
#' coef(msaenet.fit)[msaenet.nzv(msaenet.fit)]
msaenet.nzv <- function(object) {
  if (!.is.msaenet(object)) {
    stop('object class must be "msaenet"')
  }

  idx <- which(abs(as.vector(object$"beta")) > .Machine$double.eps)
  idx
}

#' Get Indices of Non-Zero Variables in All Steps
#'
#' Get the indices of non-zero variables in all steps from msaenet model objects.
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#'
#' @return List containing indices vectors of non-zero variables in all steps.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.nzv.all
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001
#' )
#' 
#' msaenet.fit <- msaenet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2),
#'   nsteps = 3L, seed = 1003
#' )
#' 
#' msaenet.nzv.all(msaenet.fit)
msaenet.nzv.all <- function(object) {
  if (!.is.msaenet(object)) {
    stop('object class must be "msaenet"')
  }

  if (.is.multistep(object)) {
    n <- length(object$"beta.list")
    idx <- vector("list", n)
    for (i in 1L:n) idx[[i]] <-
        which(abs(as.vector(object[["beta.list"]][[i]])) > .Machine$double.eps)
  }

  if (.is.adaptive(object)) {
    idx <- vector("list", 2L)
    idx[[1L]] <- which(abs(as.vector(object[["beta.first"]])) > .Machine$double.eps)
    idx[[2L]] <- which(abs(as.vector(object[["beta"]])) > .Machine$double.eps)
  }

  idx
}

#' Extract Model Coefficients
#'
#' Extract model coefficients from the final model in msaenet model objects.
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param ... Additional parameters for \code{\link{coef}} (not used).
#'
#' @return A numerical vector of model coefficients.
#'
#' @method coef msaenet
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001
#' )
#' 
#' msaenet.fit <- msaenet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2),
#'   nsteps = 3L, seed = 1003
#' )
#' 
#' coef(msaenet.fit)
coef.msaenet <- function(object, ...) {
  if (!.is.msaenet(object)) {
    stop('object class must be "msaenet"')
  }

  bhat <- as.vector(object$"beta")
  bhat
}

#' Get the Number of False Positive Selections
#'
#' Get the number of false positive selections from msaenet model objects,
#' given the indices of true variables (if known).
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param true.idx Vector. Indices of true variables.
#'
#' @return Number of false positive variables in the model.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.fp
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001
#' )
#' 
#' msaenet.fit <- msaenet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2),
#'   nsteps = 3L, seed = 1003
#' )
#' 
#' msaenet.fp(msaenet.fit, 1:5)
msaenet.fp <- function(object, true.idx) {
  if (!.is.msaenet(object)) {
    stop('object class must be "msaenet"')
  }

  length(setdiff(msaenet.nzv(object), true.idx))
}

#' Get the Number of False Negative Selections
#'
#' Get the number of false negative selections from msaenet model objects,
#' given the indices of true variables (if known).
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param true.idx Vector. Indices of true variables.
#'
#' @return Number of false negative variables in the model.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.fn
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001
#' )
#' 
#' msaenet.fit <- msaenet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2),
#'   nsteps = 3L, seed = 1003
#' )
#' 
#' msaenet.fn(msaenet.fit, 1:5)
msaenet.fn <- function(object, true.idx) {
  if (!.is.msaenet(object)) {
    stop('object class must be "msaenet"')
  }

  length(setdiff(true.idx, msaenet.nzv(object)))
}

#' Get the Number of True Positive Selections
#'
#' Get the number of true positive selections from msaenet model objects,
#' given the indices of true variables (if known).
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param true.idx Vector. Indices of true variables.
#'
#' @return Number of true positive variables in the model.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.tp
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001
#' )
#' 
#' msaenet.fit <- msaenet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2),
#'   nsteps = 3L, seed = 1003
#' )
#' 
#' msaenet.tp(msaenet.fit, 1:5)
msaenet.tp <- function(object, true.idx) {
  if (!.is.msaenet(object)) {
    stop('object class must be "msaenet"')
  }

  length(intersect(msaenet.nzv(object), true.idx))
}
