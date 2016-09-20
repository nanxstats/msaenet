#' Get Indices of Non-Zero Variables
#'
#' Get the indices of non-zero variables from AENet/MSAENet model objects.
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}} or \code{\link{msaenet}}.
#'
#' @return Indices vector of non-zero variables in the model.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export msaenet.nzv
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' msaenet.fit = msaenet(dat$x.tr, dat$y.tr,
#'                       alphas = seq(0.2, 0.8, 0.2),
#'                       nsteps = 3L, seed = 1003)
#'
#' msaenet.nzv(msaenet.fit)

msaenet.nzv = function(object) {

  if (!('msaenet' %in% class(object)))
    stop('object class must be "msaenet"')

  return(which(as.vector(object$'beta') != 0))

}

#' Get the Number of False Positive Selections
#'
#' Get the number of false positive selections from AENet/MSAENet model objects,
#' given the indices of true variables (if known).
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}} or \code{\link{msaenet}}.
#' @param true.idx Vector. Indices of true variables.
#'
#' @return Number of false positive variables in the model.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export msaenet.fp
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' msaenet.fit = msaenet(dat$x.tr, dat$y.tr,
#'                       alphas = seq(0.2, 0.8, 0.2),
#'                       nsteps = 3L, seed = 1003)
#'
#' msaenet.fp(msaenet.fit, 1:5)

msaenet.fp = function(object, true.idx) {

  if (!('msaenet' %in% class(object)))
    stop('object class must be "msaenet"')

  return(length(setdiff(msaenet.nzv(object), true.idx)))

}

#' Get the Number of True Positive Selections
#'
#' Get the number of true positive selections from AENet/MSAENet model objects,
#' given the indices of true variables (if known).
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}} or \code{\link{msaenet}}.
#' @param true.idx Vector. Indices of true variables.
#'
#' @return Number of true positive variables in the model.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export msaenet.tp
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' msaenet.fit = msaenet(dat$x.tr, dat$y.tr,
#'                       alphas = seq(0.2, 0.8, 0.2),
#'                       nsteps = 3L, seed = 1003)
#'
#' msaenet.tp(msaenet.fit, 1:5)

msaenet.tp = function(object, true.idx) {

  if (!('msaenet' %in% class(object)))
    stop('object class must be "msaenet"')

  return(length(intersect(msaenet.nzv(object), true.idx)))

}
