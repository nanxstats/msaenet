#' Make Predictions from an msaenet Model
#'
#' Make predictions on new data by a msaenet model object.
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param newx New data to predict with.
#' @param ... Additional parameters, particularly prediction \code{type} in
#' \code{\link[glmnet]{predict.glmnet}}, \code{\link[ncvreg]{predict.ncvreg}},
#' or \code{\link[ncvreg]{predict.ncvsurv}}.
#'
#' @return Numeric matrix of the predicted values.
#'
#' @method predict msaenet
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom stats predict
#'
#' @export
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
#' msaenet.pred = predict(msaenet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, msaenet.pred)

predict.msaenet = function(object, newx, ...) {

  if (missing(newx)) stop('Please specify newx to predict on')

  if (!.is.msaenet(object))
    stop('object class must be "msaenet"')

  if (.is.glmnet(object$'model'))
    pred = predict(object$'model', newx = newx, ...)

  if (.is.ncvreg(object$'model'))
    pred = predict(object$'model', X = newx, ...)

  pred

}
