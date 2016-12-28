#' Make Predictions from an AENet/MSAENet Model
#'
#' Make predictions on new data by an AENet/MSAENet model object.
#'
#' @param object An object of class \code{msaenet} produced
#' by \code{\link{aenet}} or \code{\link{msaenet}}.
#' @param newx New data to predict with.
#' @param ... Additional parameters for \code{\link{predict}}.
#'
#' @return Numeric matrix of the predicted values.
#'
#' @method predict msaenet
#'
#' @author Nan Xiao <\url{http://nanx.me}>
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

  if (missing(newx)) stop('Please specify newx')

  if (!('msaenet' %in% class(object)))
    stop('object class must be "msaenet"')

  pred = predict(object$'model', newx = newx)
  return(pred)

}
