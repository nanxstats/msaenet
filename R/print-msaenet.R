#' Print msaenet Model Information
#'
#' Print msaenet model objects (currently, only
#' printing the model information of the final step).
#'
#' @param x An object of class \code{msaenet}.
#' @param ... Additional parameters for \code{\link{print}} (not used).
#'
#' @method print msaenet
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @importFrom utils capture.output
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
#' print(msaenet.fit)

print.msaenet = function(x, ...) {

  cat(paste('Call:', paste(capture.output(x$'call'),
                           collapse = '\n')), '\n')

  if (.is.ncvreg(x$'model')) {
    model.info = data.frame(.ncvdf(x$'model'),
                            x$'model'$'lambda',
                            x$'model'$'gamma',
                            x$'model'$'alpha')
    names(model.info) = c('Df', 'Lambda', 'Gamma', 'Alpha')
    print(model.info)
  }

  if (.is.glmnet(x$'model')) {
    model.info = data.frame(x$'model'$'df',
                            x$'model'$'dev.ratio',
                            x$'model'$'lambda')
    names(model.info) = c('Df', '%Dev', 'Lambda')
    print(model.info)
  }

}
