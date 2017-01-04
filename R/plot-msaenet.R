#' Plot msaenet Model Objects
#'
#' Plot msaenet model objects.
#'
#' @param x An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param ... Other parameters (not used).
#'
#' @method plot msaenet
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @importFrom lattice parallelplot
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
#' plot(msaenet.fit)

plot.msaenet = function(x, ...) {

  if (!.is.msaenet(x))
    stop('object class must be "msaenet"')

  if (.is.multistep(x))
    beta.mat = do.call(cbind, x$'beta.list')

  if (.is.adaptive(x))
    beta.mat = as.matrix(cbind(x$'beta.first', x$'beta'))

  colnames(beta.mat) = paste('Step', 1L:ncol(beta.mat))

  p = parallelplot(~ beta.mat, horizontal.axis = FALSE)
  print(p)

}
