#' Plot msaenet Model Objects
#'
#' Plot msaenet model objects.
#'
#' @param x An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param type Plot type, \code{"coef"} for a coefficient path plot
#' across all estimation steps; \code{"criterion"} for a scree plot of
#' the model evaluation criterion used (CV error, AIC, BIC, or EBIC).
#' @param nsteps Maximum number of estimation steps to plot.
#' Default is to plot all steps.
#' @param highlight Should we highlight the "optimal" step
#' according to the criterion? Default is \code{TRUE}.
#' @param col Color palette to use for the coefficient paths.
#' If it is \code{NULL}, a default color palette will be assigned.
#' @param label Should we label all the non-zero variables of the
#' optimal step in the coefficient plot? Default is \code{FALSE}.
#' If \code{TRUE} and \code{label.vars = NULL}, the index of the
#' variables will be used as labels.
#' @param label.vars Labels to use for all the variables
#' if \code{label = "TRUE"}.
#' @param ... Other parameters (not used).
#'
#' @method plot msaenet
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom graphics axis lines matlines matplot plot text
#'
#' @export
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' msasnet.fit = msasnet(dat$x.tr, dat$y.tr,
#'                       alphas = seq(0.2, 0.8, 0.2),
#'                       nsteps = 5L, tune.nsteps = "ebic",
#'                       seed = 1003)
#'
#' plot(msasnet.fit)
#' plot(msasnet.fit, label = TRUE)
#' plot(msasnet.fit, label = TRUE, nsteps = 5)
#' plot(msasnet.fit, type = "criterion")
#' plot(msasnet.fit, type = "criterion", nsteps = 5)

plot.msaenet = function(x, type = c('coef', 'criterion'), nsteps = NULL,
                        highlight = TRUE, col = NULL,
                        label = FALSE, label.vars = NULL, ...) {

  type = match.arg(type)

  if (!.is.msaenet(x))
    stop('object class must be "msaenet"')

  if (.is.multistep(x)) {
    beta.mat = do.call(cbind, x$'beta.list')
    best.step = x$'best.step'
  }

  if (.is.adaptive(x)) {
    beta.mat = as.matrix(cbind(x$'beta.first', x$'beta'))
    best.step = 2L
  }

  if (is.null(nsteps))
    nsteps = ncol(beta.mat)

  if (type == 'coef') {
    nzv.idx = msaenet.nzv(x)
    .parcor(beta.mat, nsteps, best.step, nzv.idx,
            highlight, col, label, label.vars)
  }

  if (type == 'criterion') {
    if (is.null(x$'post.criterion'))
      stop('No post selection ICs available, since `tune.nsteps = "max"
           or it is a one-step-only adaptive model object`') else
             .scree(x$'post.criterion', nsteps, best.step, highlight)
  }

  invisible()

}

# parallel coordinates plot
.parcor = function(x, nsteps, best.step, nzv.idx,
                   highlight, col, label, label.vars) {

  x = x[, 1L:nsteps]
  xmin = min(x)
  xmax = max(x)

  # box
  matplot(1L:nsteps, t(x), xlab = '', ylab = 'Coefficients',
          xaxt = 'n', yaxt = 'n', type = 'n', axes = TRUE)

  # axes with only ticks
  axis(2, at = c(xmin, 0, xmax), labels = c('Min', '0', 'Max'), lwd = 0, lwd.ticks = 1)
  axis(1, at = 1L:nsteps, labels = paste('Step', 1L:nsteps), lwd = 0, lwd.ticks = 1)

  # step lines
  for (i in 1L:nsteps) lines(x = c(i, i), y = c(xmin, xmax), col = 'grey90')

  # highlight optimal step
  if (highlight) {
    lines(x = c(best.step, best.step), y = c(xmin, xmax),
          col = 'white', lty = 1, lwd = 1.5)
    lines(x = c(best.step, best.step), y = c(xmin, xmax),
          col = 'darkred', lty = 2, lwd = 1.5)
  }

  # coefficient paths
  if (is.null(col))  # 10-color palette from D3 (v3)
    col = c('#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
            '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF')
  matlines(1L:nsteps, t(x), lty = 1, lwd = 1.2, col = col)

  # zero line
  lines(x = c(1L, nsteps), y = c(0, 0))

  # label variables
  if (label & is.null(label.vars))
    text(x = best.step, y = x[nzv.idx, best.step],
         labels = as.character(nzv.idx),
         offset = 0.3, pos = 2, cex = 0.7)

  if (label & !is.null(label.vars)) {

    if (length(label.vars) != nrow(x))
      stop('Length of `label.vars` should be the same as the number of variables') else
        text(x = best.step, y = x[nzv.idx, best.step],
             labels = label.vars[nzv.idx],
             offset = 0.3, pos = 2, cex = 0.7)

  }

  invisible()

}

.scree = function(x, nsteps, best.step, highlight) {

  x = x[1L:nsteps]

  plot(1L:length(x), x, type = 'b', xaxt = 'n',
       xlab = '', ylab = 'Post Selection Criterion')

  axis(1, at = 1L:nsteps, labels = paste('Step', 1L:nsteps),
       lwd = 0, lwd.ticks = 1)

  if (highlight) lines(x = c(best.step, best.step),
                       y = c(min(x) - 0.5, max(x) + 0.5),
                       col = 'darkred', lty = 2, lwd = 1.5)

}
