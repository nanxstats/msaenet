#' Plot msaenet Model Objects
#'
#' Plot msaenet model objects.
#'
#' @param x An object of class \code{msaenet} produced
#' by \code{\link{aenet}}, \code{amnet}, \code{asnet},
#' \code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.
#' @param type Plot type, \code{"coef"} for a coefficient path plot
#' across all estimation steps; \code{"criterion"} for a scree plot of
#' the model evaluation criterion used (CV error, AIC, BIC, or EBIC);
#' \code{"dotplot"} for a Cleveland dot plot of the coefficients
#' estimated by the model at the optimal step.
#' @param nsteps Maximum number of estimation steps to plot.
#' Default is to plot all steps.
#' @param highlight Should we highlight the "optimal" step
#' according to the criterion? Default is \code{TRUE}.
#' @param col Color palette to use for the coefficient paths.
#' If it is \code{NULL}, a default color palette will be assigned.
#' @param label Should we label all the non-zero variables of the
#' optimal step in the coefficient plot or the dot plot?
#' Default is \code{FALSE}. If \code{TRUE} and \code{label.vars = NULL},
#' the index of the non-zero variables will be used as labels.
#' @param label.vars Labels to use for all the variables
#' if \code{label = "TRUE"}.
#' @param label.pos Position of the labels. See argument
#' \code{pos} in \code{\link[graphics]{text}} for details.
#' @param label.offset Offset of the labels. See argument
#' \code{offset} in \code{\link[graphics]{text}} for details.
#' @param label.cex Character expansion factor of the labels.
#' See argument \code{cex} in \code{\link[graphics]{text}} for details.
#' @param label.srt Label rotation in degrees for the Cleveland dot plot.
#' Default is \code{90}. See argument \code{srt} in
#' \code{\link[graphics]{par}} for details.
#' @param xlab Title for x axis. If is \code{NULL}, will use the default title.
#' @param ylab Title for y axis. If is \code{NULL}, will use the default title.
#' @param abs Should we plot the absolute values of the coefficients
#' instead of the raw coefficients in the Cleveland dot plot?
#' Default is \code{FALSE}.
#' @param ... Other parameters (not used).
#'
#' @method plot msaenet
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom graphics axis lines matlines matplot plot text abline points
#' @importFrom stats coef
#'
#' @export
#'
#' @examples
#' dat = msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001)
#'
#' msasnet.fit = msasnet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2),
#'   nsteps = 5L, tune.nsteps = "ebic",
#'   seed = 1003)
#'
#' plot(msasnet.fit)
#' plot(msasnet.fit, label = TRUE)
#' plot(msasnet.fit, label = TRUE, nsteps = 5)
#' plot(msasnet.fit, type = "criterion")
#' plot(msasnet.fit, type = "criterion", nsteps = 5)
#' plot(msasnet.fit, type = "dotplot", label = TRUE)
#' plot(msasnet.fit, type = "dotplot", label = TRUE, abs = TRUE)

plot.msaenet = function(
  x, type = c('coef', 'criterion', 'dotplot'), nsteps = NULL,
  highlight = TRUE, col = NULL,
  label = FALSE, label.vars = NULL,
  label.pos = 2, label.offset = 0.3, label.cex = 0.7,
  label.srt = 90,
  xlab = NULL, ylab = NULL,
  abs = FALSE, ...) {

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
    .parcor(
      beta.mat, nsteps, best.step, nzv.idx,
      highlight, col,
      label, label.vars,
      label.pos, label.offset, label.cex,
      xlab, ylab)
  }

  if (type == 'criterion') {
    if (is.null(x$'post.criterion'))
      stop('No post selection ICs available, since `tune.nsteps = "max"
           or it is a one-step-only adaptive model object`') else
             .scree(x$'post.criterion',
                    nsteps, best.step, highlight,
                    xlab, ylab)
  }

  if (type == 'dotplot') {
    .dotplot(
      x, abs,
      label, label.vars, label.cex, label.srt,
      xlab, ylab)
  }

  invisible()

}

# parallel coordinates plot
.parcor = function(
  x, nsteps, best.step, nzv.idx,
  highlight, col,
  label, label.vars,
  label.pos, label.offset, label.cex,
  xlab, ylab) {

  x = x[, 1L:nsteps]
  xmin = min(x)
  xmax = max(x)

  if (is.null(xlab)) xlab = 'Number of Estimation Steps'
  if (is.null(ylab)) ylab = 'Coefficients'

  # box
  matplot(1L:nsteps, t(x), xlab = xlab, ylab = ylab,
          xaxt = 'n', yaxt = 'n', type = 'n', axes = TRUE)

  # axes with only ticks
  axis(2, at = c(xmin, 0, xmax), labels = c('Min', '0', 'Max'), lwd = 0, lwd.ticks = 1)
  axis(1, at = 1L:nsteps, labels = as.character(1L:nsteps), lwd = 0, lwd.ticks = 1)

  # step lines
  for (i in 1L:nsteps)
    lines(x = c(i, i), y = c(xmin - 42, xmax + 42), lty = 3, col = 'grey70')

  # highlight optimal step
  if (highlight) {
    lines(x = c(best.step, best.step), y = c(xmin - 42, xmax + 42),
          col = 'white', lty = 1, lwd = 1.5)
    lines(x = c(best.step, best.step), y = c(xmin - 42, xmax + 42),
          col = 'darkred', lty = 2, lwd = 1.5)
  }

  # coefficient paths
  if (is.null(col))  # 10-color palette from D3 (v3)
    col = c('#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
            '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF')
  matlines(1L:nsteps, t(x), lty = 1, lwd = 1.2, col = col)

  # zero line
  lines(x = c(0L, nsteps + 1L), y = c(0, 0), lwd = 1.2)

  # label variables
  if (label & is.null(label.vars))
    text(x = best.step, y = x[nzv.idx, best.step],
         labels = as.character(nzv.idx),
         pos = label.pos, offset = label.offset, cex = label.cex)

  if (label & !is.null(label.vars)) {

    if (length(label.vars) != nrow(x))
      stop('Length of `label.vars` should be the same as the number of variables') else
        text(x = best.step, y = x[nzv.idx, best.step],
             labels = label.vars[nzv.idx],
             pos = label.pos, offset = label.offset, cex = label.cex)

  }

  invisible()

}

.scree = function(x, nsteps, best.step, highlight, xlab, ylab) {

  x = x[1L:nsteps]

  if (is.null(xlab)) xlab = 'Number of Estimation Steps'
  if (is.null(ylab)) ylab = 'Model Selection Criterion'

  plot(1L:length(x), x, type = 'b', xaxt = 'n', xlab = xlab, ylab = ylab)

  axis(1, at = 1L:nsteps, labels = as.character(1L:nsteps),
       lwd = 0, lwd.ticks = 1)

  if (highlight) lines(
    x = c(best.step, best.step),
    y = c(min(x) - 0.5, max(x) + 0.5),
    col = 'darkred', lty = 2, lwd = 1.5)

}

# Cleveland dot plot for model coefficients at the optimal step
.dotplot = function(
  x, abs, label, label.vars, label.cex, label.srt, xlab, ylab) {

  idx.nzv = msaenet.nzv(x)

  if (is.null(xlab)) xlab = 'Selected Variables'
  if (is.null(ylab)) ylab = 'Coefficients'

  if (label & is.null(label.vars)) {
    label.nzv = as.character(idx.nzv)
  } else if (label & !is.null(label.vars)) {
    label.nzv = label.vars[idx.nzv]
  } else {
    label.nzv = rep('', length(idx.nzv))
  }

  coef.nzv = coef(x)[idx.nzv]
  if (abs) coef.nzv = abs(coef.nzv)

  ord.nzv = order(coef.nzv, decreasing = TRUE)

  if (all(coef.nzv > 0)) {
    plot(coef.nzv[ord.nzv], type = 'h', xaxt = 'n', xlab = xlab, ylab = ylab,
         ylim = c(-0.3, max(coef.nzv) + 0.05))
  } else if (all(coef.nzv < 0)) {
    plot(coef.nzv[ord.nzv], type = 'h', xaxt = 'n', xlab = xlab, ylab = ylab,
         ylim = c(min(coef.nzv) - 0.05, 0.3))
  } else {
    plot(coef.nzv[ord.nzv], type = 'h', xaxt = 'n', xlab = xlab, ylab = ylab)
  }

  abline(a = 0, b = 0)

  col.vec = ifelse(coef.nzv > 0, '#BC3C29', '#0072B5')
  points(
    1L:length(coef.nzv), coef.nzv[ord.nzv], pch = 21,
    bg = col.vec[ord.nzv], col = col.vec[ord.nzv])

  idx.pos = which(coef.nzv[ord.nzv] > 0)
  idx.neg = which(coef.nzv[ord.nzv] <= 0)

  if (length(idx.pos) > 0L)
    text(
      idx.pos, -0.05, labels = label.nzv[ord.nzv][idx.pos],
      srt = label.srt, cex = label.cex, adj = 1)

  if (length(idx.neg) > 0L)
    text(
      idx.neg, 0.05, labels = label.nzv[ord.nzv][idx.neg],
      srt = label.srt, cex = label.cex, adj = 0)

}
