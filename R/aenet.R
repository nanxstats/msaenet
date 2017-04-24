#' Adaptive Elastic-Net
#'
#' Adaptive Elastic-Net
#'
#' @param x Data matrix.
#' @param y Response vector if \code{family} is \code{"gaussian"},
#' \code{"binomial"}, or \code{"poisson"}. If \code{family} is
#' \code{"cox"}, a response matrix created by \code{\link[survival]{Surv}}.
#' @param family Model family, can be \code{"gaussian"},
#' \code{"binomial"}, \code{"poisson"}, or \code{"cox"}.
#' @param init Type of the penalty used in the initial
#' estimation step. Can be \code{"enet"} or \code{"ridge"}.
#' @param alphas Vector of candidate \code{alpha}s to use in
#' \code{\link[glmnet]{cv.glmnet}}.
#' @param tune Parameter tuning method for each estimation step.
#' Possible options are \code{"cv"}, \code{"ebic"}, \code{"bic"},
#' and \code{"aic"}. Default is \code{"cv"}.
#' @param nfolds Fold numbers of cross-validation when \code{tune = "cv"}.
#' @param rule Lambda selection criterion when \code{tune = "cv"},
#' can be \code{"lambda.min"} or \code{"lambda.1se"}.
#' See \code{\link[glmnet]{cv.glmnet}} for details.
#' @param ebic.gamma Parameter for Extended BIC penalizing
#' size of the model space when \code{tune = "ebic"},
#' default is \code{1}. For details, see Chen and Chen (2008).
#' @param scale Scaling factor for adaptive weights:
#' \code{weights = coefficients^(-scale)}.
#' @param lower.limits Lower limits for coefficients.
#' Default is \code{-Inf}. For details, see \code{\link[glmnet]{glmnet}}.
#' @param upper.limits Upper limits for coefficients.
#' Default is \code{Inf}. For details, see \code{\link[glmnet]{glmnet}}.
#' @param seed Random seed for cross-validation fold division.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#' @param verbose Should we print out the estimation progress?
#'
#' @return List of model coefficients, \code{glmnet} model object,
#' and the optimal parameter set.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @references
#' Zou, Hui, and Hao Helen Zhang. (2009).
#' On the adaptive elastic-net with a diverging number of parameters.
#' \emph{The Annals of Statistics} 37(4), 1733--1751.
#'
#' @importFrom glmnet glmnet
#' @importFrom Matrix Matrix
#'
#' @export aenet
#'
#' @examples
#' dat = msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001)
#'
#' aenet.fit = aenet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2), seed = 1002)
#'
#' print(aenet.fit)
#' msaenet.nzv(aenet.fit)
#' msaenet.fp(aenet.fit, 1:5)
#' msaenet.tp(aenet.fit, 1:5)
#' aenet.pred = predict(aenet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, aenet.pred)
#' plot(aenet.fit)

aenet = function(
  x, y,
  family = c('gaussian', 'binomial', 'poisson', 'cox'),
  init   = c('enet', 'ridge'),
  alphas = seq(0.05, 0.95, 0.05),
  tune   = c('cv', 'ebic', 'bic', 'aic'),
  nfolds = 5L, rule = c('lambda.min', 'lambda.1se'),
  ebic.gamma = 1,
  scale  = 1,
  lower.limits = -Inf, upper.limits = Inf,
  seed   = 1001, parallel = FALSE, verbose = FALSE) {

  family = match.arg(family)
  init   = match.arg(init)
  tune   = match.arg(tune)
  rule   = match.arg(rule)
  call   = match.call()

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'enet') {
    enet.cv = msaenet.tune.glmnet(
      x = x, y = y, family = family,
      alphas = alphas,
      tune = tune,
      nfolds = nfolds, rule = rule,
      ebic.gamma = ebic.gamma,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      seed = seed, parallel = parallel)
  }

  if (init == 'ridge') {
    enet.cv = msaenet.tune.glmnet(
      x = x, y = y, family = family,
      alphas = 0,
      tune = tune,
      nfolds = nfolds, rule = rule,
      ebic.gamma = ebic.gamma,
      lower.limits = lower.limits,
      upper.limits = upper.limits,
      seed = seed, parallel = parallel)
  }

  best.alpha.enet     = enet.cv$'best.alpha'
  best.lambda.enet    = enet.cv$'best.lambda'
  step.criterion.enet = enet.cv$'step.criterion'

  enet.full = glmnet(
    x = x, y = y, family = family,
    alpha  = best.alpha.enet,
    lambda = best.lambda.enet,
    lower.limits = lower.limits,
    upper.limits = upper.limits)

  bhat = as.matrix(enet.full$'beta')
  if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  adpen = (pmax(abs(bhat), .Machine$double.eps))^(-scale)

  if (verbose) cat('Starting step 2 ...\n')

  aenet.cv = msaenet.tune.glmnet(
    x = x, y = y, family = family,
    alphas = alphas,
    tune = tune,
    nfolds = nfolds, rule = rule,
    ebic.gamma = ebic.gamma,
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    seed = seed + 1L, parallel = parallel,
    penalty.factor = adpen)

  best.alpha.aenet     = aenet.cv$'best.alpha'
  best.lambda.aenet    = aenet.cv$'best.lambda'
  step.criterion.aenet = aenet.cv$'step.criterion'

  aenet.full = glmnet(
    x = x, y = y, family = family,
    alpha  = best.alpha.aenet,
    lambda = best.lambda.aenet,
    lower.limits = lower.limits,
    upper.limits = upper.limits,
    penalty.factor = adpen)

  # final beta stored as sparse matrix
  bhat.full = Matrix(aenet.full$'beta', sparse = TRUE)

  aenet.model = list(
    'beta'              = bhat.full,
    'model'             = aenet.full,
    'beta.first'        = enet.full$'beta',
    'model.first'       = enet.full,
    'best.alpha.enet'   = best.alpha.enet,
    'best.alpha.aenet'  = best.alpha.aenet,
    'best.lambda.enet'  = best.lambda.enet,
    'best.lambda.aenet' = best.lambda.aenet,
    'step.criterion'    = c(step.criterion.enet, step.criterion.aenet),
    'adpen'             = adpen,
    'seed'              = seed,
    'call'              = call)

  class(aenet.model) = c('msaenet', 'msaenet.aenet')
  aenet.model

}
