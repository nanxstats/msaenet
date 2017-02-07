#' Adaptive SCAD-Net
#'
#' Adaptive SCAD-Net
#'
#' @param x Data matrix.
#' @param y Response vector if \code{family} is \code{"gaussian"},
#' \code{"binomial"}, or \code{"poisson"}. If \code{family} is
#' \code{"cox"}, a response matrix created by \code{\link[survival]{Surv}}.
#' @param family Model family, can be \code{"gaussian"},
#' \code{"binomial"}, \code{"poisson"}, or \code{"cox"}.
#' @param init Type of the penalty used in the initial
#' estimation step. Can be \code{"snet"} or \code{"ridge"}.
#' @param gammas Vector of candidate \code{gamma}s (the concavity parameter)
#' to use in SCAD-Net. Default is \code{3.7}.
#' @param alphas Vector of candidate \code{alpha}s to use in SCAD-Net.
#' @param tune Parameter tuning method for each estimation step.
#' Possible options are \code{"cv"}, \code{"ebic"}, \code{"bic"},
#' and \code{"aic"}. Default is \code{"cv"}.
#' @param nfolds Fold numbers of cross-validation when \code{tune = "cv"}.
#' @param ebic.gamma Parameter for Extended BIC penalizing
#' size of the model space when \code{tune = "ebic"},
#' default is \code{1}. For details, see Chen and Chen (2008).
#' @param scale Scaling factor for adaptive weights:
#' \code{weights = coefficients^(-scale)}.
#' @param eps Convergence threshhold to use in SCAD-net.
#' @param max.iter Maximum number of iterations to use in SCAD-net.
#' @param seed Random seed for cross-validation fold division.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#' @param verbose Should we print out the estimation progress?
#'
#' @return List of model coefficients, \code{ncvreg} model object,
#' and the optimal parameter set.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom ncvreg ncvreg ncvsurv
#' @importFrom Matrix Matrix
#'
#' @export asnet
#'
#' @examples
#' dat = msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001)
#'
#' asnet.fit = asnet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2), seed = 1002)
#'
#' print(asnet.fit)
#' msaenet.nzv(asnet.fit)
#' msaenet.fp(asnet.fit, 1:5)
#' msaenet.tp(asnet.fit, 1:5)
#' asnet.pred = predict(asnet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, asnet.pred)
#' plot(asnet.fit)

asnet = function(x, y,
                 family = c('gaussian', 'binomial', 'poisson', 'cox'),
                 init = c('snet', 'ridge'),
                 gammas = 3.7, alphas = seq(0.05, 0.95, 0.05),
                 tune = c('cv', 'ebic', 'bic', 'aic'),
                 nfolds = 5L,
                 ebic.gamma = 1,
                 scale = 1,
                 eps = 1e-4, max.iter = 10000L,
                 seed = 1001, parallel = FALSE, verbose = FALSE) {

  family = match.arg(family)
  init   = match.arg(init)
  tune   = match.arg(tune)
  call   = match.call()
  nvar   = ncol(x)

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'snet') {

    snet.cv = msaenet.tune.ncvreg(
      x = x, y = y, family = family, penalty = 'SCAD',
      gammas = gammas, alphas = alphas,
      tune = tune,
      nfolds = nfolds,
      ebic.gamma = ebic.gamma,
      eps = eps, max.iter = max.iter,
      seed = seed, parallel = parallel)

    best.gamma.snet     = snet.cv$'best.gamma'
    best.alpha.snet     = snet.cv$'best.alpha'
    best.lambda.snet    = snet.cv$'best.lambda'
    step.criterion.snet = snet.cv$'step.criterion'

    snet.full = .ncvnet(
      x = x, y = y, family = family, penalty = 'SCAD',
      gamma  = best.gamma.snet,
      alpha  = best.alpha.snet,
      lambda = best.lambda.snet,
      eps = eps, max.iter = max.iter)

    bhat = .coef.ncvreg(snet.full, nvar)

  }

  if (init == 'ridge') {

    snet.cv = msaenet.tune.glmnet(
      x = x, y = y, family = family,
      alphas = 0,
      tune = tune,
      nfolds = nfolds, rule = 'lambda.min',
      ebic.gamma = ebic.gamma,
      seed = seed, parallel = parallel)

    best.gamma.snet     = NA
    best.alpha.snet     = snet.cv$'best.alpha'
    best.lambda.snet    = snet.cv$'best.lambda'
    step.criterion.snet = snet.cv$'step.criterion'

    snet.full = glmnet(
      x = x, y = y, family = family,
      alpha  = best.alpha.snet,
      lambda = best.lambda.snet)

    bhat = as.matrix(snet.full$'beta')

  }

  if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  adpen = (pmax(abs(bhat), .Machine$double.eps))^(-scale)

  if (verbose) cat('Starting step 2 ...\n')

  asnet.cv = msaenet.tune.ncvreg(
    x = x, y = y, family = family, penalty = 'SCAD',
    gammas = gammas, alphas = alphas,
    tune = tune,
    nfolds = nfolds,
    ebic.gamma = ebic.gamma,
    eps = eps, max.iter = max.iter,
    seed = seed + 1L, parallel = parallel,
    penalty.factor = adpen)

  best.gamma.asnet     = asnet.cv$'best.gamma'
  best.alpha.asnet     = asnet.cv$'best.alpha'
  best.lambda.asnet    = asnet.cv$'best.lambda'
  step.criterion.asnet = asnet.cv$'step.criterion'

  asnet.full = .ncvnet(
    x = x, y = y, family = family, penalty = 'SCAD',
    gamma  = best.gamma.asnet,
    alpha  = best.alpha.asnet,
    lambda = best.lambda.asnet,
    eps = eps, max.iter = max.iter,
    penalty.factor = adpen)

  # final beta stored as sparse matrix
  bhat.full  = Matrix(.coef.ncvreg(asnet.full, nvar), sparse = TRUE)
  bhat.first = Matrix(.coef.ncvreg(snet.full,  nvar), sparse = TRUE)

  asnet.model = list(
    'beta'              = bhat.full,
    'model'             = asnet.full,
    'beta.first'        = bhat.first,
    'model.first'       = snet.full,
    'best.alpha.snet'   = best.alpha.snet,
    'best.alpha.asnet'  = best.alpha.asnet,
    'best.lambda.snet'  = best.lambda.snet,
    'best.lambda.asnet' = best.lambda.asnet,
    'best.gamma.snet'   = best.gamma.snet,
    'best.gamma.asnet'  = best.gamma.asnet,
    'step.criterion'    = c(step.criterion.snet, step.criterion.asnet),
    'adpen'             = adpen,
    'seed'              = seed,
    'call'              = call)

  class(asnet.model) = c('msaenet', 'msaenet.asnet')
  asnet.model

}
