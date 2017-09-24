#' Multi-Step Adaptive MCP-Net
#'
#' Multi-Step Adaptive MCP-Net
#'
#' @param x Data matrix.
#' @param y Response vector if \code{family} is \code{"gaussian"},
#' \code{"binomial"}, or \code{"poisson"}. If \code{family} is
#' \code{"cox"}, a response matrix created by \code{\link[survival]{Surv}}.
#' @param family Model family, can be \code{"gaussian"},
#' \code{"binomial"}, \code{"poisson"}, or \code{"cox"}.
#' @param init Type of the penalty used in the initial
#' estimation step. Can be \code{"mnet"} or \code{"ridge"}.
#' @param gammas Vector of candidate \code{gamma}s (the concavity parameter)
#' to use in MCP-Net. Default is \code{3}.
#' @param alphas Vector of candidate \code{alpha}s to use in MCP-Net.
#' @param tune Parameter tuning method for each estimation step.
#' Possible options are \code{"cv"}, \code{"ebic"}, \code{"bic"},
#' and \code{"aic"}. Default is \code{"cv"}.
#' @param nfolds Fold numbers of cross-validation when \code{tune = "cv"}.
#' @param ebic.gamma Parameter for Extended BIC penalizing
#' size of the model space when \code{tune = "ebic"},
#' default is \code{1}. For details, see Chen and Chen (2008).
#' @param nsteps Maximum number of adaptive estimation steps.
#' At least \code{2}, assuming adaptive MCP-net has only
#' one adaptive estimation step.
#' @param tune.nsteps Optimal step number selection method
#' (aggregate the optimal model from the each step and compare).
#' Options include \code{"max"} (select the final-step model directly),
#' or compare these models using \code{"ebic"}, \code{"bic"}, or \code{"aic"}.
#' Default is \code{"max"}.
#' @param ebic.gamma.nsteps Parameter for Extended BIC penalizing
#' size of the model space when \code{tune.nsteps = "ebic"},
#' default is \code{1}.
#' @param scale Scaling factor for adaptive weights:
#' \code{weights = coefficients^(-scale)}.
#' @param eps Convergence threshhold to use in MCP-net.
#' @param max.iter Maximum number of iterations to use in MCP-net.
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
#' @export msamnet
#'
#' @examples
#' dat = msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001)
#'
#' msamnet.fit = msamnet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.3, 0.9, 0.3),
#'   nsteps = 3L, seed = 1003)
#'
#' print(msamnet.fit)
#' msaenet.nzv(msamnet.fit)
#' msaenet.fp(msamnet.fit, 1:5)
#' msaenet.tp(msamnet.fit, 1:5)
#' msamnet.pred = predict(msamnet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, msamnet.pred)
#' plot(msamnet.fit)

msamnet = function(
  x, y,
  family = c('gaussian', 'binomial', 'poisson', 'cox'),
  init   = c('mnet', 'ridge'),
  gammas = 3, alphas = seq(0.05, 0.95, 0.05),
  tune   = c('cv', 'ebic', 'bic', 'aic'),
  nfolds = 5L,
  ebic.gamma = 1,
  nsteps = 2L,
  tune.nsteps = c('max', 'ebic', 'bic', 'aic'),
  ebic.gamma.nsteps = 1,
  scale  = 1,
  eps    = 1e-4, max.iter = 10000L,
  seed   = 1001, parallel = FALSE, verbose = FALSE) {

  if (nsteps < 2L) stop('nsteps must be an integer >= 2')

  family      = match.arg(family)
  init        = match.arg(init)
  tune        = match.arg(tune)
  tune.nsteps = match.arg(tune.nsteps)
  call        = match.call()
  nvar        = ncol(x)

  best.gammas    = rep(NA, nsteps + 1L)
  best.alphas    = rep(NA, nsteps + 1L)
  best.lambdas   = rep(NA, nsteps + 1L)
  step.criterion = rep(NA, nsteps + 1L)
  beta.list      = vector('list', nsteps + 1L)
  model.list     = vector('list', nsteps + 1L)
  adapen.list    = vector('list', nsteps)

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'mnet') {

    model.cv = msaenet.tune.ncvreg(
      x = x, y = y, family = family, penalty = 'MCP',
      gammas = gammas, alphas = alphas,
      tune = tune,
      nfolds = nfolds,
      ebic.gamma = ebic.gamma,
      eps = eps, max.iter = max.iter,
      seed = seed, parallel = parallel)

    best.gammas[[1L]]    = model.cv$'best.gamma'
    best.alphas[[1L]]    = model.cv$'best.alpha'
    best.lambdas[[1L]]   = model.cv$'best.lambda'
    step.criterion[[1L]] = model.cv$'step.criterion'

    model.list[[1L]] = .ncvnet(
      x = x, y = y, family = family, penalty = 'MCP',
      gamma  = best.gammas[[1L]],
      alpha  = best.alphas[[1L]],
      lambda = best.lambdas[[1L]],
      eps = eps, max.iter = max.iter)

    if (.df(model.list[[1L]]) < 0.5)
      stop('Null model produced by the full fit (all coefficients are zero).
           Please try a different parameter setting.')

    bhat = .coef.ncvreg(model.list[[1L]], nvar)

  }

  if (init == 'ridge') {

    model.cv = msaenet.tune.glmnet(
      x = x, y = y, family = family,
      alphas = 0,
      tune = tune,
      nfolds = nfolds, rule = 'lambda.min',
      ebic.gamma = ebic.gamma,
      lower.limits = -Inf, upper.limits = Inf,
      seed = seed, parallel = parallel)

    best.gammas[[1L]]    = NA
    best.alphas[[1L]]    = model.cv$'best.alpha'
    best.lambdas[[1L]]   = model.cv$'best.lambda'
    step.criterion[[1L]] = model.cv$'step.criterion'

    model.list[[1L]] = glmnet(
      x = x, y = y, family = family,
      alpha  = best.alphas[[1L]],
      lambda = best.lambdas[[1L]])

    if (.df(model.list[[1L]]) < 0.5)
      stop('Null model produced by the full fit (all coefficients are zero).
           Please try a different parameter setting.')

    bhat = as.matrix(model.list[[1L]][['beta']])

  }

  if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))
  beta.list[[1L]] = bhat

  # MSAMNet steps
  for (i in 1L:nsteps) {

    adpen.raw = (pmax(abs(beta.list[[i]]), .Machine$double.eps))^(-scale)
    adapen.list[[i]] = as.vector(adpen.raw)

    if (verbose) cat('Starting step', i + 1, '...\n')

    model.cv = msaenet.tune.ncvreg(
      x = x, y = y, family = family, penalty = 'MCP',
      gammas = gammas, alphas = alphas,
      tune = tune,
      nfolds = nfolds,
      ebic.gamma = ebic.gamma,
      eps = eps, max.iter = max.iter,
      seed = seed + i, parallel = parallel,
      penalty.factor = adapen.list[[i]])

    best.gammas[[i + 1L]]    = model.cv$'best.gamma'
    best.alphas[[i + 1L]]    = model.cv$'best.alpha'
    best.lambdas[[i + 1L]]   = model.cv$'best.lambda'
    step.criterion[[i + 1L]] = model.cv$'step.criterion'

    model.list[[i + 1L]] = .ncvnet(
      x = x, y = y, family = family, penalty = 'MCP',
      gamma  = best.gammas[[i + 1L]],
      alpha  = best.alphas[[i + 1L]],
      lambda = best.lambdas[[i + 1L]],
      eps = eps, max.iter = max.iter,
      penalty.factor = adapen.list[[i]])

    if (.df(model.list[[i + 1L]]) < 0.5)
      stop('Null model produced by the full fit (all coefficients are zero).
           Please try a different parameter setting.')

    bhat = .coef.ncvreg(model.list[[i + 1L]], nvar)
    if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))
    beta.list[[i + 1L]] = bhat

  }

  # select optimal step
  post.ics = msaenet.tune.nsteps.ncvreg(
    model.list, tune.nsteps, ebic.gamma.nsteps)

  best.step = post.ics$'best.step'
  post.criterion = post.ics$'ics'

  msamnet.model = list(
    'beta'           = Matrix(beta.list[[best.step]], sparse = TRUE),
    'model'          = model.list[[best.step]],
    'best.step'      = best.step,
    'best.alphas'    = best.alphas,
    'best.gammas'    = best.gammas,
    'best.lambdas'   = best.lambdas,
    'step.criterion' = step.criterion,
    'post.criterion' = post.criterion,
    'beta.list'      = beta.list,
    'model.list'     = model.list,
    'adapen.list'    = adapen.list,
    'seed'           = seed,
    'call'           = call)

  class(msamnet.model) = c('msaenet', 'msaenet.msamnet')
  msamnet.model

}
