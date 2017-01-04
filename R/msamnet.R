#' Multi-Step Adaptive MCP-Net
#'
#' Multi-Step Adaptive MCP-Net
#'
#' @param x Data matrix.
#' @param y Response vector if \code{family} is \code{"gaussian"},
#' \code{"binomial"} or \code{"poisson"}.
#' If \code{family} is \code{"cox"}, a response matrix made by
#' \code{\link[survival]{Surv}}.
#' @param family Model family, can be \code{"gaussian"},
#' \code{"binomial"}, \code{"poisson"}, or \code{"cox"}.
#' @param init Type of the penalty used in the initial
#' estimation step. Can be \code{"mnet"} or \code{"ridge"}.
#' @param nsteps How many adaptive estimation steps? At least 2.
#' We assume adaptive MCP-net has only 1 adaptive estimation step.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Vector of candidate \code{gamma}s to use in MCP-Net.
#' @param alphas Vector of candidate \code{alpha}s to use in MCP-Net.
#' @param gamma Scaling factor for adaptive weights:
#' \code{weights = coefs^(-gamma)}.
#' @param seed Two random seeds for cross-validation fold division
#' in two estimation steps.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#' @param verbose Should we print out the estimation progress?
#'
#' @return List of coefficients \code{beta} and
#' \code{ncvreg} model object \code{model}.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @importFrom ncvreg ncvreg ncvsurv
#' @importFrom Matrix Matrix
#'
#' @export msamnet
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' msamnet.fit = msamnet(dat$x.tr, dat$y.tr,
#'                       gammas = 3, alphas = seq(0.3, 0.9, 0.3),
#'                       nsteps = 3L, seed = 1003)
#'
#' print(msamnet.fit)
#' msaenet.nzv(msamnet.fit)
#' msaenet.fp(msamnet.fit, 1:5)
#' msaenet.tp(msamnet.fit, 1:5)
#' msamnet.pred = predict(msamnet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, msamnet.pred)
#' plot(msamnet.fit)

msamnet = function(x, y,
                   family = c('gaussian', 'binomial', 'poisson', 'cox'),
                   init = c('mnet', 'ridge'),
                   nsteps = 2L, nfolds = 5L,
                   gammas = c(1.01, 1.7, 3, 100), alphas = seq(0.05, 0.95, 0.05),
                   gamma = 1,
                   seed = 1001, parallel = FALSE, verbose = FALSE) {

  if (nsteps < 2L) stop('nsteps must be an integer >= 2')

  family = match.arg(family)
  init = match.arg(init)
  call = match.call()
  nvar = ncol(x)

  best.gammas  = rep(NA, nsteps + 1L)
  best.alphas  = rep(NA, nsteps + 1L)
  best.lambdas = rep(NA, nsteps + 1L)
  beta.list    = vector('list', nsteps + 1L)
  model.list   = vector('list', nsteps + 1L)
  adapen.list  = vector('list', nsteps)

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'mnet') {
    model.cv = msaenet.tune.ncvreg(x, y, penalty = 'MCP',
                                   nfolds = nfolds,
                                   family = family,
                                   gammas = gammas, alphas = alphas,
                                   seed = seed, parallel = parallel)
  }

  if (init == 'ridge') {
    model.cv = msaenet.tune.ncvreg(x, y, penalty = 'MCP',
                                   nfolds = nfolds,
                                   family = family,
                                   gammas = gammas, alphas = 1e-16,
                                   seed = seed, parallel = parallel)
  }

  best.gammas[[1L]]  = model.cv$'best.gamma'
  best.alphas[[1L]]  = model.cv$'best.alpha'
  best.lambdas[[1L]] = model.cv$'best.model'$'lambda.min'

  model.list[[1L]] = .ncvnet(x, y, penalty = 'MCP',
                             gamma  = best.gammas[[1L]],
                             alpha  = best.alphas[[1L]],
                             lambda = best.lambdas[[1L]],
                             family = family)

  if (.ncvdf(model.list[[1L]]) < 0.5)
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to change gammas, alphas, seed, nfolds, or increase sample size.')

  bhat = .ncv.coef(model.list[[1L]], nvar)
  if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))
  beta.list[[1L]] = bhat

  # MSAMNet steps
  for (i in 1L:nsteps) {

    adpen.raw = (pmax(abs(beta.list[[i]]), .Machine$double.eps))^(-gamma)
    adapen.list[[i]] = as.vector(adpen.raw)

    if (verbose) cat('Starting step', i + 1, '...\n')

    model.cv = msaenet.tune.ncvreg(x, y, penalty = 'MCP',
                                   nfolds = nfolds,
                                   family = family,
                                   penalty.factor = adapen.list[[i]],
                                   gammas = gammas, alphas = alphas,
                                   seed = seed + i, parallel = parallel)

    best.gammas[[i + 1L]]  = model.cv$'best.gamma'
    best.alphas[[i + 1L]]  = model.cv$'best.alpha'
    best.lambdas[[i + 1L]] = model.cv$'best.model'$'lambda.min'

    model.list[[i + 1L]] = .ncvnet(x, y, penalty = 'MCP',
                                   penalty.factor = adapen.list[[i]],
                                   gamma  = best.gammas[[i + 1L]],
                                   alpha  = best.alphas[[i + 1L]],
                                   lambda = best.lambdas[[i + 1L]],
                                   family = family)

    if (.ncvdf(model.list[[i + 1L]]) < 0.5)
      stop('Null model produced by the full fit (all coefficients are zero).
           Please try to change gammas, alphas, seed, nfolds, or increase sample size.')

    bhat = .ncv.coef(model.list[[i + 1L]], nvar)
    if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))
    beta.list[[i + 1L]] = bhat

  }

  msamnet.model = list('beta' = Matrix(beta.list[[nsteps + 1L]], sparse = TRUE),
                       # final beta stored as sparse matrix
                       'model' = model.list[[nsteps + 1L]],
                       # final model object
                       'best.alphas'  = best.alphas,
                       'best.gammas'  = best.gammas,
                       'best.lambdas' = best.lambdas,
                       'beta.list'    = beta.list,
                       'model.list'   = model.list,
                       'adapen.list'  = adapen.list,
                       'seed' = seed,
                       'call' = call)

  class(msamnet.model) = c('msaenet', 'msaenet.msamnet')
  return(msamnet.model)

}
