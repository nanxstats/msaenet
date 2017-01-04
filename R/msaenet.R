#' Multi-Step Adaptive Elastic-Net
#'
#' Multi-Step Adaptive Elastic-Net
#'
#' @param x Data matrix.
#' @param y Response vector.
#' @param family Response type.
#' @param init Type of the penalty used in the initial
#' estimation step. Can be \code{"enet"} or \code{"ridge"}.
#' See \code{\link[glmnet]{glmnet}} for details.
#' @param nsteps How many adaptive estimation steps? At least 2.
#' We assume adaptive elastic-net has only 1 adaptive estimation step.
#' @param nfolds Fold numbers of cross-validation.
#' @param alphas Vector of candidate \code{alpha}s to use in
#' \code{\link[glmnet]{cv.glmnet}}.
#' @param gamma Scaling factor for adaptive weights:
#' \code{weights = coefs^(-gamma)}.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}} for details.
#' @param seed Two random seeds for cross-validation fold division
#' in two estimation steps.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#' @param verbose Should we print out the estimation progress?
#'
#' @return List of coefficients \code{beta} and
#' \code{glmnet} model object \code{model}.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @references
#' Nan Xiao and Qing-Song Xu. (2015). Multi-step adaptive elastic-net:
#' reducing false positives in high-dimensional variable selection.
#' \emph{Journal of Statistical Computation and Simulation} 85(18), 3755--3765.
#'
#' @importFrom glmnet glmnet
#' @importFrom Matrix Matrix
#'
#' @export msaenet
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
#' msaenet.nzv(msaenet.fit)
#' msaenet.fp(msaenet.fit, 1:5)
#' msaenet.tp(msaenet.fit, 1:5)
#' msaenet.pred = predict(msaenet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, msaenet.pred)
#' plot(msaenet.fit)

msaenet = function(x, y,
                   family = c('gaussian', 'binomial', 'poisson',
                              'multinomial', 'cox', 'mgaussian'),
                   init = c('enet', 'ridge'),
                   nsteps = 2L, nfolds = 5L,
                   alphas = seq(0.05, 0.95, 0.05), gamma = 1,
                   rule = c('lambda.min', 'lambda.1se'),
                   seed = 1001, parallel = FALSE, verbose = FALSE) {

  if (nsteps < 2L) stop('nsteps must be an integer >= 2')

  family = match.arg(family)
  rule = match.arg(rule)
  init = match.arg(init)
  call = match.call()

  best.alphas  = rep(NA, nsteps + 1L)
  best.lambdas = rep(NA, nsteps + 1L)
  beta.list    = vector('list', nsteps + 1L)
  model.list   = vector('list', nsteps + 1L)
  adapen.list  = vector('list', nsteps)

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'enet') {
    model.cv = msaenet.tune.glmnet(x, y, family = family,
                                   nfolds = nfolds, alphas = alphas,
                                   seed = seed, parallel = parallel)
  }

  if (init == 'ridge') {
    model.cv = msaenet.tune.glmnet(x, y, family = family,
                                   nfolds = nfolds, alphas = 0,
                                   seed = seed, parallel = parallel)
  }

  best.alphas[[1L]] = model.cv$'best.alpha'

  if (rule == 'lambda.min') {
    best.lambdas[[1L]] = model.cv$'best.model'$'lambda.min'
  } else if (rule == 'lambda.1se') {
    best.lambdas[[1L]] = model.cv$'best.model'$'lambda.1se'
  }

  model.list[[1L]] = glmnet(x, y, family = family,
                            lambda = best.lambdas[[1L]],
                            alpha  = best.alphas[[1L]])

  if (model.list[[1L]][['df']] < 0.5)
    stop('Null model produced by the full fit (all coefficients are zero).
         Please try to change rule, alphas, seed, nfolds, or increase sample size.')

  bhat = as.matrix(model.list[[1L]][['beta']])
  if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))
  beta.list[[1L]] = bhat

  # MSAEnet steps
  for (i in 1L:nsteps) {

    adpen.raw  = (pmax(abs(beta.list[[i]]), .Machine$double.eps))^(-gamma)
    adapen.list[[i]] = as.vector(adpen.raw)
    adpen.name = rownames(adpen.raw)
    names(adapen.list[[i]]) = adpen.name

    if (verbose) cat('Starting step', i + 1, '...\n')

    model.cv = msaenet.tune.glmnet(x, y, family = family, nfolds = nfolds,
                                   penalty.factor = adapen.list[[i]],
                                   alphas = alphas,
                                   seed = seed + i,
                                   parallel = parallel)

    best.alphas[[i + 1L]] = model.cv$'best.alpha'

    if (rule == 'lambda.min') {
      best.lambdas[[i + 1L]] = model.cv$'best.model'$'lambda.min'
    } else if (rule == 'lambda.1se') {
      best.lambdas[[i + 1L]] = model.cv$'best.model'$'lambda.1se'
    }

    model.list[[i + 1L]] = glmnet(x, y, family = family,
                                  lambda = best.lambdas[[i + 1L]],
                                  alpha = best.alphas[[i + 1L]],
                                  penalty.factor = adapen.list[[i]])

    if (model.list[[i + 1L]][['df']] < 0.5)
      stop('Null model produced by the full fit (all coefficients are zero).
           Please try to change rule, alphas, seed, nfolds, or increase sample size.')

    bhat = as.matrix(model.list[[i + 1L]][['beta']])
    if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))
    beta.list[[i + 1L]] = bhat

  }

  msaenet.model = list('beta' = Matrix(beta.list[[nsteps + 1L]], sparse = TRUE),
                       # final beta stored as sparse matrix
                       'model' = model.list[[nsteps + 1L]],
                       # final model object
                       'best.alphas'  = best.alphas,
                       'best.lambdas' = best.lambdas,
                       'beta.list'    = beta.list,
                       'model.list'   = model.list,
                       'adapen.list'  = adapen.list,
                       'seed' = seed,
                       'call' = call)

  class(msaenet.model) = c('msaenet', 'msaenet.msaenet')
  return(msaenet.model)

}
