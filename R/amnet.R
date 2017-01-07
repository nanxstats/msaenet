#' Adaptive MCP-Net
#'
#' Adaptive MCP-Net
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
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Vector of candidate \code{gamma}s to use in MCP-Net.
#' @param alphas Vector of candidate \code{alpha}s to use in MCP-Net.
#' @param eps Convergence threshhold to use in MCP-net.
#' @param max.iter Maximum number of iterations to use in MCP-net.
#' @param gamma Scaling factor for adaptive weights:
#' \code{weights = coefs^(-gamma)}.
#' @param seed Random seed for cross-validation fold division.
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
#' @export amnet
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' amnet.fit = amnet(dat$x.tr, dat$y.tr,
#'                   gammas = 3, alphas = seq(0.2, 0.8, 0.2), seed = 1002)
#'
#' print(amnet.fit)
#' msaenet.nzv(amnet.fit)
#' msaenet.fp(amnet.fit, 1:5)
#' msaenet.tp(amnet.fit, 1:5)
#' amnet.pred = predict(amnet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, amnet.pred)
#' plot(amnet.fit)

amnet = function(x, y,
                 family = c('gaussian', 'binomial', 'poisson', 'cox'),
                 init = c('mnet', 'ridge'),
                 nfolds = 5L,
                 gammas = c(1.01, 1.7, 3, 100), alphas = seq(0.05, 0.95, 0.05),
                 eps = 1e-4, max.iter = 10000L,
                 gamma = 1,
                 seed = 1001, parallel = FALSE, verbose = FALSE) {

  family = match.arg(family)
  init = match.arg(init)
  call = match.call()
  nvar = ncol(x)

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'mnet') {
    mnet.cv = msaenet.tune.ncvreg(x, y, penalty = 'MCP',
                                  nfolds = nfolds,
                                  family = family,
                                  gammas = gammas, alphas = alphas,
                                  eps = eps, max.iter = max.iter,
                                  seed = seed, parallel = parallel)
  }

  if (init == 'ridge') {
    mnet.cv = msaenet.tune.ncvreg(x, y, penalty = 'MCP',
                                  nfolds = nfolds,
                                  family = family,
                                  gammas = gammas, alphas = 1e-16,
                                  eps = eps, max.iter = max.iter,
                                  seed = seed, parallel = parallel)
  }

  best.gamma.mnet  = mnet.cv$'best.gamma'
  best.alpha.mnet  = mnet.cv$'best.alpha'
  best.lambda.mnet = mnet.cv$'best.model'$'lambda.min'

  mnet.full = .ncvnet(x, y, penalty = 'MCP',
                      gamma  = best.gamma.mnet,
                      alpha  = best.alpha.mnet,
                      lambda = best.lambda.mnet,
                      family = family,
                      eps = eps, max.iter = max.iter)

  bhat = .ncv.coef(mnet.full, nvar)
  if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  adpen = (pmax(abs(bhat), .Machine$double.eps))^(-gamma)

  if (verbose) cat('Starting step 2 ...\n')

  amnet.cv = msaenet.tune.ncvreg(x, y, penalty = 'MCP',
                                 nfolds = nfolds,
                                 family = family,
                                 penalty.factor = adpen,
                                 gammas = gammas, alphas = alphas,
                                 eps = eps, max.iter = max.iter,
                                 seed = seed + 1L, parallel = parallel)

  best.gamma.amnet  = amnet.cv$'best.gamma'
  best.alpha.amnet  = amnet.cv$'best.alpha'
  best.lambda.amnet = amnet.cv$'best.model'$'lambda.min'

  amnet.full = .ncvnet(x, y, penalty = 'MCP',
                       penalty.factor = adpen,
                       gamma  = best.gamma.amnet,
                       alpha  = best.alpha.amnet,
                       lambda = best.lambda.amnet,
                       family = family,
                       eps = eps, max.iter = max.iter)

  # final beta stored as sparse matrix
  bhat.full  = Matrix(.ncv.coef(amnet.full, nvar), sparse = TRUE)
  bhat.first = Matrix(.ncv.coef(mnet.full,  nvar), sparse = TRUE)

  amnet.model = list('beta'  = bhat.full,
                     'model' = amnet.full,
                     'beta.first'  = bhat.first,
                     'model.first' = mnet.full,
                     'best.alpha.mnet'   = best.alpha.mnet,
                     'best.alpha.amnet'  = best.alpha.amnet,
                     'best.lambda.mnet'  = best.lambda.mnet,
                     'best.lambda.amnet' = best.lambda.amnet,
                     'best.gamma.mnet'   = best.gamma.mnet,
                     'best.gamma.amnet'  = best.gamma.amnet,
                     'adpen' = adpen,
                     'seed'  = seed,
                     'call'  = call)

  class(amnet.model) = c('msaenet', 'msaenet.amnet')
  return(amnet.model)

}
