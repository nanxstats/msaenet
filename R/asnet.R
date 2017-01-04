#' Adaptive SCAD-Net
#'
#' Adaptive SCAD-Net
#'
#' @param x Data matrix.
#' @param y Response vector if \code{family} is \code{"gaussian"},
#' \code{"binomial"} or \code{"poisson"}.
#' If \code{family} is \code{"cox"}, a response matrix made by
#' \code{\link[survival]{Surv}}.
#' @param family Model family, can be \code{"gaussian"},
#' \code{"binomial"}, \code{"poisson"}, or \code{"cox"}.
#' @param init Type of the penalty used in the initial
#' estimation step. Can be \code{"snet"} or \code{"ridge"}.
#' @param nfolds Fold numbers of cross-validation.
#' @param gammas Vector of candidate \code{gamma}s to use in SCAD-Net.
#' @param alphas Vector of candidate \code{alpha}s to use in SCAD-Net.
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
#' @export asnet
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' asnet.fit = asnet(dat$x.tr, dat$y.tr,
#'                   gammas = 3.7, alphas = seq(0.2, 0.8, 0.2), seed = 1002)
#'
#' print(asnet.fit)
#' msaenet.nzv(asnet.fit)
#' msaenet.fp(asnet.fit, 1:5)
#' msaenet.tp(asnet.fit, 1:5)
#' asnet.pred = predict(asnet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, asnet.pred)

asnet = function(x, y,
                 family = c('gaussian', 'binomial', 'poisson', 'cox'),
                 init = c('snet', 'ridge'),
                 nfolds = 5L,
                 gammas = c(2.01, 2.3, 3.7, 200), alphas = seq(0.05, 0.95, 0.05),
                 gamma = 1,
                 seed = 1001, parallel = FALSE, verbose = FALSE) {

  family = match.arg(family)
  init = match.arg(init)
  call = match.call()
  nvar = ncol(x)

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'snet') {
    snet.cv = msaenet.tune.ncvreg(x, y, penalty = 'SCAD',
                                  nfolds = nfolds,
                                  family = family,
                                  gammas = gammas, alphas = alphas,
                                  seed = seed, parallel = parallel)
  }

  if (init == 'ridge') {
    snet.cv = msaenet.tune.ncvreg(x, y, penalty = 'SCAD',
                                  nfolds = nfolds,
                                  family = family,
                                  gammas = gammas, alphas = 1e-16,
                                  seed = seed, parallel = parallel)
  }

  best.gamma.snet  = snet.cv$'best.gamma'
  best.alpha.snet  = snet.cv$'best.alpha'
  best.lambda.snet = snet.cv$'best.model'$'lambda.min'

  snet.full = .ncvnet(x, y, penalty = 'SCAD',
                      gamma  = best.gamma.snet,
                      alpha  = best.alpha.snet,
                      lambda = best.lambda.snet,
                      family = family)

  bhat = .ncv.coef(snet.full, nvar)
  if (all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  adpen = (pmax(abs(bhat), .Machine$double.eps))^(-gamma)

  if (verbose) cat('Starting step 2 ...\n')

  asnet.cv = msaenet.tune.ncvreg(x, y, penalty = 'SCAD',
                                 nfolds = nfolds,
                                 family = family,
                                 penalty.factor = adpen,
                                 gammas = gammas, alphas = alphas,
                                 seed = seed + 1L, parallel = parallel)

  best.gamma.asnet  = asnet.cv$'best.gamma'
  best.alpha.asnet  = asnet.cv$'best.alpha'
  best.lambda.asnet = asnet.cv$'best.model'$'lambda.min'

  asnet.full = .ncvnet(x, y, penalty = 'SCAD',
                       penalty.factor = adpen,
                       gamma  = best.gamma.asnet,
                       alpha  = best.alpha.asnet,
                       lambda = best.lambda.asnet,
                       family = family)

  # final beta stored as sparse matrix
  bhat.full  = Matrix(.ncv.coef(asnet.full, nvar), sparse = TRUE)
  bhat.first = Matrix(.ncv.coef(snet.full,  nvar), sparse = TRUE)

  asnet.model = list('beta'  = bhat.full,
                     'model' = asnet.full,
                     'beta.first'  = bhat.first,
                     'model.first' = snet.full,
                     'best.alpha.snet'   = best.alpha.snet,
                     'best.alpha.asnet'  = best.alpha.asnet,
                     'best.lambda.snet'  = best.lambda.snet,
                     'best.lambda.asnet' = best.lambda.asnet,
                     'best.gamma.snet'   = best.gamma.snet,
                     'best.gamma.asnet'  = best.gamma.asnet,
                     'adpen' = adpen,
                     'seed'  = seed,
                     'call'  = call)

  class(asnet.model) = c('msaenet', 'msaenet.asnet')
  return(asnet.model)

}
