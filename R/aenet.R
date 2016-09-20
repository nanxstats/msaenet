#' Adaptive Elastic-Net
#'
#' Adaptive Elastic-Net
#'
#' @param x Data matrix.
#' @param y Response vector.
#' @param family Response type.
#' See \code{\link[glmnet]{glmnet}} for details.
#' @param init Type of the penalty used in the initial
#' estimation step. Could be \code{"enet"} or \code{"ridge"}.
#' @param nfolds Fold numbers of cross-validation.
#' @param alphas Vector of alphas to use in \code{\link[glmnet]{cv.glmnet}}.
#' @param gamma Scaling factor for adaptive weights:
#' \code{weights = coefs^(-gamma)}.
#' @param rule Model selection criterion, \code{"lambda.min"} or
#' \code{"lambda.1se"}. See \code{\link[glmnet]{cv.glmnet}} for details.
#' @param seed Random seed for cross-validation fold division.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is {FALSE}. To enable parallel tuning, load the
#' \code{doParallel} package and run \code{registerDoParallel()}
#' with the number of CPU cores before calling this function.
#' @param verbose Should we print out the progress?
#'
#' @return List of coefficients \code{beta} and
#' \code{glmnet} model object \code{model}.
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @references
#' Zou, Hui, and Hao Helen Zhang. (2009).
#' On the Adaptive Elastic-Net with a Diverging Number of Parameters.
#' \emph{The Annals of Statistics} 37(4), 1733--51.
#'
#' @importFrom glmnet glmnet
#'
#' @export aenet
#'
#' @examples
#' dat = msaenet.sim.gaussian(n = 150, p = 500, rho = 0.6,
#'                            coef = rep(1, 5), snr = 2, p.train = 0.7,
#'                            seed = 1001)
#'
#' aenet.fit = aenet(dat$x.tr, dat$y.tr,
#'                   alphas = seq(0.2, 0.8, 0.2), seed = 1002)
#'
#' print(aenet.fit)
#' msaenet.nzv(aenet.fit)
#' msaenet.fp(aenet.fit, 1:5)
#' msaenet.tp(aenet.fit, 1:5)
#' aenet.pred = predict(aenet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, aenet.pred)

aenet = function(x, y,
                 family = c('gaussian', 'binomial', 'poisson',
                            'multinomial', 'cox', 'mgaussian'),
                 init = c('enet', 'ridge'),
                 nfolds = 5L,
                 alphas = seq(0.05, 0.95, 0.05), gamma = 1,
                 rule = c('lambda.min', 'lambda.1se'),
                 seed = 1001, parallel = FALSE, verbose = FALSE) {

  family = match.arg(family)
  rule = match.arg(rule)
  init = match.arg(init)
  call = match.call()

  if (verbose) cat('Starting step 1 ...\n')

  if (init == 'enet') {
    enet.y = msaenet.tune.glmnet.alpha(x, y, family = family,
                                       nfolds = nfolds, alphas = alphas,
                                       seed = seed, parallel = parallel)
  }

  if (init == 'ridge') {
    enet.y = msaenet.tune.glmnet.alpha(x, y, family = family,
                                       nfolds = nfolds, alphas = 0,
                                       seed = seed, parallel = parallel)
  }

  best.alpha.enet = enet.y$'best.alpha'

  if (rule == 'lambda.min') {
    best.lambda.enet = enet.y$'best.model'$'lambda.min'
  } else if (rule == 'lambda.1se') {
    best.lambda.enet = enet.y$'best.model'$'lambda.1se'
  }

  enet.all = glmnet(x, y, family = family,
                    lambda = best.lambda.enet,
                    alpha  = best.alpha.enet)

  bhat = as.matrix(enet.all$beta)
  if(all(bhat == 0)) bhat = rep(.Machine$double.eps * 2, length(bhat))

  adpen = (pmax(abs(bhat), .Machine$double.eps))^(-gamma)

  if (verbose) cat('Starting step 2 ...\n')

  aenet.y = msaenet.tune.glmnet.alpha(x, y, family = family, nfolds = nfolds,
                                      penalty.factor = adpen,
                                      alphas = alphas,
                                      seed = seed + 1L,
                                      parallel = parallel)

  best.alpha.aenet = aenet.y$'best.alpha'

  if (rule == 'lambda.min') {
    best.lambda.aenet = aenet.y$'best.model'$'lambda.min'
  } else if (rule == 'lambda.1se') {
    best.lambda.aenet = aenet.y$'best.model'$'lambda.1se'
  }

  aenet.all = glmnet(x, y, family = family,
                     lambda = best.lambda.aenet,
                     penalty.factor = adpen,
                     alpha  = best.alpha.aenet)

  aenet.model = list('beta' = aenet.all$'beta',
                     'model' = aenet.all,
                     'best.alpha.enet' = best.alpha.enet,
                     'best.alpha.aenet' = best.alpha.aenet,
                     'best.lambda.enet' = best.lambda.enet,
                     'best.lambda.aenet' = best.lambda.aenet,
                     'beta.first' = enet.all$'beta',
                     'model.first' = enet.all,
                     'adpen' = adpen,
                     'seed' = seed,
                     'call' = call)

  class(aenet.model) = c('msaenet', 'msaenet.aenet')
  return(aenet.model)

}
