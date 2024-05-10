#' Adaptive MCP-Net
#'
#' Adaptive MCP-Net
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
#' @param scale Scaling factor for adaptive weights:
#' \code{weights = coefficients^(-scale)}.
#' @param eps Convergence threshold to use in MCP-net.
#' @param max.iter Maximum number of iterations to use in MCP-net.
#' @param penalty.factor.init The multiplicative factor for the penalty
#' applied to each coefficient in the initial estimation step. This is
#' useful for incorporating prior information about variable weights,
#' for example, emphasizing specific clinical variables. To make certain
#' variables more likely to be selected, assign a smaller value.
#' Default is \code{rep(1, ncol(x))}.
#' @param seed Random seed for cross-validation fold division.
#' @param parallel Logical. Enable parallel parameter tuning or not,
#' default is \code{FALSE}. To enable parallel tuning, load the
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
#' @export amnet
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 150, p = 500, rho = 0.6,
#'   coef = rep(1, 5), snr = 2, p.train = 0.7,
#'   seed = 1001
#' )
#'
#' amnet.fit <- amnet(
#'   dat$x.tr, dat$y.tr,
#'   alphas = seq(0.2, 0.8, 0.2), seed = 1002
#' )
#'
#' print(amnet.fit)
#' msaenet.nzv(amnet.fit)
#' msaenet.fp(amnet.fit, 1:5)
#' msaenet.tp(amnet.fit, 1:5)
#' amnet.pred <- predict(amnet.fit, dat$x.te)
#' msaenet.rmse(dat$y.te, amnet.pred)
#' plot(amnet.fit)
amnet <- function(
    x, y,
    family = c("gaussian", "binomial", "poisson", "cox"),
    init = c("mnet", "ridge"),
    gammas = 3, alphas = seq(0.05, 0.95, 0.05),
    tune = c("cv", "ebic", "bic", "aic"),
    nfolds = 5L,
    ebic.gamma = 1,
    scale = 1,
    eps = 1e-4, max.iter = 10000L,
    penalty.factor.init = rep(1, ncol(x)),
    seed = 1001, parallel = FALSE, verbose = FALSE) {
  family <- match.arg(family)
  init <- match.arg(init)
  tune <- match.arg(tune)
  call <- match.call()
  nvar <- ncol(x)

  if (verbose) cat("Starting step 1 ...\n")

  if (init == "mnet") {
    mnet.cv <- msaenet.tune.ncvreg(
      x = x, y = y, family = family, penalty = "MCP",
      gammas = gammas, alphas = alphas,
      tune = tune,
      nfolds = nfolds,
      ebic.gamma = ebic.gamma,
      eps = eps, max.iter = max.iter,
      penalty.factor = penalty.factor.init,
      seed = seed, parallel = parallel
    )

    best.gamma.mnet <- mnet.cv$"best.gamma"
    best.alpha.mnet <- mnet.cv$"best.alpha"
    best.lambda.mnet <- mnet.cv$"best.lambda"
    step.criterion.mnet <- mnet.cv$"step.criterion"

    mnet.full <- .ncvnet(
      x = x, y = y, family = family, penalty = "MCP",
      gamma = best.gamma.mnet,
      alpha = best.alpha.mnet,
      lambda = best.lambda.mnet,
      eps = eps, max.iter = max.iter,
      penalty.factor = penalty.factor.init
    )

    bhat <- .coef.ncvreg(mnet.full, nvar)
  }

  if (init == "ridge") {
    mnet.cv <- msaenet.tune.glmnet(
      x = x, y = y, family = family,
      alphas = 0,
      tune = tune,
      nfolds = nfolds, rule = "lambda.min",
      ebic.gamma = ebic.gamma,
      lower.limits = -Inf, upper.limits = Inf,
      penalty.factor = penalty.factor.init,
      seed = seed, parallel = parallel
    )

    best.gamma.mnet <- NA
    best.alpha.mnet <- mnet.cv$"best.alpha"
    best.lambda.mnet <- mnet.cv$"best.lambda"
    step.criterion.mnet <- mnet.cv$"step.criterion"

    mnet.full <- glmnet(
      x = x, y = y, family = family,
      alpha = best.alpha.mnet,
      lambda = best.lambda.mnet,
      penalty.factor = penalty.factor.init
    )

    bhat <- as.matrix(mnet.full$"beta")
  }

  if (all(bhat == 0)) bhat <- rep(.Machine$double.eps * 2, length(bhat))

  adpen <- (pmax(abs(bhat), .Machine$double.eps))^(-scale)

  if (verbose) cat("Starting step 2 ...\n")

  amnet.cv <- msaenet.tune.ncvreg(
    x = x, y = y, family = family, penalty = "MCP",
    gammas = gammas, alphas = alphas,
    tune = tune,
    nfolds = nfolds,
    ebic.gamma = ebic.gamma,
    eps = eps, max.iter = max.iter,
    seed = seed + 1L, parallel = parallel,
    penalty.factor = adpen
  )

  best.gamma.amnet <- amnet.cv$"best.gamma"
  best.alpha.amnet <- amnet.cv$"best.alpha"
  best.lambda.amnet <- amnet.cv$"best.lambda"
  step.criterion.amnet <- amnet.cv$"step.criterion"

  amnet.full <- .ncvnet(
    x = x, y = y, family = family, penalty = "MCP",
    gamma = best.gamma.amnet,
    alpha = best.alpha.amnet,
    lambda = best.lambda.amnet,
    eps = eps, max.iter = max.iter,
    penalty.factor = adpen
  )

  # final beta stored as sparse matrix
  bhat.full <- Matrix(.coef.ncvreg(amnet.full, nvar), sparse = TRUE)
  bhat.first <- Matrix(.coef.ncvreg(mnet.full, nvar), sparse = TRUE)

  amnet.model <- list(
    "beta" = bhat.full,
    "model" = amnet.full,
    "beta.first" = bhat.first,
    "model.first" = mnet.full,
    "best.alpha.mnet" = best.alpha.mnet,
    "best.alpha.amnet" = best.alpha.amnet,
    "best.lambda.mnet" = best.lambda.mnet,
    "best.lambda.amnet" = best.lambda.amnet,
    "best.gamma.mnet" = best.gamma.mnet,
    "best.gamma.amnet" = best.gamma.amnet,
    "step.criterion" = c(step.criterion.mnet, step.criterion.amnet),
    "adpen" = adpen,
    "seed" = seed,
    "call" = call
  )

  class(amnet.model) <- c("msaenet", "msaenet.amnet")
  amnet.model
}
