#' Automatic (parallel) parameter tuning for ncvreg models
#'
#' @return Optimal model object, parameter set, and criterion value
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom ncvreg cv.ncvreg
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' @importFrom foreach foreach
#'
#' @references
#' Chen, Jiahua, and Zehua Chen. (2008).
#' Extended Bayesian information criteria for model selection with
#' large model spaces. \emph{Biometrika} 95(3), 759--771.
#'
#' @keywords internal

msaenet.tune.ncvreg = function(
  x, y, family, penalty,
  gammas, alphas,
  tune,
  nfolds,
  ebic.gamma,
  eps, max.iter,
  seed, parallel, ...) {

  if (tune == 'cv') {

    if (!parallel) {

      model.list = vector('list', length(gammas))
      for (k in 1L:length(model.list)) {
        model.list[[k]] = vector('list', length(alphas))
      }

      for (i in 1L:length(gammas)) {
        for (j in 1L:length(alphas)) {
          set.seed(seed)
          if (family == "cox") {
            model.list[[i]][[j]] = cv.ncvsurv(
              X = x, y = y, penalty = penalty,
              gamma = gammas[i], alpha = alphas[j],
              nfolds = nfolds,
              eps = eps, max.iter = max.iter, ...)
          } else {
            model.list[[i]][[j]] = cv.ncvreg(
              X = x, y = y, family = family, penalty = penalty,
              gamma = gammas[i], alpha = alphas[j],
              nfolds = nfolds,
              eps = eps, max.iter = max.iter, ...)
          }
        }
      }

    } else {

      if (family == "cox") {
        model.list <- foreach(gammas = gammas) %:%
          foreach(alphas = alphas) %dopar% {
            set.seed(seed)
            cv.ncvsurv(
              X = x, y = y, penalty = penalty,
              gamma = gammas, alpha = alphas,
              nfolds = nfolds,
              eps = eps, max.iter = max.iter, ...)
          }
      } else {
        model.list <- foreach(gammas = gammas) %:%
          foreach(alphas = alphas) %dopar% {
            set.seed(seed)
            cv.ncvreg(
              X = x, y = y, family = family, penalty = penalty,
              gamma = gammas, alpha = alphas,
              nfolds = nfolds,
              eps = eps, max.iter = max.iter, ...)
          }
      }

    }

    simple.model.list = unlist(model.list, recursive = FALSE)

    errors = unlist(lapply(simple.model.list, function(x) min(sqrt(x$'cve'))))
    errors.min.idx = which.min(errors)
    best.model = simple.model.list[[errors.min.idx]]

    best.gamma  = best.model$'fit'$'gamma'
    best.alpha  = best.model$'fit'$'alpha'
    best.lambda = best.model$'lambda.min'

    step.criterion = errors[errors.min.idx]

  } else {

    if (!parallel) {

      model.list = vector('list', length(gammas))
      for (k in 1L:length(model.list)) {
        model.list[[k]] = vector('list', length(alphas))
      }

      for (i in 1L:length(gammas)) {
        for (j in 1L:length(alphas)) {
          set.seed(seed)
          if (family == "cox") {
            model.list[[i]][[j]] = ncvsurv(
              X = x, y = y, penalty = penalty,
              gamma = gammas[i], alpha = alphas[j],
              eps = eps, max.iter = max.iter, ...)
          } else {
            model.list[[i]][[j]] = ncvreg(
              X = x, y = y, family = family, penalty = penalty,
              gamma = gammas[i], alpha = alphas[j],
              eps = eps, max.iter = max.iter, ...)
          }
        }
      }

    } else {

      if (family == "cox") {
        model.list <- foreach(gammas = gammas) %:%
          foreach(alphas = alphas) %dopar% {
            set.seed(seed)
            ncvsurv(
              X = x, y = y, penalty = penalty,
              gamma = gammas, alpha = alphas,
              eps = eps, max.iter = max.iter, ...)
          }
      } else {
        model.list <- foreach(gammas = gammas) %:%
          foreach(alphas = alphas) %dopar% {
            set.seed(seed)
            ncvreg(
              X = x, y = y, family = family, penalty = penalty,
              gamma = gammas, alpha = alphas,
              eps = eps, max.iter = max.iter, ...)
          }
      }

    }

    simple.model.list = unlist(model.list, recursive = FALSE)

    if (tune == 'aic') {

      ics.list = mapply(
        .aic,
        deviance = lapply(simple.model.list, .deviance),
        df       = lapply(simple.model.list, .df),
        SIMPLIFY = FALSE)

    }

    if (tune == 'bic') {

      ics.list = mapply(
        .bic,
        deviance = lapply(simple.model.list, .deviance),
        df       = lapply(simple.model.list, .df),
        nobs     = lapply(simple.model.list, .nobs),
        SIMPLIFY = FALSE)

    }

    if (tune == 'ebic') {

      ics.list = mapply(
        .ebic,
        deviance = lapply(simple.model.list, .deviance),
        df       = lapply(simple.model.list, .df),
        nobs     = lapply(simple.model.list, .nobs),
        nvar     = lapply(simple.model.list, .nvar),
        gamma    = ebic.gamma,
        SIMPLIFY = FALSE)

    }

    ics = sapply(ics.list, function(x) min(x))
    ics.min.idx = which.min(ics)
    best.model  = simple.model.list[[ics.min.idx]]

    best.gamma  = best.model$'gamma'
    best.alpha  = best.model$'alpha'

    best.ic.min.idx = which.min(ics.list[[ics.min.idx]])
    best.lambda = best.model$'lambda'[[best.ic.min.idx]]

    step.criterion = ics.list[[ics.min.idx]][[best.ic.min.idx]]

  }

  list('best.model'     = best.model,
       'best.gamma'     = best.gamma,
       'best.alpha'     = best.alpha,
       'best.lambda'    = best.lambda,
       'step.criterion' = step.criterion)

}

#' Select the number of adaptive estimation steps
#'
#' @return optimal step number
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @keywords internal

msaenet.tune.nsteps.ncvreg = function(
  model.list, tune.nsteps, ebic.gamma.nsteps) {

  nmods = length(model.list)

  if (tune.nsteps == 'max') {

    ics = NULL
    best.step = nmods

  } else {

    if (tune.nsteps == 'aic')
      ics = .aic(
        deviance  = sapply(model.list, .deviance),
        df        = sapply(model.list, .df))

    if (tune.nsteps == 'bic')
      ics = .bic(
        deviance  = sapply(model.list, .deviance),
        df        = sapply(model.list, .df),
        nobs      = sapply(model.list, .nobs))

    if (tune.nsteps == 'ebic')
      ics = .ebic(
        deviance = sapply(model.list, .deviance),
        df       = sapply(model.list, .df),
        nobs     = sapply(model.list, .nobs),
        nvar     = sapply(model.list, .nvar),
        gamma    = ebic.gamma.nsteps)

    best.step = which.min(ics)

  }

  list('best.step' = best.step, 'ics' = ics)

}

# wrapper for ncvreg::ncvreg and ncvreg::ncvsurv with two hotfixes
.ncvnet = function(
  x, y, family, penalty,
  gamma, alpha, lambda,
  eps, max.iter, ...) {

  if (family == 'cox') {
    fit = ncvreg::ncvsurv(
      X = x, y = y, penalty = penalty,
      gamma = gamma, alpha = alpha, lambda = lambda,
      eps = eps, max.iter = max.iter, ...)
  } else {
    fit = ncvreg::ncvreg(
      X = x, y = y, family = family, penalty = penalty,
      gamma = gamma, alpha = alpha, lambda = lambda,
      eps = eps, max.iter = max.iter, ...)
  }

  fit

}

# check if ncvreg model object has an intercept
# return a coef vector without intercept
.coef.ncvreg = function(model, nvar) {

  nvar.model = length(as.vector(model$'beta'))

  bhat = if (nvar.model == nvar + 1L)
    as.vector(model$'beta')[-1L] else
      as.vector(model$'beta')

  bhat

}
