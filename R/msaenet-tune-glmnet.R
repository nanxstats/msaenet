#' Automatic (parallel) parameter tuning for glmnet models
#'
#' @return Optimal model object, parameter set, and criterion value
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @references
#' Chen, Jiahua, and Zehua Chen. (2008).
#' Extended Bayesian information criteria for model selection with
#' large model spaces. \emph{Biometrika} 95(3), 759--771.
#'
#' @keywords internal

msaenet.tune.glmnet = function(x, y, family,
                               alphas,
                               tune,
                               nfolds, rule,
                               ebic.gamma,
                               seed, parallel, ...) {

  if (tune == 'cv') {

    if (!parallel) {
      model.list = vector('list', length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] = cv.glmnet(x = x, y = y, family = family,
                                    nfolds = nfolds, alpha = alphas[i], ...)
      }
    } else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        cv.glmnet(x = x, y = y, family = family,
                  nfolds = nfolds, alpha = alphas, ...)
      }
    }

    errors = unlist(lapply(model.list, function(x) min(sqrt(x$'cvm'))))
    errors.min.idx = which.min(errors)
    best.model = model.list[[errors.min.idx]]

    best.alpha = alphas[errors.min.idx]

    if (rule == 'lambda.min') best.lambda = best.model$'lambda.min'
    if (rule == 'lambda.1se') best.lambda = best.model$'lambda.1se'

    step.criterion = errors[errors.min.idx]

  } else {

    if (!parallel) {
      model.list = vector('list', length(alphas))
      for (i in 1L:length(alphas)) {
        set.seed(seed)
        model.list[[i]] = glmnet(x = x, y = y, family = family,
                                 alpha = alphas[i], ...)
      }
    } else {
      model.list <- foreach(alphas = alphas) %dopar% {
        set.seed(seed)
        glmnet(x = x, y = y, family = family,
               alpha = alphas, ...)
      }
    }

    if (tune == 'aic') {

      ics.list = mapply(
        .aic,
        deviance = lapply(model.list, .deviance),
        df       = lapply(model.list, .df),
        SIMPLIFY = FALSE)

    }

    if (tune == 'bic') {

      ics.list = mapply(
        .bic,
        deviance = lapply(model.list, .deviance),
        df       = lapply(model.list, .df),
        nobs     = lapply(model.list, .nobs),
        SIMPLIFY = FALSE)

    }

    if (tune == 'ebic') {

      ics.list = mapply(
        .ebic,
        deviance = lapply(model.list, .deviance),
        df       = lapply(model.list, .df),
        nobs     = lapply(model.list, .nobs),
        nvar     = lapply(model.list, .nvar),
        gamma    = ebic.gamma,
        SIMPLIFY = FALSE)

    }

    ics = sapply(ics.list, function(x) min(x))
    ics.min.idx = which.min(ics)
    best.model  = model.list[[ics.min.idx]]

    best.alpha = alphas[ics.min.idx]

    best.ic.min.idx = which.min(ics.list[[ics.min.idx]])
    best.lambda = best.model$'lambda'[[best.ic.min.idx]]

    step.criterion = ics.list[[ics.min.idx]][[best.ic.min.idx]]

  }

  list('best.model'     = best.model,
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

msaenet.tune.nsteps.glmnet = function(model.list,
                                      tune.nsteps, ebic.gamma.nsteps) {

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
