#' Automatic parameter tuning for ncvreg by k-fold cross-validation
#'
#' @return best model object, best gamma, and best alpha
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom ncvreg cv.ncvreg
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' @importFrom foreach foreach
#'
#' @keywords internal

msaenet.tune.ncvreg = function(x, y, family, penalty,
                               gammas, alphas,
                               nfolds,
                               eps, max.iter,
                               seed, parallel, ...) {

  if (!parallel) {

    model.list = vector('list', length(gammas))
    for (k in 1L:length(model.list)) {
      model.list[[k]] = vector('list', length(alphas))
    }

    for (i in 1L:length(gammas)) {
      for (j in 1L:length(alphas)) {
        set.seed(seed)
        if (family == "cox") {
          model.list[[i]][[j]] =
            cv.ncvsurv(X = x, y = y, penalty = penalty,
                       gamma = gammas[i], alpha = alphas[j],
                       nfolds = nfolds,
                       eps = eps, max.iter = max.iter, ...)
        } else {
          model.list[[i]][[j]] =
            cv.ncvreg(X = x, y = y, family = family, penalty = penalty,
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
          cv.ncvsurv(X = x, y = y, penalty = penalty,
                     gamma = gammas, alpha = alphas,
                     nfolds = nfolds,
                     eps = eps, max.iter = max.iter, ...)
        }
    } else {
      model.list <- foreach(gammas = gammas) %:%
        foreach(alphas = alphas) %dopar% {
          set.seed(seed)
          cv.ncvreg(X = x, y = y, family = family, penalty = penalty,
                    gamma = gammas, alpha = alphas,
                    nfolds = nfolds,
                    eps = eps, max.iter = max.iter, ...)
        }
    }

  }

  simple.model.list = unlist(model.list, recursive = FALSE)

  # select model for best gamma/alpha first (then lambda)
  # criterion: minimal cross-validation error
  errors = unlist(lapply(simple.model.list, function(x) min(sqrt(x$'cve'))))
  best.model = simple.model.list[[which.min(errors)]]

  list('best.model'  = best.model,
       'best.gamma'  = best.model$'fit'$'gamma',
       'best.alpha'  = best.model$'fit'$'alpha',
       'best.lambda' = best.model$'lambda.min')

}

# wrapper for ncvreg::ncvreg and ncvreg::ncvsurv with two hotfixes
.ncvnet = function (x, y, family, penalty,
                    gamma, alpha, lambda,
                    eps, max.iter, ...) {

  if (family == 'cox') {
    fit = ncvreg::ncvsurv(X = x, y = y, penalty = penalty,
                          gamma = gamma, alpha = alpha, lambda = lambda,
                          eps = eps, max.iter = max.iter, ...)
  } else {
    fit = ncvreg::ncvreg(X = x, y = y, family = family, penalty = penalty,
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
