#' Automatic parameter tuning for ncvreg by k-fold cross-validation
#'
#' @return best model object, best gamma, and best alpha
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @importFrom ncvreg cv.ncvreg
#' @importFrom ncvreg cv.ncvsurv
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' @importFrom foreach foreach
#'
#' @keywords internal

msaenet.tune.ncvreg = function(..., family, gammas, alphas,
                               eps, max.iter, seed, parallel) {

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
            cv.ncvsurv(..., gamma = gammas[i], alpha = alphas[j],
                       eps = eps, max.iter = max.iter)
        } else {
          model.list[[i]][[j]] =
            cv.ncvreg(..., family = family, gamma = gammas[i], alpha = alphas[j],
                      eps = eps, max.iter = max.iter)
        }
      }
    }

  } else {

    if (family == "cox") {
      model.list <- foreach(gammas = gammas) %:%
        foreach(alphas = alphas) %dopar% {
          set.seed(seed)
          cv.ncvsurv(..., gamma = gammas, alpha = alphas,
                     eps = eps, max.iter = max.iter)
        }
    } else {
      model.list <- foreach(gammas = gammas) %:%
        foreach(alphas = alphas) %dopar% {
          set.seed(seed)
          cv.ncvreg(..., family = family, gamma = gammas, alpha = alphas,
                    eps = eps, max.iter = max.iter)
        }
    }

  }

  simple.model.list = unlist(model.list, recursive = FALSE)

  # select model for best lambda first (then gamma/alpha)
  # criterion: minimal cross-validation error
  errors = unlist(lapply(simple.model.list, function(x) min(sqrt(x$'cve'))))

  return(list('best.model' = simple.model.list[[which.min(errors)]],
              'best.gamma' = simple.model.list[[which.min(errors)]]$'fit'$'gamma',
              'best.alpha' = simple.model.list[[which.min(errors)]]$'fit'$'alpha'))

}

# wrapper for ncvreg::ncvreg and ncvreg::ncvsurv with two hotfixes
.ncvnet = function (..., lambda, family, eps, max.iter) {

  if (family == 'cox') {
    fit = ncvreg::ncvsurv(..., lambda = lambda,
                          eps = eps, max.iter = max.iter)
  } else {
    fit = ncvreg::ncvreg(..., family = family, lambda = lambda,
                         eps = eps, max.iter = max.iter)
  }

  fit

}

# returns "df" for ncvreg model objects
.ncvdf = function(model)
  sum(abs(as.vector(model[['beta']])[-1L]) > .Machine$double.eps)

# check if ncvreg model object has an intercept
# return a coef vector without intercept
.ncv.coef = function(model, nvar) {

  nvar.model = length(as.vector(model$'beta'))

  bhat = if (nvar.model == nvar + 1L)
    as.vector(model$'beta')[-1L] else
      as.vector(model$'beta')

  bhat

}
