#' Automatic parameter tuning for glmnet by k-fold cross-validation
#'
#' @return best model object and best alpha value
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal

msaenet.tune.glmnet = function(..., alphas, seed, parallel) {

  if (!parallel) {
    model.list = vector('list', length(alphas))
    for (i in 1L:length(alphas)) {
      set.seed(seed)
      model.list[[i]] = cv.glmnet(..., alpha = alphas[i])
    }
  } else {
    model.list <- foreach(alphas = alphas) %dopar% {
      set.seed(seed)
      cv.glmnet(..., alpha = alphas)
    }
  }

  # select model for best lambda first (then alpha)
  # criterion: minimal cross-validation error
  errors = unlist(lapply(model.list, function(x) min(sqrt(x$'cvm'))))

  return(list('best.model' = model.list[[which.min(errors)]],
              'best.alpha' = alphas[which.min(errors)]))

}
