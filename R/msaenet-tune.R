#' Automatic alpha tuning by k-fold cross-validation
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

msaenet.tune.glmnet.alpha = function(..., alphas, seed, parallel) {

  if (!parallel) {
    modelList = vector('list', length(alphas))
    for (i in 1L:length(alphas)) {
      set.seed(seed)
      modelList[[i]] = cv.glmnet(..., alpha = alphas[i])
    }
  } else {
    modelList <- foreach(alphas = alphas) %dopar% {
      set.seed(seed)
      cv.glmnet(..., alpha = alphas)
    }
  }

  # Choose model for best lambda first (then alpha)
  errors = unlist(lapply(modelList, function(x) min(sqrt(x$cvm))))

  return(list('best.model' = modelList[[which.min(errors)]],
              'best.alpha' = alphas[which.min(errors)]))

}
