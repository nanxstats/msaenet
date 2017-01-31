#' Automatic parameter tuning for glmnet models
#'
#' @return best model object and best alpha value
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom glmnet cv.glmnet
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#'
#' @keywords internal

msaenet.tune.glmnet = function(x, y, family,
                               alphas,
                               nfolds, rule,
                               seed, parallel, ...) {

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

  # select model for best alpha (then lambda)
  # criterion: minimal cross-validation error
  errors = unlist(lapply(model.list, function(x) min(sqrt(x$'cvm'))))
  best.model = model.list[[which.min(errors)]]

  best.alpha = alphas[which.min(errors)]

  if (rule == 'lambda.min') best.lambda = best.model$'lambda.min'
  if (rule == 'lambda.1se') best.lambda = best.model$'lambda.1se'

  list('best.model'  = best.model,
       'best.alpha'  = best.alpha,
       'best.lambda' = best.lambda)

}
