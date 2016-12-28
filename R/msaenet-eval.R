#' Root Mean Squared Error (RMSE)
#'
#' Compute root mean squared error (RMSE).
#'
#' @param yreal Vector. True response.
#' @param ypred Vector. Predicted response.
#'
#' @return RMSE
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export msaenet.rmse
msaenet.rmse = function(yreal, ypred) {

  rmse = sqrt(mean((yreal - ypred)^2))
  return(rmse)

}

#' Mean Absolute Error (MAE)
#'
#' Compute mean absolute error (MAE).
#'
#' @param yreal Vector. True response.
#' @param ypred Vector. Predicted response.
#'
#' @return MAE
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export msaenet.mae
msaenet.mae = function(yreal, ypred) {

  mae = mean(abs(yreal - ypred))
  return(mae)

}

#' Root Mean Squared Logarithmic Error (RMSLE)
#'
#' Compute root mean squared logarithmic error (RMSLE).
#'
#' @param yreal Vector. True response.
#' @param ypred Vector. Predicted response.
#'
#' @return RMSLE
#'
#' @author Nan Xiao <\url{http://nanx.me}>
#'
#' @export msaenet.rmsle
msaenet.rmsle = function(yreal, ypred) {

  rmsle = sqrt(mean((log(ypred + 1) - log(yreal + 1))^2))
  return(rmsle)

}
