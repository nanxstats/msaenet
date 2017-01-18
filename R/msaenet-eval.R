#' Mean Squared Error (MSE)
#'
#' Compute mean squared error (MSE).
#'
#' @param yreal Vector. True response.
#' @param ypred Vector. Predicted response.
#'
#' @return MSE
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.mse
msaenet.mse = function(yreal, ypred)
  mean((yreal - ypred)^2)

#' Root Mean Squared Error (RMSE)
#'
#' Compute root mean squared error (RMSE).
#'
#' @param yreal Vector. True response.
#' @param ypred Vector. Predicted response.
#'
#' @return RMSE
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.rmse
msaenet.rmse = function(yreal, ypred)
  sqrt(mean((yreal - ypred)^2))

#' Mean Absolute Error (MAE)
#'
#' Compute mean absolute error (MAE).
#'
#' @param yreal Vector. True response.
#' @param ypred Vector. Predicted response.
#'
#' @return MAE
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.mae
msaenet.mae = function(yreal, ypred)
  mean(abs(yreal - ypred))

#' Root Mean Squared Logarithmic Error (RMSLE)
#'
#' Compute root mean squared logarithmic error (RMSLE).
#'
#' @param yreal Vector. True response.
#' @param ypred Vector. Predicted response.
#'
#' @return RMSLE
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @export msaenet.rmsle
msaenet.rmsle = function(yreal, ypred)
  sqrt(mean((log(ypred + 1) - log(yreal + 1))^2))
