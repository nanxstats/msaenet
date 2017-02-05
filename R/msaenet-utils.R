.is.msaenet = function(x) 'msaenet' %in% class(x)

.is.glmnet = function(x) 'glmnet' %in% class(x)

.is.ncvreg = function(x) 'ncvreg' %in% class(x)

.is.adaptive = function(x)
  any(class(x) %in% c('msaenet.aenet', 'msaenet.amnet', 'msaenet.asnet'))

.is.multistep = function(x)
  any(class(x) %in% c('msaenet.msaenet', 'msaenet.msamnet', 'msaenet.msasnet'))

# AIC
.aic = function(deviance, df)
  deviance + (2L * df)

# BIC
.bic = function(deviance, df, nobs)
  deviance + (log(nobs) * df)

# Extended BIC
.ebic = function(deviance, df, nobs, nvar, gamma)
  deviance + (log(nobs) * df) + (2L * gamma * lchoose(nvar, df))

# deviance vector
.deviance = function(model) {
  if (.is.glmnet(model)) return((1L - model$'dev.ratio') * model$'nulldev')
  if (.is.ncvreg(model)) return(model$'loss')
}

# degree of freedom vector
.df = function(model) {
  if (.is.glmnet(model)) return(model$'df')
  if (.is.ncvreg(model)) return(unname(colSums(
    as.matrix(abs(model$'beta'[-1L, ])) > .Machine$double.eps)))
}

# number of observations in predictor matrix
.nobs = function(model) {
  if (.is.glmnet(model)) return(model$'nobs')
  if (.is.ncvreg(model)) return(model$'n')
}

# dimensionality of predictor matrix
.nvar = function(model) {
  if (.is.glmnet(model)) return(model$'dim'[[1L]])
  if (.is.ncvreg(model)) return(length(model$'penalty.factor'))
}
