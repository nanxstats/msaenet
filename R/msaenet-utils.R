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

# deviance vector for glmnet model objects
.deviance.glmnet = function(model) (1L - model$'dev.ratio') * model$'nulldev'

# degree of freedom vector for glmnet model objects
.df.glmnet = function(model) model$'df'

# number of observations in X from glmnet model objects
.nobs.glmnet = function(model) model$'nobs'

# dimensionality of X from glmnet model objects
.nvar.glmnet = function(model) model$'dim'[[1L]]

# deviance vector for ncvreg model objects
.deviance.ncvreg = function(model) model$'loss'

# degree of freedom vector for ncvreg model objects
.df.ncvreg = function(model)
  unname(colSums(as.matrix(abs(model$'beta'[-1L, ])) > .Machine$double.eps))

# number of observations in X from ncvreg model objects
.nobs.ncvreg = function(model) model$'n'

# dimensionality of X from ncvreg model objects
.nvar.ncvreg = function(model) length(model$'penalty.factor')
