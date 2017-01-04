.is.msaenet = function(x) 'msaenet' %in% class(x)

.is.glmnet = function(x) 'glmnet' %in% class(x)

.is.ncvreg = function(x) 'ncvreg' %in% class(x)

.is.adaptive = function(x)
  any(class(x) %in% c('msaenet.aenet', 'msaenet.amnet', 'msaenet.asnet'))

.is.multistep = function(x)
  any(class(x) %in% c('msaenet.msaenet', 'msaenet.msamnet', 'msaenet.msasnet'))
