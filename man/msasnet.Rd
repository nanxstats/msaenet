% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msasnet.R
\name{msasnet}
\alias{msasnet}
\title{Multi-Step Adaptive SCAD-Net}
\usage{
msasnet(
  x,
  y,
  family = c("gaussian", "binomial", "poisson", "cox"),
  init = c("snet", "ridge"),
  gammas = 3.7,
  alphas = seq(0.05, 0.95, 0.05),
  tune = c("cv", "ebic", "bic", "aic"),
  nfolds = 5L,
  ebic.gamma = 1,
  nsteps = 2L,
  tune.nsteps = c("max", "ebic", "bic", "aic"),
  ebic.gamma.nsteps = 1,
  scale = 1,
  eps = 1e-04,
  max.iter = 10000L,
  penalty.factor.init = rep(1, ncol(x)),
  seed = 1001,
  parallel = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{x}{Data matrix.}

\item{y}{Response vector if \code{family} is \code{"gaussian"},
\code{"binomial"}, or \code{"poisson"}. If \code{family} is
\code{"cox"}, a response matrix created by \code{\link[survival]{Surv}}.}

\item{family}{Model family, can be \code{"gaussian"},
\code{"binomial"}, \code{"poisson"}, or \code{"cox"}.}

\item{init}{Type of the penalty used in the initial
estimation step. Can be \code{"snet"} or \code{"ridge"}.}

\item{gammas}{Vector of candidate \code{gamma}s (the concavity parameter)
to use in SCAD-Net. Default is \code{3.7}.}

\item{alphas}{Vector of candidate \code{alpha}s to use in SCAD-Net.}

\item{tune}{Parameter tuning method for each estimation step.
Possible options are \code{"cv"}, \code{"ebic"}, \code{"bic"},
and \code{"aic"}. Default is \code{"cv"}.}

\item{nfolds}{Fold numbers of cross-validation when \code{tune = "cv"}.}

\item{ebic.gamma}{Parameter for Extended BIC penalizing
size of the model space when \code{tune = "ebic"},
default is \code{1}. For details, see Chen and Chen (2008).}

\item{nsteps}{Maximum number of adaptive estimation steps.
At least \code{2}, assuming adaptive SCAD-net has only
one adaptive estimation step.}

\item{tune.nsteps}{Optimal step number selection method
(aggregate the optimal model from the each step and compare).
Options include \code{"max"} (select the final-step model directly),
or compare these models using \code{"ebic"}, \code{"bic"}, or \code{"aic"}.
Default is \code{"max"}.}

\item{ebic.gamma.nsteps}{Parameter for Extended BIC penalizing
size of the model space when \code{tune.nsteps = "ebic"},
default is \code{1}.}

\item{scale}{Scaling factor for adaptive weights:
\code{weights = coefficients^(-scale)}.}

\item{eps}{Convergence threshold to use in SCAD-net.}

\item{max.iter}{Maximum number of iterations to use in SCAD-net.}

\item{penalty.factor.init}{The multiplicative factor for the penalty
applied to each coefficient in the initial estimation step. This is
useful for incorporating prior information about variable weights,
for example, emphasizing specific clinical variables. To make certain
variables more likely to be selected, assign a smaller value.
Default is \code{rep(1, ncol(x))}.}

\item{seed}{Random seed for cross-validation fold division.}

\item{parallel}{Logical. Enable parallel parameter tuning or not,
default is \code{FALSE}. To enable parallel tuning, load the
\code{doParallel} package and run \code{registerDoParallel()}
with the number of CPU cores before calling this function.}

\item{verbose}{Should we print out the estimation progress?}
}
\value{
List of model coefficients, \code{ncvreg} model object,
and the optimal parameter set.
}
\description{
Multi-Step Adaptive SCAD-Net
}
\examples{
dat <- msaenet.sim.gaussian(
  n = 150, p = 500, rho = 0.6,
  coef = rep(1, 5), snr = 2, p.train = 0.7,
  seed = 1001
)

msasnet.fit <- msasnet(
  dat$x.tr, dat$y.tr,
  alphas = seq(0.3, 0.9, 0.3),
  nsteps = 3L, seed = 1003
)

print(msasnet.fit)
msaenet.nzv(msasnet.fit)
msaenet.fp(msasnet.fit, 1:5)
msaenet.tp(msasnet.fit, 1:5)
msasnet.pred <- predict(msasnet.fit, dat$x.te)
msaenet.rmse(dat$y.te, msasnet.pred)
plot(msasnet.fit)
}
\author{
Nan Xiao <\url{https://nanx.me}>
}
