% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msaenet-inspect.R
\name{msaenet.tp}
\alias{msaenet.tp}
\title{Get the Number of True Positive Selections}
\usage{
msaenet.tp(object, true.idx)
}
\arguments{
\item{object}{An object of class \code{msaenet} produced
by \code{\link{aenet}}, \code{amnet}, \code{asnet},
\code{\link{msaenet}}, \code{\link{msamnet}}, or \code{\link{msasnet}}.}

\item{true.idx}{Vector. Indices of true variables.}
}
\value{
Number of true positive variables in the model.
}
\description{
Get the number of true positive selections from msaenet model objects,
given the indices of true variables (if known).
}
\examples{
dat <- msaenet.sim.gaussian(
  n = 150, p = 500, rho = 0.6,
  coef = rep(1, 5), snr = 2, p.train = 0.7,
  seed = 1001
)

msaenet.fit <- msaenet(
  dat$x.tr, dat$y.tr,
  alphas = seq(0.2, 0.8, 0.2),
  nsteps = 3L, seed = 1003
)

msaenet.tp(msaenet.fit, 1:5)
}
\author{
Nan Xiao <\url{https://nanx.me}>
}
