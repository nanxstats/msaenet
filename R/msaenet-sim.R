#' Generate Simulation Data for Benchmarking Sparse Regressions
#' (Gaussian Response)
#'
#' Generate simulation data (Gaussian case) following the
#' settings in Xiao and Xu (2015).
#'
#' @param n Number of observations.
#' @param p Number of variables.
#' @param rho Correlation base for generating correlated variables.
#' @param coef Vector of non-zero coefficients.
#' @param snr Signal-to-noise ratio (SNR).
#' @param p.train Percentage of training set.
#' @param seed Random seed for reproducibility.
#'
#' @return List of \code{x.tr}, \code{x.te}, \code{y.tr}, and \code{y.te}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @references
#' Nan Xiao and Qing-Song Xu. (2015). Multi-step adaptive elastic-net:
#' reducing false positives in high-dimensional variable selection.
#' \emph{Journal of Statistical Computation and Simulation} 85(18), 3755--3765.
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#'
#' @export msaenet.sim.gaussian
#'
#' @examples
#' dat <- msaenet.sim.gaussian(
#'   n = 300, p = 500, rho = 0.6,
#'   coef = rep(1, 10), snr = 3, p.train = 0.7,
#'   seed = 1001
#' )
#'
#' dim(dat$x.tr)
#' dim(dat$x.te)
msaenet.sim.gaussian <- function(
  n = 300, p = 500,
  rho = 0.5, coef = rep(0.2, 50), snr = 1,
  p.train = 0.7, seed = 1001) {

  call <- match.call()

  set.seed(seed)

  sigma <- matrix(0, p, p)
  corvec <- function(i, p, rho) rho^(abs(i - 1L:p))
  for (i in 1:p) sigma[i, ] <- corvec(i, p, rho)

  X <- rmvnorm(n, rep(0, p), sigma)

  # non-zero coefficients
  beta0 <- matrix(c(coef, rep(0, (p - length(coef)))))

  snr.numerator <- as.vector(t(beta0) %*% sigma %*% beta0)
  snr.denominator <- snr.numerator / snr
  sd <- sqrt(snr.denominator)
  eps <- matrix(rnorm(n, 0, sd))

  y <- as.matrix((X %*% beta0) + eps)

  # training / test set splitting
  tr.row <- sample(1L:n, round(n * p.train), replace = FALSE)

  x.tr <- X[tr.row, , drop = FALSE]
  y.tr <- y[tr.row, , drop = FALSE]
  x.te <- X[-tr.row, , drop = FALSE]
  y.te <- y[-tr.row, , drop = FALSE]

  list(
    "x.tr" = x.tr, "y.tr" = y.tr,
    "x.te" = x.te, "y.te" = y.te,
    "call" = call
  )
}

#' Generate Simulation Data for Benchmarking Sparse Regressions
#' (Binomial Response)
#'
#' Generate simulation data for benchmarking sparse logistic regression models.
#'
#' @param n Number of observations.
#' @param p Number of variables.
#' @param rho Correlation base for generating correlated variables.
#' @param coef Vector of non-zero coefficients.
#' @param snr Signal-to-noise ratio (SNR).
#' @param p.train Percentage of training set.
#' @param seed Random seed for reproducibility.
#'
#' @return List of \code{x.tr}, \code{x.te}, \code{y.tr}, and \code{y.te}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @importFrom stats rbinom
#'
#' @export msaenet.sim.binomial
#'
#' @examples
#' dat <- msaenet.sim.binomial(
#'   n = 300, p = 500, rho = 0.6,
#'   coef = rep(1, 10), snr = 3, p.train = 0.7,
#'   seed = 1001
#' )
#'
#' dim(dat$x.tr)
#' dim(dat$x.te)
#' table(dat$y.tr)
#' table(dat$y.te)
msaenet.sim.binomial <- function(
  n = 300, p = 500,
  rho = 0.5, coef = rep(0.2, 50), snr = 1,
  p.train = 0.7, seed = 1001) {

  call <- match.call()

  set.seed(seed)

  sigma <- matrix(0, p, p)
  corvec <- function(i, p, rho) rho^(abs(i - 1L:p))
  for (i in 1:p) sigma[i, ] <- corvec(i, p, rho)

  X <- rmvnorm(n, rep(0, p), sigma)

  # non-zero coefficients
  beta0 <- matrix(c(coef, rep(0, (p - length(coef)))))

  snr.numerator <- as.vector(t(beta0) %*% sigma %*% beta0)
  snr.denominator <- snr.numerator / snr
  sd <- sqrt(snr.denominator)
  eps <- matrix(rnorm(n, 0, sd))

  f.eta <- (X %*% beta0) + eps
  f.prob <- exp(f.eta) / (1L + exp(f.eta))
  y <- as.matrix(rbinom(n, 1L, f.prob))

  # training / test set splitting
  tr.row <- sample(1L:n, round(n * p.train), replace = FALSE)

  x.tr <- X[tr.row, , drop = FALSE]
  y.tr <- y[tr.row, , drop = FALSE]
  x.te <- X[-tr.row, , drop = FALSE]
  y.te <- y[-tr.row, , drop = FALSE]

  list(
    "x.tr" = x.tr, "y.tr" = y.tr,
    "x.te" = x.te, "y.te" = y.te,
    "call" = call
  )
}

#' Generate Simulation Data for Benchmarking Sparse Regressions
#' (Poisson Response)
#'
#' Generate simulation data for benchmarking sparse Poisson regression models.
#'
#' @param n Number of observations.
#' @param p Number of variables.
#' @param rho Correlation base for generating correlated variables.
#' @param coef Vector of non-zero coefficients.
#' @param snr Signal-to-noise ratio (SNR).
#' @param p.train Percentage of training set.
#' @param seed Random seed for reproducibility.
#'
#' @return List of \code{x.tr}, \code{x.te}, \code{y.tr}, and \code{y.te}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @importFrom stats rpois
#'
#' @export msaenet.sim.poisson
#'
#' @examples
#' dat <- msaenet.sim.poisson(
#'   n = 300, p = 500, rho = 0.6,
#'   coef = rep(1, 10), snr = 3, p.train = 0.7,
#'   seed = 1001
#' )
#'
#' dim(dat$x.tr)
#' dim(dat$x.te)
msaenet.sim.poisson <- function(
  n = 300, p = 500,
  rho = 0.5, coef = rep(0.2, 50), snr = 1,
  p.train = 0.7, seed = 1001) {

  call <- match.call()

  set.seed(seed)

  sigma <- matrix(0, p, p)
  corvec <- function(i, p, rho) rho^(abs(i - 1L:p))
  for (i in 1:p) sigma[i, ] <- corvec(i, p, rho)

  X <- rmvnorm(n, rep(0, p), sigma)

  # non-zero coefficients
  beta0 <- matrix(c(coef, rep(0, (p - length(coef)))))

  snr.numerator <- as.vector(t(beta0) %*% sigma %*% beta0)
  snr.denominator <- snr.numerator / snr
  sd <- sqrt(snr.denominator)
  eps <- matrix(rnorm(n, 0, sd))

  rates <- exp((X %*% beta0) + eps)
  y <- as.matrix(rpois(n, rates))

  # training / test set splitting
  tr.row <- sample(1L:n, round(n * p.train), replace = FALSE)

  x.tr <- X[tr.row, , drop = FALSE]
  y.tr <- y[tr.row, , drop = FALSE]
  x.te <- X[-tr.row, , drop = FALSE]
  y.te <- y[-tr.row, , drop = FALSE]

  list(
    "x.tr" = x.tr, "y.tr" = y.tr,
    "x.te" = x.te, "y.te" = y.te,
    "call" = call
  )
}

#' Generate Simulation Data for Benchmarking Sparse Regressions (Cox Model)
#'
#' Generate simulation data for benchmarking sparse Cox regression models.
#'
#' @param n Number of observations.
#' @param p Number of variables.
#' @param rho Correlation base for generating correlated variables.
#' @param coef Vector of non-zero coefficients.
#' @param snr Signal-to-noise ratio (SNR).
#' @param p.train Percentage of training set.
#' @param seed Random seed for reproducibility.
#'
#' @return List of \code{x.tr}, \code{x.te}, \code{y.tr}, and \code{y.te}.
#'
#' @author Nan Xiao <\url{https://nanx.me}>
#'
#' @references
#' Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2011).
#' Regularization Paths for Cox's Proportional Hazards Model via
#' Coordinate Descent. \emph{Journal of Statistical Software}, 39(5), 1--13.
#'
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @importFrom survival Surv
#'
#' @export msaenet.sim.cox
#'
#' @examples
#' dat <- msaenet.sim.cox(
#'   n = 300, p = 500, rho = 0.6,
#'   coef = rep(1, 10), snr = 3, p.train = 0.7,
#'   seed = 1001
#' )
#'
#' dim(dat$x.tr)
#' dim(dat$x.te)
#' dim(dat$y.tr)
#' dim(dat$y.te)
msaenet.sim.cox <- function(
  n = 300, p = 500,
  rho = 0.5, coef = rep(0.2, 50), snr = 1,
  p.train = 0.7, seed = 1001) {

  call <- match.call()

  set.seed(seed)

  sigma <- matrix(0, p, p)
  corvec <- function(i, p, rho) rho^(abs(i - 1L:p))
  for (i in 1:p) sigma[i, ] <- corvec(i, p, rho)

  X <- rmvnorm(n, rep(0, p), sigma)

  # non-zero coefficients
  beta0 <- matrix(c(coef, rep(0, (p - length(coef)))))

  snr.numerator <- as.vector(t(beta0) %*% sigma %*% beta0)
  snr.denominator <- snr.numerator / snr
  k <- sqrt(snr.denominator)
  Z <- rnorm(n, 0, 1)
  eps <- matrix(k * Z)

  y.true <- exp((X %*% beta0) + eps)
  time.censor <- exp(k * Z)
  time.record <- pmin(y.true, time.censor)
  event <- as.integer(y.true <= time.censor)

  y <- Surv(time = time.record, event = event, type = "right")
  colnames(y) <- c("time", "status") # for glmnet

  # training / test set splitting
  tr.row <- sample(1L:n, round(n * p.train), replace = FALSE)

  x.tr <- X[tr.row, , drop = FALSE]
  y.tr <- y[tr.row, , drop = FALSE]
  x.te <- X[-tr.row, , drop = FALSE]
  y.te <- y[-tr.row, , drop = FALSE]

  list(
    "x.tr" = x.tr, "y.tr" = y.tr,
    "x.te" = x.te, "y.te" = y.te,
    "call" = call
  )
}
