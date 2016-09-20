#' Generate Simulation Data (Gaussian Case)
#'
#' Generate simulation data (Gaussian) following the
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
#' @author Nan Xiao <\url{http://nanx.me}>
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
#' dat = msaenet.sim.gaussian(n = 300, p = 500, rho = 0.6,
#'                            coef = rep(1, 10), snr = 3, p.train = 0.7,
#'                            seed = 1001)
#' dim(dat$x.tr)
#' dim(dat$x.te)

msaenet.sim.gaussian = function(n = 300, p = 500,
                                rho = 0.5, coef = rep(0.2, 50), snr = 1,
                                p.train = 0.7, seed = 1001) {

  call = match.call()

  set.seed(seed)

  sigma = matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      sigma[i, j] = rho^(abs(i - j))
    }
  }

  X = rmvnorm(n = n, mean = rep(0, p), sigma = sigma)

  # non-zero coefficients
  beta0 = matrix(c(coef, rep(0, (p - length(coef)))))

  snr.numerator = as.vector(t(beta0) %*% sigma %*% beta0)
  snr.denominator = snr.numerator/snr
  sd = sqrt(snr.denominator)
  eps = matrix(rnorm(n, 0, sd))

  y = (X %*% beta0) + eps

  # training / test set splitting
  tr.row = sample(1L:n, round(n * p.train), replace = FALSE)

  x.tr = X[tr.row, ]
  y.tr = as.matrix(y[tr.row, ])
  x.te = X[-tr.row, ]
  y.te = as.matrix(y[-tr.row, ])

  return(list('x.tr' = x.tr, 'y.tr' = y.tr,
              'x.te' = x.te, 'y.te' = y.te,
              'call' = call))

}
