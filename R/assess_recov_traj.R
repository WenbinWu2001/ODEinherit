#' Assess trajectory recovery performance for a single cell.
#'
#' Computes \eqn{R^2} metrics comparing observed trajectories with those
#' recovered by Kernel ODE.
#'
#'
#' @param Y A numeric matrix of dimension (`n`, `p`) containing the observed
#'   trajectories for a cell.
#' @param Y_est A numeric matrix of dimension (`n`, `p`) containing the
#'   trajectories recovered from Kernel ODE.
#'
#' @return A list containing:
#' \describe{
#'   \item{`R2`}{An overall \eqn{R^2} value (a scalar) as the average of the variable-specific \eqn{R^2} values (stored in `R2_per_var_vec`), representing the overall proportion of variance explained by the recovered trajectories.}
#'   \item{`R2_per_var_vec`}{A numeric vector of length `p` giving the variable-specific \eqn{R^2} values.}
#' }
#'
#' @details The \eqn{R^2} metric for variable \eqn{j} is defined as
#' \deqn{R_j^{2} = \max \left\{ 1 -
#' \frac{\mathrm{MSS}_{j,\mathrm{res}}}{\mathrm{MSS}_{j,\mathrm{tot}}}, \; 0
#' \right\},} where \deqn{\mathrm{MSS}_{j,\mathrm{res}} =
#' \frac{1}{n}\sum_{i=1}^n \left( y_{ij} - \hat{y}_{ij} \right)^2,}
#' \deqn{\mathrm{MSS}_{j,\mathrm{tot}} = \frac{1}{n}\sum_{i=1}^n \left( y_{ij} -
#' \bar{y}_j \right)^2.}
#'
#' The overall \eqn{R^2} is the average of the variable-specific \eqn{R^2}
#' values: \deqn{R^{2} = \frac{1}{p}\sum_{j=1}^{p} R_j^{2}.}
#'
#' @examples
#' set.seed(1)
#' obs_time <- seq(0, 1, length.out = 10)
#' Y <- cbind(sin(2 * pi * obs_time), cos(4 * pi * obs_time)) + 0.1 * matrix(rnorm(20), 10, 2)  # each col is a variable
#' Y_est <- sapply(1:ncol(Y), function(j){lm(Y[,j] ~ obs_time)$fitted.values})  # linear regression
#' assess_recov_traj(Y = Y, Y_est = Y_est)
#'
#' @export
assess_recov_traj <- function(Y,
                              Y_est){
  if (any(dim(Y) != dim(Y_est))) {stop("Y and Y_est should share the same dimension.")}

  # reshape Y if it is a numeric vector
  if (is.vector(Y)) {Y <- matrix(Y, ncol = 1)}
  if (is.vector(Y_est)) {Y_est <- matrix(Y_est, ncol = 1)}

  p <- ncol(Y)

  R2_per_var_vec <- sapply(1:p, function(j){
    MS_res_est <- mean((Y[,j] - Y_est[,j])^2)
    MS_total <- mean((Y[,j] - mean(Y[,j]))^2)

    # If R2_est < 0 (i.e., the reconstructed traj. is no better than the sample mean), truncate at 0.
    R2_est <- max(0, 1 - MS_res_est/MS_total)
    R2_est
  })

  R2 <- mean(R2_per_var_vec)  # overall R2 for a cell

  return (list(R2 = R2,
               R2_per_var_vec = R2_per_var_vec))
}
