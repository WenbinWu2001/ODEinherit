#' Smoothing Spline Estimates (Step 1 of Kernel ODE)
#'
#' Computes smoothed trajectories and their derivatives using cubic smoothing splines.
#' This function serves as Step 1 in the Kernel ODE pipeline.
#'
#'
#' @param Y A numeric matrix of dimension (`n`, `p`), where each column corresponds to the observed trajectory of a variable. Rows align with `obs_time`.
#' @param obs_time A numeric vector of length `n` representing observation time points.
#' @param tt A numeric vector representing a finer time grid used for evaluating the smoothed trajectories and their derivatives.
#'
#' @return A list with components:
#' \describe{
#'   \item{`yy_smth`}{A numeric matrix of dimension (`length(tt)`, `p`), where each column contains the smoothed trajectory of a variable evaluated on `tt`.}
#'   \item{`init_vals_smth`}{A numeric vector of length `p` containing the estimated initial values (at time 0) for each variable.}
#'   \item{`deriv_smth`}{A numeric matrix of dimension (`length(tt)`, `p`), where each column contains the smoothed first order derivative of a variable evaluated on `tt`.}
#' }
#'
#' @references
#' Original implementation adapted from <https://github.com/ChenShizhe/GRADE>
#'
#' @export
#'
#' @examples
#' # Example usage:
#' set.seed(1)
#' obs_time <- seq(0, 1, length.out = 10)
#' Y <- cbind(sin(2 * pi * obs_time), cos(4 * pi * obs_time)) + 0.1 * matrix(rnorm(20), 10, 2)  # each col is a variable
#' tt <- seq(0, 1, length.out = 100)
#' result <- smoother_SS(obs_time, Y, tt)
#' matplot(tt, result$yy_smth, type = "l", lty = 1, col = 1:2)
kernelODE_step1 <- function(Y, obs_time, tt){
  p <- ncol(Y)
  times_e <- c(0,tt)  # also fit an initial value
  init_vals_smth <- rep(NA, p)
  yy_smth <- matrix(NA, length(tt), p)  # on tt
  deriv_smth <- matrix(NA, length(tt), p)  # on tt
  for(j in 1:p){
    SS_model <- stats::smooth.spline(obs_time, Y[,j], all.knots = T)  # a cubic smoothing spline where all points in `obs_time` are used as knots
    gcvscore <- SS_model$cv.crit
    pred_times_e <- stats::predict(SS_model, times_e)$y

    init_vals_smth[j] <- pred_times_e[1]
    yy_smth[,j] <- pred_times_e[-1]
    deriv_smth[,j] <- (stats::predict(SS_model,times_e,deriv=1)$y)[-1]
  }

  return (list(yy_smth = yy_smth,
               init_vals_smth = init_vals_smth,
               deriv_smth = deriv_smth))
}


