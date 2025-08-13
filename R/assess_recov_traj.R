#' Title
#'
#' @param Y
#' @param Y_est
#'
#' @returns
#' @export
#'
#' @examples
assess_recov_traj <- function(Y,
                              Y_est){
  # Assess the trajectory recovery performance for one cell. Obtain the R2 metrics.
  # Specifically, compute the proportion of variation in the observed and smoothed trajectories that is explained by the reconstructed trajectories, denoted as `R2_est` and `R2_smth`, respectively.
  # `Y` is the observed trajectories or raw trajectories on the observed time grid for this cell. It has dimension (n, p).
  # `Y_est` is the reconstructed trajectories on the observed time grid from kernel ODE. It has dimension (n, p).

  if (any(dim(Y) != dim(Y_est))) {stop("Y and Y_est should share the same dimension.")}

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
