#' Refit Kernel ODE with a Given Regulatory Network
#'
#' Refits the Kernel ODE model using a pre-specified regulatory network, bypassing the network estimation step.
#' @inheritParams kernelODE_step2
#' @param adj_matrix An adjacency matrix (`p` × `p`) representing the regulatory
#'   network to use for refitting. Entry `(k, j) = 1` indicates variable `k`
#'   regulates variable `j`. Typically obtained from a previous Kernel ODE
#'   estimation. If `NULL`, defaults to a fully connected network.
#'   Although technically possible, this function is *not* intended for
#'   estimating \eqn{F_j} or inferring a regulatory network.
#'
#' @returns A list with components:
#' \describe{
#'   \item{`metrics`}{Recovery metrics for the trajectories under the given
#'     network, including the overall \eqn{R^2} and variable-specific \eqn{R^2}
#'     values. See [assess_recov_traj()].}
#'   \item{`Y_refit`}{Recovered trajectories from the refitted model, in the same
#'     format as `Y`.}
#'   }
#' @export
#'
#' @examples
refit_kernel_ODE <- function(Y,
                             obs_time,
                             yy_smth,
                             tt,
                             kernel,
                             kernel_params,
                             interaction_term,
                             adj_matrix,
                             theta_initial = NULL,
                             nzero_thres = NULL,
                             tol = 0.001,
                             max_iter = 10,
                             verbose = 0){
  if (is.null(adj_matrix)) {
    warning("No network is given. Set to fully connected network.")
    adj_matrix <- matrix(1, nrow = ncol(Y), ncol = ncol(Y))
  }

  # refit KODE using the given network
  res_KODE_refit <- kernelODE_step2(Y = Y,
                                    obs_time = obs_time,
                                    yy_smth = yy_smth,
                                    tt = tt,
                                    kernel = kernel,
                                    kernel_params = kernel_params,
                                    interaction_term = interaction_term,
                                    theta_initial = theta_initial,
                                    adj_matrix = adj_matrix,
                                    nzero_thres = nzero_thres,
                                    tol = tol,
                                    max_iter = max_iter,
                                    verbose = verbose)

  # recover trajectory and compute the metrics
  obs_idx <- .map_to(obs_time, tt)
  Y_refit <- sapply(1:ncol(Y), function(j){
    res_traj_j <- evaluate_Fj(bj = res_KODE_refit$res_bj[j],
                              cj = res_KODE_refit$res_cj[,j],
                              interaction_term = interaction_term,
                              kernel = kernel,
                              kernel_params = kernel_params,
                              obs_time = obs_time,
                              theta_j = res_KODE_refit$res_theta[,j],
                              tt = tt,
                              Yj = Y[,j],
                              yy_smth = yy_smth)
    res_traj_j$yy_est[obs_idx]
  })

  metrics <- assess_recov_traj(Y = Y,
                               Y_est = Y_refit)

  return (list(metrics = metrics,
               Y_refit = Y_refit))
}


