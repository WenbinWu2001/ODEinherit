refit_kernel_ODE <- function(Y,
                             obs_time,
                             yy_smth,
                             tt,
                             kernel,
                             kernel_params,
                             interaction_term,
                             theta_initial =NULL,
                             adj_matrix = NULL,
                             nzero_thres = NULL,
                             tol = 0.001,
                             max_iter = 10,
                             verbose = 0){
  # refit KODE using the given network
  res_KODE_refit <- LiveImagingCells::kernelODE_step2(Y = Y,
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
  obs_idx <- .map_to(res_KODE_refit$config$obs_time, tt)
  Y_est_refit <- sapply(1:ncol(Y), function(j){
    res_traj_j <- LiveImagingCells::evaluate_Fj(bj = res_KODE_refit$res_bj[j],
                                                cj = res_KODE_refit$res_cj[,j],
                                                interaction_term = interaction_term,
                                                kernel = kernel,
                                                kernel_params = kernel_params,
                                                obs_time = res_KODE_refit$config$obs_time,
                                                theta_j = res_KODE_refit$res_theta[,j],
                                                tt = tt,
                                                Yj = Y[,j],
                                                yy_smth = yy_smth)
    res_traj_j$yy_est[obs_idx]
  })

  metrics <- assess_recov_traj(Y = Y,
                               Y_est = Y_est_refit,
                               Y_smth = yy_smth[obs_idx,])

  return (list(metrics = metrics,
               Y_est_refit = Y_est_refit))
}


