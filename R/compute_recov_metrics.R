compute_recov_metrics <- function(Y_list,
                                  Y_est_list,
                                  Y_smth_list,
                                  Y_raw_list){
  p <- ncol(Y_list[[1]])
  K <- length(Y_list)
  metrics_list <- list()
  for (k in 1:K){
    metrics_list[[k]] <- list()

    for (j in 1:p){
      # w.r.t. denoised series
      MS_res_est <- mean((Y_list[[k]][,j] - Y_est_list[[k]][,j])^2)  # mean squared residuals of est. trajectory w.r.t. observations (variance in our estimated trajectory / "recovered signals")
      MS_res_smth <- mean((Y_list[[k]][,j] - Y_smth_list[[k]][,j])^2)  # mean squared residuals of smoothing estimates w.r.t. observations (variance of "true signals")
      MS_total <- mean(Y_list[[k]][,j]^2)  # mean squared total of the observations. Note that `Y_list` has centered columns.

      # w.r.t. raw series
      MS_res_est_raw <- mean((Y_raw_list[[k]][,j] - Y_est_list[[k]][,j])^2)  # mean squared residuals of est. trajectory w.r.t. observations (variance in our estimated trajectory / "recovered signals")
      MS_res_smth_raw <- mean((Y_raw_list[[k]][,j] - Y_smth_list[[k]][,j])^2)  # mean squared residuals of smoothing estimates w.r.t. observations (variance of "true signals")
      MS_total_raw <- mean(Y_raw_list[[k]][,j]^2)  # mean squared total of the observations. Note that `Y_list` has centered columns.

      # R2 for all variables of cell k
      R2_est <- max(0, 1 - MS_res_est/MS_total)  # if the reconstructed trajectory is nonsense (so MS_res_est > MS_total), truncate at 0
      R2_smth <- max(0, 1 - MS_res_smth/MS_total)
      R2_est_raw <- max(0, 1 - MS_res_est_raw/MS_total_raw)
      R2_smth_raw <- max(0, 1 - MS_res_smth_raw/MS_total_raw)

      metrics_list[[k]][[j]] <- list(MS_res_est = MS_res_est,
                                     MS_res_smth = MS_res_smth,
                                     MS_total = MS_total,
                                     R2_est = R2_est,
                                     R2_smth = R2_smth,
                                     MS_res_est_raw = MS_res_est_raw,
                                     MS_res_smth_raw = MS_res_smth_raw,
                                     MS_total_raw = MS_total_raw,
                                     R2_est_raw = R2_est_raw,
                                     R2_smth_raw = R2_smth_raw)
    }
  }

  return (metrics_list)
}



create_contingency_table <- function(network_est, network_true) {
  # Calculate true positives, false positives, true negatives, false negatives
  TP <- sum(network_est == 1 & network_true == 1)
  TN <- sum(network_est == 0 & network_true == 0)
  FP <- sum(network_est == 1 & network_true == 0)
  FN <- sum(network_est == 0 & network_true == 1)

  # Create contingency table
  contingency_table <- matrix(c(TP, FP, FN, TN), nrow = 2, byrow = TRUE,
                              dimnames = list(c("Prediction: Pos", "Prediction: Neg"),
                                              c("Actual: Pos", "Actual: Neg")))

  # Calculate rates
  FPR <- FP / (FP + TN)  # False Positive Rate
  FNR <- FN / (FN + TP)  # False Negative Rate

  return(list(contingency_table = contingency_table,
              FPR = FPR,
              FNR = FNR))
}




### deprecated version
assess_recov_traj <- function(Y,
                              Y_est,
                              Y_smth,
                              obs_time = NULL){
  # Assess the trajectory recovery performance for one cell.
  # Specifically, compute the proportion of variation in the observed and smoothed trajectories that is explained by the reconstructed trajectories, denoted as `R2_est` and `R2_smth`, respectively.
  # `Y` is the observed trajectories or raw trajectories on the observed time grid for this cell. It has dimension (n, p).
  # `Y_est` is the reconstructed trajectories on the observed time grid from kernel ODE. It has dimension (n, p).
  # `Y_smth` is the smoothed trajectories on the observed time grid from step 1, i.e. `yy_smth` evaluated on `obs_time`. It has dimension (n, p).
  # `obs_time` is the n observed time points (aligned with the trajectories above), of length n. By default it is set to be the evenly spaced grid on [0,1] (`1/n * (1:n)`).

  if (any(dim(Y) != dim(Y_est)) | any(dim(Y) != dim(Y_smth))) {stop("Y, Y_est, and Y_smth should have the same dimensions.")}

  p <- dim(Y)[2]
  R2_per_var <- list()
  for (j in 1:p){
    MS_res_est <- mean((Y[,j] - Y_est[,j])^2)  # mean squared residuals of est. trajectory w.r.t. observations (variance in our estimated trajectory / "recovered signals")
    MS_res_smth <- mean((Y[,j] - Y_smth[,j])^2)  # mean squared residuals of smoothing estimates w.r.t. observations (variance of "true signals")
    MS_total <- mean((Y[,j] - mean(Y[,j]))^2)  # mean squared total of the observations.

    # R2 for all variables for this cell
    R2_est <- max(0, 1 - MS_res_est/MS_total)  # if the reconstructed trajectory is nonsense (so MS_res_est > MS_total), truncate at 0
    R2_smth <- max(0, 1 - MS_res_smth/MS_total)

    R2_per_var[[j]] <- list(MS_res_est = MS_res_est,
                            MS_res_smth = MS_res_smth,
                            MS_total = MS_total,
                            R2_est = R2_est,
                            R2_smth = R2_smth)
  }

  R2_overall <- list(
    R2_est = mean(sapply(R2_per_var, function(obj){obj$R2_est})),
    R2_smth = mean(sapply(R2_per_var, function(obj){obj$R2_smth})),
    R2_est_per_var = sapply(R2_per_var, function(obj){obj$R2_est}),
    R2_smth_per_var = sapply(R2_per_var, function(obj){obj$R2_smth})
  )

  # return (list(R2_per_var = R2_per_var,
  #              R2_overall = R2_overall))

  return (R2_per_var)
}
