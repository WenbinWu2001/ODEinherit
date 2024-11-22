kernelODE_step1 <- function(Y,
                            obs_time,
                            tt,
                            type_smooth = "smoothspline",
                            type_data = "Perturbations"){
  ## Kernel ODE Step 1: Smoothing spline estimation (single sample version).
  ## `Y` is a matrix or data.frame of dimension (n, p) with n samples and p variables. Columns of Y are centered to 0 (IMPORTANT).
  ## `obs_time` is a numeric vector specifying the time points when samples in `Y` are observed. Should lie in the range [0,1].
  ## `tt` is a numeric vector specifying the integration grid on [0, 1] on which the smoothing splines will be evaluated. The results are fed into step 2 to approximate the integrals in the algorithm.
  ## Returns `yy_smth`: interpolated values of the estimated smoothing splines on the integration grid `tt`.
  # if (!(is.matrix(Y) | is.data.frame(Y))) {stop("Y should be a nxp matrix or dataframe.")}
  # if (nrow(Y) < 2) {stop("Y should contain at least 2 observations.")}
  # if ((obs_time[1] < 0) | obs_time[length(obs_time)] > 1) {stop("obs_time should lie in the range [0,1].")}
  # if ((tt[1] < 0) | tt[length(tt)] > 1) {stop("tt should lie in the range [0,1].")}

  # Y <- scale(Y, center = TRUE, scale = FALSE)  # ensure Y is centered at 0
  # SS_fit_list <- .SS_estimation(obs_time, Y)  # smoothing spline estimation
  # yy_smth <- .SS_interpolation(SS_fit_list, tt)  # interpolation of variables values

  # observations <- list(cbind(obs_time, Y))  # some specific form of input for smoothX
  # smthed<-smoothX(observations=observations,
  #                 times_e = c(0,tt),  # interpolate on time 0 (for initial value) and tt
  #                 type_smooth = type_smooth,
  #                 type_data = type_data)
  #
  # init_vals_smth <- smthed$Xhat[[1]][1,]  # p initial values of the smoothing trajectories (one for each variable).
  #
  # yy_smth <- smthed$Xhat[[1]][-1,]  # smoothing trajectories on `tt`
  # # DO NOT CENTER yy_smth!

  # observations <- list(cbind(obs_time, Y))  # some specific form of input for smoothX
  # smthed<-smoothX(observations=observations,
  #                 times_e = tt,  # interpolate on time 0 (for initial value) and tt
  #                 type_smooth = type_smooth,
  #                 type_data = type_data)
  # yy_smth <- smthed$Xhat[[1]]
  # # DO NOT CENTER yy_smth!
  # init_vals_smth <- NULL

  smthed <- smoother_SS(obs_time = obs_time,
                        Y = Y,
                        tt = tt)
  yy_smth <- smthed$yy_smth
  init_vals_smth <- smthed$init_vals_smth
  Fj_smth <- smthed$deriv_smth  # est. deriv of the smoothing trajectories on `tt`
  ## if GRADE somehow does not return it, compute manually
  # array_lag1_temp <- array(NA, dim = c(len, p))
  # array_lag1_temp[1,] <- init_vals_smth
  # array_lag1_temp[2:len,] <- yy_smth[1:(len-1),]
  # Fj_smth <- (yy_smth - array_lag1_temp) / delta  # (len, p)

  smoother_config <- list(type_smooth = type_smooth,
                          type_data = type_data)

  res_smth <- list(init_vals_smth = init_vals_smth,
                   yy_smth = yy_smth,
                   tt = tt,
                   Fj_smth = Fj_smth,
                   smoother_config = smoother_config)
  return(res_smth)
}

