kernelODE <- function(Y,
                      kernel,
                      kernel_params,
                      interaction_term = FALSE,
                      theta_initial = NULL,
                      adj_matrix = NULL,
                      nzero_thres = NULL,
                      eval_edge_R2 = FALSE,
                      type_smooth = "smoothspline",
                      type_data = "Perturbations",
                      tol = 0.001,
                      max_iter = 10,
                      verbose = 0,
                      plot_theta_traj = FALSE){
  ## Kernel ODE algorithm
  ## Dai X, Li L. Kernel Ordinary Differential Equations. J Am Stat Assoc. 2022;117(540):1711-1725. doi: 10.1080/01621459.2021.1882466. Epub 2021 Apr 27. PMID: 36845295; PMCID: PMC9949731.
  ## `Y` (nxp) is the a matrix or data.frame of p variables observed at n time points. These data are assumed to lie in the range [0,1].
  ## `kernel` specifies the main kernel function that is used for all variables.
  ## `kernel_params` is a list of p lists (or a list of one list, and this kernel parameter is used for all variables).
  ## `kernel_params[[k]]` contains the parameters of variable k for the specified kernel. Check `kernel_main` to see the parameters needed for each kernel.
  ## `nzero_thres` either `NULL` or a scalar between 0 and 1 (both inclusive).
  ## If `nzero_thres` is `NULL`, regular non-negative Lasso is implemented as in the paper.
  ## Otherwise, it controls the proportion of nonzero coefficients in theta_j to be <= `nzero_thres` during theta_j estimation using Lasso.
  ## `interaction_term`: logical, whether to include two-way interaction effects.
  ## `theta_initial`: the initialization of theta. It should be a p^2 x p matrix if `interaction_term` = TRUE (default) and pxp if `interaction_term` = FALSE.
  ##  `type_smooth` and `type_data` are configurations for smoothing the observed trajectories in `Y` in step 1. Details are given by GRADE.
  ## `tol` and `max_iter` specify the stopping criteria
  ## `verbose`: either 0, 1, or 2. Controls the output messages indicating the process of the algorithm. A larger value gives more detailed outputs.
  ## `plot_theta_traj`: logical, whether to plot the theta_j trajectories (and selected lambda's) from Lasso, for ALL (variables j's) and during ALL iterations. Note that this significantly slows the speed of the algorithm.


  ### Start: commented for debug purpose ###
  # # check inputs
  # if (!(is.matrix(Y) | is.data.frame(Y))) {stop("Y should be a nxp matrix or data.frame.")}
  # if (nrow(Y) < 2) {stop("Y should contain at least 2 observations.")}
  # if (!is.null(adj_matrix)) {
  #   if (any(dim(adj_matrix) != c(ncol(Y), ncol(Y)))) {stop("adj_matrix should be pxp.")}
  #   if (!all(adj_matrix %in% c(0,1))) {stop("adj_matrix should be binary.")}
  # }
  # if (!is.null(nzero_thres)){
  #   if (!(nzero_thres>=0 & nzero_thres<=1)){
  #     stop("nzero_thres should either be NULL or a scalar between 0 and 1 (both inclusive).")
  #   }
  # }
  #
  # # dimensionality check of theta_initial
  # if (!is.null(theta_initial)){  # if user specifies theta_initial
  #   if (interaction_term){
  #     if (!(is.matrix(theta_initial) & all(dim(theta_initial) == c(ncol(Y)^2, ncol(Y))))) {stop("Invalid theta_initial")}
  #   } else {
  #     if (!(is.matrix(theta_initial) & all(dim(theta_initial) == c(ncol(Y), ncol(Y))))) {stop("Invalid theta_initial")}
  #   }
  # }
  # # structure check of kernel_params
  # if (!(is.list(kernel_params) & all(sapply(kernel_params, is.list)))) {stop("kernel_params should be of a list of lists.")}  # must be a list of lists
  # if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y), kernel_params[[1]], simplify = FALSE)}  # if only one parameter set is specified, it is recycled for all variables
  # if (length(kernel_params) != ncol(Y)) {stop("kernel_params should be of length p.")}
  #
  # # center the data
  # Y <- scale(Y, center = TRUE, scale = FALSE)
  #
  ### End: commented for debug purpose ###



  # time grids
  n <- nrow(Y)
  obs_time <- 1/n * (1:n)
  delta <- 0.001
  tt <- seq(delta, 1, by = delta)  # time grid for numerical integration, does not include 0

  # # Step 1: Obtaining the smoothed trajectories
  # res_smth <- kernelODE_step1(Y = Y,
  #                             obs_time = obs_time,
  #                             tt = tt)
  # yy_smth <- res_smth$yy_smth
  # # DO NOT CENTER yy_smth!

  smthed <- smoother_SS(obs_time = obs_time,
                        Y = Y,
                        tt = tt)
  yy_smth <- smthed$yy_smth

  # Step 2: Iterative algorithm for fitting the ODE system and inferring the regulatory network
  res_step2 <- kernelODE_step2(Y = Y,
                               obs_time = obs_time,
                               yy_smth = yy_smth,
                               tt = tt,
                               kernel = kernel,
                               kernel_params = kernel_params,
                               interaction_term = interaction_term,
                               theta_initial = theta_initial,
                               nzero_thres = nzero_thres,
                               tol = 0.001,
                               max_iter = 10,
                               verbose = 1)

  # res_step2$res_smth <- res_smth  # add the smoothing results from step 1



  return (res_step2)
}
