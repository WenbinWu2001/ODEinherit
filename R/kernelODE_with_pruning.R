kernelODE_step2_with_pruning <- function(Y,
                                         obs_time,
                                         yy_smth,
                                         tt,
                                         kernel,
                                         kernel_params,
                                         interaction_term = FALSE,
                                         theta_initial = NULL,
                                         adj_matrix = NULL,
                                         nzero_thres = NULL,
                                         prune = TRUE,
                                         prune_thres = 0.05,
                                         eval_edge_R2 = FALSE,
                                         tol = 0.001,
                                         max_iter = 10,
                                         verbose = 0,
                                         plot_theta_traj = FALSE){

  ## kernel ODE Step 2: Iterative optimization algorithm (single sample version)
  # structure check of kernel_params
  if (!(is.list(kernel_params) & all(sapply(kernel_params, is.list)))) {stop("kernel_params should be of a list of lists.")}  # must be a list of lists
  if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y), kernel_params[[1]], simplify = FALSE)}  # if only one parameter set is specified, it is used for all variables
  if (length(kernel_params) != ncol(Y)) {stop("kernel_params should be of length p.")}

  n <- nrow(Y)
  p <- ncol(Y)

  obs_idx <- .map_to(obs_time, tt)

  # just in case `Y` and `yy_smth` are not centered at 0
  Y <- scale(Y, center = TRUE, scale = FALSE)
  # DO NOT CENTER yy_smth!

  # plot config
  if (plot_theta_traj) {
    graphics::par(mfrow = c(max_iter + 1, ifelse(interaction_term, p^2, p) + 1))

    # print the label for each column (element of theta) corresponding to either a main effect of an interaction effect
    # and its index in theta_j (so it corresponds to theta_j[index])
    if (!interaction_term){
      col_name_values <- paste0("main_", 1:p)
    } else {
      col_name_values <- outer(1:p, 1:p, paste, sep = "_")
      col_name_values <- apply(col_name_values, c(1, 2), function(name) paste0("inter_", name))
      diag(col_name_values) <- paste0("main_", 1:p)
      col_name_values <- c(t(col_name_values))
    }

    graphics::plot(1, type='n', xlim=c(0, 1), ylim=c(0, 1), xlab='', ylab='', axes=FALSE)
    for (col_name in col_name_values){
      graphics::plot(1, type='n', xlim=c(0, 1), ylim=c(0, 1), xlab='', ylab='', axes=FALSE)
      graphics::text(x=0.5, y=0.5, paste0(col_name, " (", which(col_name_values == col_name), ")"), cex=3)
    }
  }


  # recyclable quantities
  B <- .construct_B(obs_time)
  qr_decomp <- qr(B)
  Q <- qr.Q(qr_decomp, complete = TRUE)
  Q1 <- matrix(Q[,1], nrow = n, ncol = 1)
  Q2 <- Q[,2:n]
  R <- qr.R(qr_decomp)[1,1]

  if (verbose > 0) {cat("-------- constructing Sigma_k_kl --------\n")}
  Sigma_k_kl <- .construct_Sigma_k_kl(interaction_term = interaction_term,
                                      kernel = kernel,
                                      kernel_params = kernel_params,
                                      obs_time = obs_time,
                                      tt = tt,
                                      yy_smth = yy_smth)

  # Iterative optimization
  res_best_eta <- rep(NA, p)
  res_bj <- rep(NA, p)
  res_cj <- matrix(NA, nrow = n, ncol = p)
  res_best_kappa <- rep(NA, p)  # in theta_j estimation step

  # initialization: if no specified theta_initial, then initialize to 1 for all entries.
  if (is.null(theta_initial)){
    num_row <- ifelse(interaction_term, yes = p^2, no = p)
    theta_initial <- matrix(1, nrow = num_row, ncol = p)  # p^2 x p or p x p
  }
  res_theta <- theta_initial

  # iteration begins
  for (num_iter in 1:max_iter){
    if (verbose > 0) {cat("\n-------- iteration", num_iter, " --------\n")}

    res_theta_prev <- res_theta

    # given theta_j, estimate F_j
    if (verbose > 0) {cat("-------- estimating Fj's --------\n")}

    for (j in 1:p) {
      if (verbose > 0) {cat("Working on variable", j, "\n")}

      Yj <- Y[,j]
      theta_j <- res_theta[,j]
      Sigma <- .construct_Sigma(Sigma_k_kl = Sigma_k_kl,
                                theta_j = theta_j)  # update Sigma using new theta_j
      res_fun_est <- function_estimation(Yj = Yj,
                                         Sigma = Sigma,
                                         Q1 = Q1,
                                         Q2 = Q2,
                                         R = R,
                                         verbose = verbose)

      res_best_eta[j] <- res_fun_est$best_eta
      res_bj[j] <- res_fun_est$bj
      res_cj[,j] <- res_fun_est$cj
    }

    # given F_j, estimate theta_j
    if (verbose > 0) {cat("-------- estimating theta_j's --------\n")}
    if (plot_theta_traj){
      graphics::plot(1, type='n', xlim=c(0, 1), ylim=c(0, 1), xlab='', ylab='', axes=FALSE)  # Create an empty plot
      graphics::text(x=0.5, y=0.5, paste0("Iter ", num_iter), cex=3)  # Add text at the center of the plot
    }
    for (j in 1:p){
      Yj <- Y[,j]
      bj <- res_bj[j]
      cj <- res_cj[,j]
      eta_j <- res_best_eta[j]
      G <- .construct_G(cj = cj,
                        interaction_term = interaction_term,
                        Sigma_k_kl = Sigma_k_kl)

      # non-negative lasso fit
      adj_col <- adj_matrix[,j]  # variable j corresponds to column j of the adjacency matrix. Note that `adj_col` is NULL if `adj_matrix` is NULL.
      res_theta_j_est <- theta_j_estimation(bj = bj,
                                            B = B,
                                            cj = cj,
                                            eta_j = eta_j,
                                            G = G,
                                            interaction_term = interaction_term,
                                            Yj = Yj,
                                            adj_col = adj_col,
                                            nzero_thres = nzero_thres,
                                            plot_theta_traj = plot_theta_traj)
      res_theta[,j] <- res_theta_j_est$theta_j
      res_best_kappa[j] <- res_theta_j_est$best_kappa
    }

    # if (any(dim(theta_j) != c(ifelse(interaction_term, p^2, p), p))) {stop("Incorrect dimension of `theta_j`.")}

    # check improvement in this iteration
    improvement <- norm(res_theta - res_theta_prev, type = "F") / norm(res_theta_prev, type = "F")

    if (verbose > 0) {
      cat("theta improvement:", round(improvement, 4), "\n")
      cat("F-norm of current theta:", round(norm(res_theta, type = "F"), 4), "\n")
    }

    if (is.na(improvement)){
      # usually because `norm(res_theta_prev, type = "F")` is 0 (i.e. `res_theta_prev` is all 0) at the first few iterations
      # OR because the jth column of the input network `adj_matrix` is all zero (i.e. no variable affects variable j) so that theta_j is all zero.
      # warning("improvement is NA, iteration terminated")
      # improvement <- 0  # break the loop by the if statement below
      warning("improvement is NA")
      improvement <- 2*tol  # continue
    }

    # check stopping criteria
    if (improvement < tol) {break}
  }  # end of for loop

  # after the iterative process, truncate small components of theta_j to 0
  res_theta[res_theta < 0.01] <- 0

  if (verbose > 0) {cat("\n -------- finished -------- \n")}

  # recover original plot config
  if (plot_theta_traj) {graphics::par(mfrow = c(1,1))}


  # obtain the estimated network
  res_network <- theta_to_adj_matrix(interaction_term = interaction_term,
                                     res_theta = res_theta)

  # pruning of the network based on the estimated Fjk components of Fj (equivalent to truncate theta_jk at some level)
  # this pruning procedure corresponds to pruning_new
  if (prune){
    if (verbose > 0) {cat("\n------ Pruning network ------\n")}
    # recover fitted trajectory and compute the R2
    Y_est <- sapply(1:p,
                    FUN = function(j){
                      res_traj_j <- evaluate_Fj(bj = res_bj[j],
                                                cj = res_cj[,j],
                                                interaction_term = interaction_term,
                                                kernel = kernel,
                                                kernel_params = kernel_params,
                                                obs_time = obs_time,
                                                theta_j = res_theta[,j],
                                                tt = tt,
                                                Yj = Y[,j],
                                                yy_smth = yy_smth)
                      res_traj_j$yy_est[obs_idx]
                      })
    metrics_unpruned <- assess_recov_traj(Y = Y,
                                          Y_est = Y_est,
                                          Y_smth = yy_smth[obs_idx,])

    # evaluate the R2 explained by each edge in the estimated network
    mat_temp <- res_network
    mat_temp[res_network == 0] <- NA  # NA indicates an edge is not selected in the input network
    R2_mat <- mat_temp  # (p, p), R2 of each selected edge (in the input network).

    # identify selected edges
    edge_idx <- which(res_network == 1, arr.ind = T)
    colnames(edge_idx) <- c("from", "to")

    # evaluate R2 of each edge and prune
    if (verbose > 0) {cat("\n--- Assessing R2 of each selected edge in the estimated network ---\n")}
    res_network_pruned <- res_network
    if (nrow(edge_idx) != 0){
      # at least one edge is selected in the network
      for (idx in 1:nrow(edge_idx)){
        var_from <- edge_idx[idx, 1]  # affecting variable
        var_to <- edge_idx[idx, 2]  # affected variable

        # for each cell, set the corresponding theta_j entry to be zero and evaluate the recovered trajectories using this theta_j
        j <- var_to  # we only need the quantities related to the affected variable (we study its Fj without the component from the affecting variable)
        theta_j_pruned <- res_theta[,j]
        if (!interaction_term){
          # without interaction term
          theta_j_pruned[var_from] <- 0  # prune this edge
        } else {
          # with interaction term, theta_j is of length p^2 and we set those theta_j entries related to var_from to 0
          temp_mat <- matrix(theta_j_pruned, nrow = p, ncol = p)
          temp_mat[var_from,] <- 0
          temp_mat[,var_from] <- 0
          theta_j_pruned <- c(temp_mat)
        }

        res_traj_j <- evaluate_Fj(bj = res_bj[j],
                                  cj = res_cj[,j],
                                  interaction_term = interaction_term,
                                  kernel = kernel,
                                  kernel_params = kernel_params,
                                  obs_time = obs_time,
                                  theta_j = theta_j_pruned,  # pruned theta_j
                                  tt = tt,
                                  Yj = Y[,j],
                                  yy_smth = yy_smth)

        Yj_est_new <- res_traj_j$yy_est[obs_idx]
        Yj <- Y[,j]

        MS_res_est <- mean((Yj - Yj_est_new)^2)  # mean squared residuals of est. trajectory w.r.t. observations (variance in our estimated trajectory / "recovered signals")
        MS_total <- mean(Yj^2)  # mean squared total of the observations. Note that `Y_list` has centered columns.
        R2_est <- max(0, 1 - MS_res_est/MS_total)  # R2 of the new reconstructed variable j trajectory with var_idx pruned
        R2_edge <- metrics_unpruned[[j]][["R2_est"]] - R2_est  # % variation explained of variable j traj by this edge for this cell w.r.t. the unpruned network

        R2_mat[var_from, var_to] <- R2_edge
      }

      # prune all edges with R2 smaller than `prune_thres`
      prune_edge_idx <- which(R2_mat < prune_thres, arr.ind = T)
      if (nrow(prune_edge_idx) >= 1) {res_network_pruned[prune_edge_idx] <- 0}
    }

    ### print results ###
    if (verbose > 0){
      cat("\n--- R2 of each selected edge in the original network ---\n")
      print(round(R2_mat, 2))

      cat("\n------ Pruning finished ------\n")
      cat("The following edges are pruned:",
          ifelse(nrow(prune_edge_idx) == 0,
                 yes = "None",
                 no = paste(prune_edge_idx[,1], prune_edge_idx[,2], sep = "-", collapse = ", ")))

      cat("\n------ Pruned network ------\n")
      print(res_network_pruned)
    }
  }  # pruning finished

  res_step2 <- list(res_network = res_network,
                    res_network_pruned = res_network_pruned,
                    res_theta = res_theta,
                    res_bj = res_bj,
                    res_cj = res_cj,
                    res_best_eta = res_best_eta,
                    num_iter = num_iter,
                    last_improvement = improvement,
                    config = list(adj_matrix = adj_matrix,
                                  interaction_term = interaction_term,
                                  kernel = kernel,
                                  kernel_params = kernel_params,
                                  max_iter = max_iter,
                                  nzero_thres = nzero_thres,
                                  prune = prune,
                                  prune_thres = prune_thres,
                                  obs_time = obs_time,
                                  theta_initial = theta_initial,
                                  tol = tol,
                                  tt = tt,
                                  Y = Y,
                                  yy_smth = yy_smth))

  return (res_step2)
}


kernelODE_with_pruning <- function(Y,
                                   kernel,
                                   kernel_params,
                                   interaction_term = FALSE,
                                   theta_initial = NULL,
                                   adj_matrix = NULL,
                                   nzero_thres = NULL,
                                   prune = TRUE,
                                   prune_thres = 0.05,
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

  # check inputs
  if (!(is.matrix(Y) | is.data.frame(Y))) {stop("Y should be a nxp matrix or data.frame.")}
  if (nrow(Y) < 2) {stop("Y should contain at least 2 observations.")}
  if (!is.null(adj_matrix)) {
    if (any(dim(adj_matrix) != c(ncol(Y), ncol(Y)))) {stop("adj_matrix should be pxp.")}
    if (!all(adj_matrix %in% c(0,1))) {stop("adj_matrix should be binary.")}
  }
  if (!is.null(nzero_thres)){
    if (!(nzero_thres>=0 & nzero_thres<=1)){
      stop("nzero_thres should either be NULL or a scalar between 0 and 1 (both inclusive).")
    }
  }

  # dimensionality check of theta_initial
  if (!is.null(theta_initial)){  # if user specifies theta_initial
    if (interaction_term){
      if (!(is.matrix(theta_initial) & all(dim(theta_initial) == c(ncol(Y)^2, ncol(Y))))) {stop("Invalid theta_initial")}
    } else {
      if (!(is.matrix(theta_initial) & all(dim(theta_initial) == c(ncol(Y), ncol(Y))))) {stop("Invalid theta_initial")}
    }
  }
  # structure check of kernel_params
  if (!(is.list(kernel_params) & all(sapply(kernel_params, is.list)))) {stop("kernel_params should be of a list of lists.")}  # must be a list of lists
  if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y), kernel_params[[1]], simplify = FALSE)}  # if only one parameter set is specified, it is recycled for all variables
  if (length(kernel_params) != ncol(Y)) {stop("kernel_params should be of length p.")}

  # center the data
  Y <- scale(Y, center = TRUE, scale = FALSE)

  # time grids
  n <- nrow(Y)
  obs_time <- 1/n * (1:n)
  delta <- 0.001
  tt <- seq(delta, 1, by = delta)  # time grid for numerical integration, does not include 0

  # Step 1: Obtaining the smoothed trajectories
  res_smth <- kernelODE_step1(Y = Y,
                              obs_time = obs_time,
                              tt = tt)
  yy_smth <- res_smth$yy_smth

  # Step 2: Iterative algorithm for fitting the ODE system and inferring the regulatory network
  res_step2 <- kernelODE_step2_with_pruning(Y = Y,
                                            obs_time = obs_time,
                                            yy_smth = yy_smth,
                                            tt = tt,
                                            kernel = kernel,
                                            kernel_params = kernel_params,
                                            interaction_term = interaction_term,
                                            theta_initial = theta_initial,
                                            adj_matrix = adj_matrix,
                                            nzero_thres = nzero_thres,
                                            prune = prune,
                                            prune_thres = prune_thres,
                                            eval_edge_R2 = eval_edge_R2,
                                            tol = tol,
                                            max_iter = max_iter,
                                            verbose = verbose,
                                            plot_theta_traj = plot_theta_traj)

  res_step2$res_smth <- res_smth  # add the smoothing results from step 1

  return (res_step2)
}
