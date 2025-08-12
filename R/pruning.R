.kernelODE_step2_single_variable <- function(adj_col,
                                             B,
                                             interaction_term,
                                             Sigma_k_kl,
                                             Yj,
                                             theta_j_initial = NULL,
                                             max_iter = 5,
                                             tol = 0.001){
  # A dedicated function for running the single-sample kernel ODE on a single variable j during the network pruning.
  # `adj_col` is the jth column of `network_est` with a specific edge deleted.
  # `theta_j_initial` is the initial value of theta_j in the iterative optimization. We suggest passing the one obtained from kernel ODE fit. By default. it is set to be all ones.
  # Note that if `adj_col` only contains zero, no variables do not affect variable j and Fj is fitted as a constant.
  n <- length(Yj)

  qr_decomp <- qr(B)
  Q <- qr.Q(qr_decomp, complete = TRUE)
  Q1 <- matrix(Q[,1], nrow = n, ncol = 1)
  Q2 <- Q[,2:n]
  R <- qr.R(qr_decomp)[1,1]

  # initialization
  if (is.null(theta_j_initial)){theta_j_initial <- rep(1, dim(Sigma_k_kl)[3])}  # of length p^2 or p (same as the 3rd dim of Sigma_k_kl)
  theta_j <- theta_j_initial

  # iteration begins
  for (num_iter in 1:max_iter){
    theta_j_prev <- theta_j

    # given theta_j, estimate F_j
    Sigma <- .construct_Sigma(Sigma_k_kl = Sigma_k_kl,
                              theta_j = theta_j)  # update Sigma using new theta_j
    res_fun_est <- function_estimation(Yj = Yj,
                                       Sigma = Sigma,
                                       Q1 = Q1,
                                       Q2 = Q2,
                                       R = R,
                                       verbose = 0)

    # given F_j, estimate theta_j
    bj <- res_fun_est$bj
    cj <- res_fun_est$cj
    eta_j <- res_fun_est$best_eta
    G <- .construct_G(cj = cj,
                      interaction_term = interaction_term,
                      Sigma_k_kl = Sigma_k_kl)

    # non-negative lasso fit
    res_theta_j_est <- theta_j_estimation(bj = bj,
                                          B = B,
                                          cj = cj,
                                          eta_j = eta_j,
                                          G = G,
                                          interaction_term = interaction_term,
                                          Yj = Yj,
                                          adj_col = adj_col)
    theta_j<- res_theta_j_est$theta_j

    # At each iteration, truncate small components of theta_j to 0
    theta_j[theta_j < 0.01] <- 0

    # check improvement in this iteration
    improvement <- sqrt(sum((theta_j - theta_j_prev)^2)) / sqrt(sum(theta_j_prev^2))

    if (is.na(improvement)){  # happens when `theta_j` is all zero.
      # warning("improvement is NA, iteration terminated")
      # improvement <- 0  # break the loop by the if statement below
      warning("improvement is NA")
      # warning(c(adj_col))  # debug
      # warning(c(theta_j))  # debug
      improvement <- 2*tol  # continue
    }

    # check stopping criteria
    if (improvement < tol) {break}
  }

  res_j <- list(theta_j = theta_j,
                bj = bj,
                cj = cj,
                eta_j = eta_j)
  return (res_j)
}


.prune_network_step <- function(network_bef,
                                prune_thres = 0.05,
                                var_to_prune = NULL,
                                Y_list,
                                yy_smth_list,
                                obs_time_list,
                                tt,
                                kernel,
                                kernel_params_list,
                                interaction_term,
                                Sigma_k_kl_list,
                                kk_array_list,
                                theta_initial_list = NULL,
                                max_iter = 5,
                                tol = 0.001,
                                verbose = 0){
  # One step (depth) of `prune_network`.
  # Prune the current network by deleting the least important (only one) edge for each variable.
  # `network_bef`: current network to prune. Variables in the 1st column affect variables in the 2nd column. If (i,j)th entry is 1, then variable i affects variable j and one row of `all_edges_idx` is (i, j).
  # `var_to_prune`: a numeric vector that specifies indices if variables to prune (we skip variable j if its affecting variables are all important). Set to NULL to prune for all variables.

  # global quantities
  p <- dim(Y_list[[1]])[2]
  K <- length(Y_list)
  n_vec <- sapply(1:K, function(k){dim(Y_list[[k]])[1]})  # number of time points of observation for each cell
  len <- dim(yy_smth_list[[1]])[1]
  delta <- 1/len
  obs_idx_list <- lapply(1:K, function(k){.map_to(obs_time_list[[k]], tt)})
  # Y_smth_list <- lapply(1:K, function(k){yy_smth_list[[k]][obs_idx_list[[k]],]})  # the smoothed trajectories from step 1. It is `yy_smth` evaluated on `obs_time`.

  B_list <- lapply(1:K, function(k){.construct_B(obs_time_list[[k]])})

  # input check
  if (!(is.null(var_to_prune) | all(var_to_prune %in% (1:p)))){stop("var_to_prune should be a vector containing values in 1,...,p")}
  if (is.null(var_to_prune)){var_to_prune <- 1:p}

  if (verbose > 0) {cat("variables to prune:", paste(var_to_prune, collapse = ", "), "\n")}

  # parallel computing config
  parallel <- TRUE
  if (parallel) {
    detectCores <- parallel::detectCores()
    mc.cores <- ifelse(detectCores > 1, detectCores - 1, 1)
  } else {  # sequential computation
    mc.cores <- 1  # This effectively disables parallel execution in mclapply
  }

  # initialization of theta_j's
  if (is.null(theta_initial_list)){
    num_row <- ifelse(interaction_term, yes = p^2, no = p)
    theta_initial <- matrix(1, nrow = num_row, ncol = p)  # p^2 x p or p x p
    theta_initial_list <- lapply(1:K, function(k){theta_initial})
  }

  # compute recovery metrics using the network before current step of pruning
  # re-fit kernel ODE using the network before pruning and compute the recovery metrics
  Y_est_bef_prune_list <- parallel::mclapply(1:K,
                                             FUN = function(k){
                                               res_KODE_refit <- kernelODE_step2(Y = Y_list[[k]],
                                                                                 obs_time = obs_time_list[[k]],
                                                                                 yy_smth = yy_smth_list[[k]],
                                                                                 tt = tt,
                                                                                 kernel = kernel,
                                                                                 kernel_params = kernel_params_list[[k]],
                                                                                 interaction_term = interaction_term,
                                                                                 theta_initial = theta_initial_list[[k]],
                                                                                 adj_matrix = network_bef,  # network before current pruning
                                                                                 nzero_thres = NULL,
                                                                                 tol = 0.001,
                                                                                 max_iter = 5,
                                                                                 verbose = 0)
                                               # recover trajectory
                                               Y_est_pruned <- sapply(1:ncol(Y_list[[k]]), function(j){
                                                 res_traj_j <- evaluate_Fj(bj = res_KODE_refit$res_bj[j],
                                                                           cj = res_KODE_refit$res_cj[,j],
                                                                           interaction_term = interaction_term,
                                                                           kernel = kernel,
                                                                           kernel_params = kernel_params_list[[k]],
                                                                           kk_array = kk_array_list[[k]],
                                                                           obs_time = res_KODE_refit$config$obs_time,
                                                                           theta_j = res_KODE_refit$res_theta[,j],
                                                                           tt = tt,
                                                                           Yj = Y_list[[k]][,j],
                                                                           yy_smth = yy_smth_list[[k]])
                                                 res_traj_j$yy_est[obs_idx_list[[k]]]
                                                 })
                                               return (Y_est_pruned)
                                               },
                                             mc.cores = mc.cores)
  metrics_bef_prune <- parallel::mclapply(1:K,
                                          FUN = function(k){
                                            assess_recov_traj(Y = Y_list[[k]],
                                                              Y_est = Y_est_bef_prune_list[[k]])
                                            },
                                          mc.cores = mc.cores)

  # pruning procedure and record the results
  network_aft <- network_bef
  prune_edge_idx <- matrix(NA, nrow = 0, ncol = 2)
  colnames(prune_edge_idx) <- c("from", "to")

  mat_temp <- network_bef
  mat_temp[,-var_to_prune] <- NA  # NA indicates either that variable column is not pruned or an edge is not selected in the input network to prune
  mat_temp[network_bef == 0] <- NA
  R2_mat_multi <- lapply(1:K, function(k){mat_temp})  # R2 (% variation explained) by each selected edge (in the input network) for each cell. Set to be NA for non-selected edges.
  R2_mat_multi <- simplify2array(R2_mat_multi)  # (p, p, K)
  R2_mat <- mat_temp  # (p, p), R2 of each selected edge (in the input network) averaged over all cells.
  # To see how much variable j (rather than an edge) explains the variation, simply take the jth row of `R2_mat` (i.e. how much does variable j explains each variable)

  for (j in var_to_prune){
    # For each variable, prune each selected edge separately to assess its influence on trajectory reconstruction
    affecting_var_idx <- which(network_bef[,j] == 1)    # obtain indices of the selected edges for this variable
    if (length(affecting_var_idx) == 0){next}  # no variables affect variable j

    for (var_idx in affecting_var_idx){
      adj_col <- network_bef[,j]
      adj_col[var_idx] <- 0  # delete this edge

      res_j_edge_list <- parallel::mclapply(1:K,
                                             FUN = function(k) {
                                               # refit Fj using the network with the specific edge pruned
                                               res_j <- .kernelODE_step2_single_variable(adj_col = adj_col,  # adj_col is the pruned column
                                                                                         B = B_list[[k]],
                                                                                         interaction_term = interaction_term,
                                                                                         Sigma_k_kl = Sigma_k_kl_list[[k]],
                                                                                         Yj = Y_list[[k]][,j],
                                                                                         theta_j_initial = theta_initial_list[[k]][,j],
                                                                                         max_iter = max_iter,
                                                                                         tol = tol)
                                               # evaluate Fj and reconstruct the trajectory
                                               res_Fj_eval <- evaluate_Fj(bj = res_j$bj,
                                                                          cj = res_j$cj,
                                                                          interaction_term = interaction_term,
                                                                          kernel = kernel,
                                                                          kernel_params = kernel_params_list[[k]],
                                                                          kk_array = kk_array_list[[k]],
                                                                          obs_time = obs_time_list[[k]],
                                                                          theta_j = res_j$theta_j,
                                                                          tt = tt,
                                                                          Yj = Y_list[[k]][,j],
                                                                          yy_smth = yy_smth_list[[k]])
                                               Yj_est_new <- res_Fj_eval$yy_est[obs_idx_list[[k]]]
                                               Yj <- Y_list[[k]][,j]

                                               # compute the R2 using the network with this edge deleted
                                               MS_est_j <- mean((Yj - Yj_est_new)^2)
                                               MS_total_j <- mean((Yj - mean(Yj))^2)
                                               R2_est_j <- max(0, 1 - MS_est_j/MS_total_j)  # R2 of the new reconstructed variable j trajectory with var_idx pruned
                                               R2_edge <- metrics_bef_prune[[k]][["R2_per_var_vec"]][j] - R2_est_j  # % variation explained of variable j traj by this edge for this cell, w.r.t. the latest network

                                               return (list(R2_edge = R2_edge,
                                                            Yj_est_new = Yj_est_new))
                                             },
                                             mc.cores = mc.cores)

      R2_j_edge_multi <- sapply(1:K, function(k){res_j_edge_list[[k]][["R2_edge"]]})  # of length K
      if (length(R2_j_edge_multi) != K) {stop("R2_j_edge_multi should be of length K.")}

      R2_mat_multi[var_idx, j, ] <- R2_j_edge_multi
      R2_mat[var_idx, j] <- mean(R2_j_edge_multi)
    }

    # prune the least important edge for this variable
    R2_j <- R2_mat[,j]
    R2_j_min <- min(R2_j, na.rm = T)  # the least important selected regulator

    if (R2_j_min < prune_thres) {
      # at least one edge can be pruned. Prune the least important regulator
      row_idx_temp <- which(R2_j == R2_j_min)[1]
      network_aft[row_idx_temp, j] <- 0
      prune_edge_idx <- rbind(prune_edge_idx, c(row_idx_temp, j))
    } else {
      # all affecting variables are important. No pruning is done.
    }

    if (verbose > 0) {cat("*", ifelse(j == var_to_prune[length(var_to_prune)], "\n", ""), sep = "")}
  }

  # prompt anomalies or interesting patterns
  neg_edge_idx <- which(R2_mat < 0, arr.ind = T)  # abnormal results (deleting this edge results in a larger % variation explained in the corresponding trajectory)
  if (nrow(neg_edge_idx) > 0) {
    cat("Note: R2 increases after deleting each of the following edges (row - col): ",
        paste(neg_edge_idx[,1], neg_edge_idx[,2], sep = "-", collapse = ", "), "\n")
  }

  # If no pruning is done, prune_edge_idx is an empty matrix with dimension (0,2).
  if (nrow(prune_edge_idx) == 0) {cat("No pruning is done. Network is already pruned to the simplest.\n")}

  var_pruned <- unique(prune_edge_idx[,2])

  return (list(network_bef = network_bef,
               network_aft = network_aft,
               prune_edge_idx = prune_edge_idx,
               var_pruned = var_pruned,
               R2_mat = R2_mat,
               R2_mat_multi = R2_mat_multi))
}




#' Title
#'
#' @param network_original
#' @param prune_thres
#' @param depth
#' @param eval_edge_R2_pruned
#' @param Y_list
#' @param yy_smth_list
#' @param obs_time_list
#' @param tt
#' @param kernel
#' @param kernel_params_list
#' @param interaction_term
#' @param theta_initial_list
#' @param max_iter
#' @param tol
#' @param verbose
#'
#' @returns
#' @export
#'
#' @examples
prune_network <- function(network_original,
                          prune_thres = 0.05,
                          depth = NULL,
                          eval_edge_R2_pruned = TRUE,  # whether to evaluate each edge in the pruned network
                          Y_list,
                          yy_smth_list,
                          obs_time_list,
                          tt,
                          kernel,
                          kernel_params_list,
                          interaction_term,
                          theta_initial_list = NULL,
                          max_iter = 5,
                          tol = 0.001,
                          verbose = 0){
  # Prune the estimated network from kernel ODE by deleting only one edge for each variable (i.e. each column).
  # The procedure assesses the % variation explained by each selected edge on the corresponding variable (averaged over all cells) and prune all edges that explains less than `prune_thres` of variation.
  # This is an iterative process. Each step uses `.prune_network_step`.
  # `network_original` is the network obtained from kernel ODE that we wish to prune.
  # `prune_thres` is the threshold for pruning.
  # `depth` is the maximum number of selected edges that can be pruned for each variable. Set `NULL` to prune until each selected edge explains at least `prune_thres` variation of the variable trajectory.
  # `Y_list` is a list of K matrices each being (n, p) and contains observed trajectories for all cells (note that n can vary for cells, depending on their birth time).
  # `yy_smth_list` is a list of K matrices each being (len, p). len=1/delta is the same across all cells.
  # `kernel` is just a character shared for all cells.
  # `kernel_params_list` is a list of K lists, giving the kernel parameters (e.g. bandwidth) for each cell.
  # `theta_initial_list` is a list of K matrices, each of shape (p^2, p) or (p, p) and contains the initial values of theta_j's (in its columns) for all cells in the iterative optimization.
  # We suggest setting them as the estimates `res_theta_j` obtained from kernel ODE fit of each cell in order to facilitate computation. By default. it is set to be all ones.
  # ONLY WORKS WHEN `interaction_term` is FALSE.

  if (length(kernel_params_list) == 1) {kernel_params_list <- replicate(length(Y_list), kernel_params_list[[1]], simplify = FALSE)}  # if only one kernel parameter set is specified, it is used for all cells
  if (length(obs_time_list) == 1) {obs_time_list <- replicate(length(Y_list), obs_time_list[[1]], simplify = FALSE)}  # similarly for obs_time

  if (verbose > 0){
    cat("\n------ Pruning network ------\n")
    cat("--- Input network to be pruned ---\n")
    print(network_original)
    cat("\n")
  }

  # global quantities
  K <- length(Y_list)
  p <- dim(Y_list[[1]])[2]
  n_vec <- sapply(1:K, function(k){dim(Y_list[[k]])[1]})  # number of time points of observation for each cell
  len <- dim(yy_smth_list[[1]])[1]
  delta <- 1/len
  obs_idx_list <- lapply(1:K, function(k){.map_to(obs_time_list[[k]], tt)})
  # Y_smth_list <- lapply(1:K, function(k){yy_smth_list[[k]][obs_idx_list[[k]],]})  # the smoothed trajectories from step 1. It is `yy_smth` evaluated on `obs_time`.

  # set depth
  if (is.null(depth)){depth <- p}  # the largest possible depth

  # parallel computing
  parallel <- TRUE
  if (parallel) {
    detectCores <- parallel::detectCores()
    mc.cores <- ifelse(detectCores > 1, detectCores - 1, 1)
  } else {  # sequential computation
    mc.cores <- 1  # This effectively disables parallel execution in mclapply
  }

  # recyclable quantities
  B_list <- lapply(1:K, function(k){.construct_B(obs_time_list[[k]])})
  Sigma_k_kl_list <- parallel::mclapply(1:K,
                                   FUN = function(k) {
                                     .construct_Sigma_k_kl(interaction_term = interaction_term,
                                                           kernel = kernel,
                                                           kernel_params = kernel_params_list[[k]],
                                                           obs_time = obs_time_list[[k]],
                                                           tt = tt,
                                                           yy_smth = yy_smth_list[[k]])
                                   },
                                   mc.cores = mc.cores)  # note dim of Sigma_k_kl may vary across cells, since different cells have different #obs tps
  kk_array_list <- parallel::mclapply(1:K,
                                       FUN = function(k){
                                         kk_array <- .construct_kk_array(Y1 = yy_smth_list[[k]],
                                                                         Y2 = yy_smth_list[[k]],
                                                                         interaction_term = interaction_term,
                                                                         kernel = kernel,
                                                                         kernel_params = kernel_params_list[[k]])
                                       },
                                       mc.cores = mc.cores)  # analog of Sigma_k_kl but on `tt`. It has the same dimension (len, len, p) or (len, len, p^2) for all cells.


  ### Assess the recovered trajectories from the unpruned (original) network ###
  # a list of K lists each containing the MS and R2 metrics for each of the p variables
  # re-fit kernel ODE using the unpruned network
  # TODO: Pack this step into a public function
  if (verbose > 0){cat("\n------ Computing unpruned network metrics ------\n")}
  Y_est_list <- parallel::mclapply(1:K,
                                   FUN = function(k){
                                     res_KODE_refit <- kernelODE_step2(Y = Y_list[[k]],
                                                                       obs_time = obs_time_list[[k]],
                                                                       yy_smth = yy_smth_list[[k]],
                                                                       tt = tt,
                                                                       kernel = kernel,
                                                                       kernel_params = kernel_params_list[[k]],
                                                                       interaction_term = interaction_term,
                                                                       theta_initial = theta_initial_list[[k]],
                                                                       adj_matrix = network_original,  # unpruned network
                                                                       nzero_thres = NULL,
                                                                       tol = 0.001,
                                                                       max_iter = 5,
                                                                       verbose = 0)
                                     # recover trajectory
                                     Y_est_unpruned <- sapply(1:p,
                                                              FUN = function(j){
                                                                res_traj_j <- evaluate_Fj(bj = res_KODE_refit$res_bj[j],
                                                                                          cj = res_KODE_refit$res_cj[,j],
                                                                                          interaction_term = interaction_term,
                                                                                          kernel = kernel,
                                                                                          kernel_params = kernel_params_list[[k]],
                                                                                          kk_array = kk_array_list[[k]],
                                                                                          obs_time = res_KODE_refit$config$obs_time,
                                                                                          theta_j = res_KODE_refit$res_theta[,j],
                                                                                          tt = tt,
                                                                                          Yj = Y_list[[k]][,j],
                                                                                          yy_smth = yy_smth_list[[k]])
                                                                res_traj_j$yy_est[obs_idx_list[[k]]]
                                                                })
                                     return (Y_est_unpruned)
                                     },
                                   mc.cores = mc.cores)
  # compute the metrics of the unpruned network
  metrics_unpruned <- parallel::mclapply(1:K,
                                               FUN = function(k){
                                                 assess_recov_traj(Y = Y_list[[k]],
                                                                   Y_est = Y_est_list[[k]])
                                                 },
                                               mc.cores = mc.cores)


  ### Pruning process ###
  res_prune_path <- list()  # record the path of the pruning

  network_bef <- network_original
  var_to_prune <- 1:p

  for (depth_idx in 1:depth){
    if (verbose > 0) {cat("\n------ Depth", depth_idx, "------\n")}

    res_prune_current <- .prune_network_step(network_bef = network_bef,
                                             prune_thres = prune_thres,
                                             var_to_prune = var_to_prune,
                                             Y_list = Y_list,
                                             yy_smth_list = yy_smth_list,
                                             obs_time_list = obs_time_list,
                                             tt = tt,
                                             kernel = kernel,
                                             kernel_params_list = kernel_params_list,
                                             interaction_term = interaction_term,
                                             Sigma_k_kl_list = Sigma_k_kl_list,
                                             kk_array_list = kk_array_list,
                                             theta_initial_list = theta_initial_list,
                                             max_iter = max_iter,
                                             tol = tol,
                                             verbose = verbose)
    R2_mat <- res_prune_current$R2_mat
    prune_edge_idx <- res_prune_current$prune_edge_idx
    network_aft <- res_prune_current$network_aft

    if (verbose > 0){
      cat("\n--- R2 of each selected edge ---\n")
      print(round(R2_mat, 2))
      cat("\n")

      cat("The following edges are pruned:",
          ifelse(nrow(prune_edge_idx) == 0,
                 yes = "None",
                 no = paste(prune_edge_idx[,1], prune_edge_idx[,2], sep = "-", collapse = ", ")))

      cat("\n--- network after pruning ---\n")
      print(network_aft)
      cat("\n")
    }

    # record into the path
    res_prune_path[[depth_idx]] <- list(depth_idx = depth_idx,
                                        network_bef = res_prune_current$network_bef,
                                        network_aft = res_prune_current$network_aft,
                                        prune_edge_idx = res_prune_current$prune_edge_idx,
                                        var_pruned = res_prune_current$var_pruned,
                                        R2_mat = res_prune_current$R2_mat,
                                        R2_mat_multi = res_prune_current$R2_mat_multi)

    # update for next iteration
    var_to_prune <- unique(prune_edge_idx[,2])  # if a variable is not pruned at current iteration (i.e. not in `prune_edge_idx[,2]`), then it is already pruned to the simplest
    network_bef <- network_aft

    if (length(var_to_prune) == 0) {
      if (verbose > 0){cat("\nPruning finished for all variables at depth", depth_idx, "\n")}
      break
    }
  }

  # save the results
  network_pruned <- network_aft


  ### re-assess importance (improvement in R2) of each edge in the final pruned network ###
  if (!eval_edge_R2_pruned){
    R2_mat_pruned <- NULL
    R2_mat_multi_pruned <- NULL
  } else {
    if (verbose > 0) {cat("\n--- Assessing R2 of each selected edge in the pruned network ---\n")}

    R2_mat_pruned <- matrix(NA, nrow = dim(network_pruned)[1], ncol = dim(network_pruned)[2])
    R2_mat_multi_pruned <- array(NA, c(dim(network_pruned)[1], dim(network_pruned)[2], K))

    for (j in 1:p){
      # For each variable, prune each selected edge separately to assess its influence on trajectory reconstruction
      affecting_var_idx <- which(network_pruned[,j] == 1)    # obtain indices of the selected edges for this variable
      if (length(affecting_var_idx) == 0){next}  # no variables affect variable j

      for (var_idx in affecting_var_idx){
        adj_col <- network_pruned[,j]
        adj_col[var_idx] <- 0  # delete this edge

        res_j_edge_list <- parallel::mclapply(1:K,
                                               FUN = function(k) {
                                                 # refit Fj using the network with the specific edge pruned
                                                 res_j <- .kernelODE_step2_single_variable(adj_col = adj_col,  # adj_col is the pruned column
                                                                                           B = B_list[[k]],
                                                                                           interaction_term = interaction_term,
                                                                                           Sigma_k_kl = Sigma_k_kl_list[[k]],
                                                                                           Yj = Y_list[[k]][,j],
                                                                                           theta_j_initial = theta_initial_list[[k]][,j],
                                                                                           max_iter = max_iter,
                                                                                           tol = tol)
                                                 # evaluate Fj and reconstruct the trajectory
                                                 res_Fj_eval <- evaluate_Fj(bj = res_j$bj,
                                                                            cj = res_j$cj,
                                                                            interaction_term = interaction_term,
                                                                            kernel = kernel,
                                                                            kernel_params = kernel_params_list[[k]],
                                                                            kk_array = kk_array_list[[k]],
                                                                            obs_time = obs_time_list[[k]],
                                                                            theta_j = res_j$theta_j,
                                                                            tt = tt,
                                                                            Yj = Y_list[[k]][,j],
                                                                            yy_smth = yy_smth_list[[k]])
                                                 Yj_est_new <- res_Fj_eval$yy_est[obs_idx_list[[k]]]
                                                 Yj <- Y_list[[k]][,j]

                                                 MS_est_j <- mean((Yj - Yj_est_new)^2)  # mean squared residuals of est. trajectory w.r.t. observations (variance in our estimated trajectory / "recovered signals")
                                                 MS_total_j <- mean(Yj^2)  # mean squared total of the observations. Note that `Y_list` has centered columns.
                                                 R2_est_j <- max(0, 1 - MS_est_j/MS_total_j)  # R2 of the new reconstructed variable j trajectory with var_idx pruned
                                                 R2_edge <- metrics_unpruned[[k]][["R2_per_var_vec"]][j] - R2_est_j  # % variation explained of variable j traj by this edge for this cell, w.r.t. the original unpruned network

                                                 return (list(R2_edge = R2_edge,
                                                              Yj_est_new = Yj_est_new))
                                               },
                                               mc.cores = mc.cores)

        R2_j_edge_multi <- sapply(1:K, function(k){res_j_edge_list[[k]][["R2_edge"]]})  # of length K
        if (length(R2_j_edge_multi) != K) {stop("R2_j_edge_multi should be of length K.")}

        R2_mat_multi_pruned[var_idx, j, ] <- R2_j_edge_multi
        R2_mat_pruned[var_idx, j] <- mean(R2_j_edge_multi)
      }
      if (verbose > 0) {cat("*", ifelse(j == p, "\n", ""), sep = "")}
    }
  }


  # print results
  if (verbose > 0){
    cat("\n------ Pruning finished ------\n")
    cat("\n------ Pruned network ------\n")
    print(network_pruned)
    cat("\n--- R2 of each selected edge in the pruned network ---\n")
    print(round(R2_mat_pruned, 2))
  }

  res_prune_network <- list(res_prune_path = res_prune_path,
                            network_pruned = network_pruned,
                            R2_mat_pruned = R2_mat_pruned,
                            R2_mat_multi_pruned = R2_mat_multi_pruned,
                            ### metrics and trajectories of old and new trajectories ###
                            config = list(
                              network_original = network_original,
                              prune_thres = prune_thres,
                              depth = depth,
                              Y_list = Y_list,
                              Y_est_list = Y_est_list,
                              yy_smth_list = yy_smth_list,
                              obs_time_list = obs_time_list,
                              tt = tt,
                              kernel = kernel,
                              kernel_params_list = kernel_params_list,
                              interaction_term = interaction_term,
                              theta_initial_list = theta_initial_list,
                              max_iter = max_iter,
                              tol = tol))

  return (res_prune_network)
}









