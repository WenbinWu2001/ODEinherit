## ---- Exported functions ----
#' Prune a Network Estimated by Kernel ODE
#'
#' @description
#' Iteratively removes edges from an estimated network by evaluating the
#' percentage of variation each edge explains in the target variable's
#' trajectory. Edges explaining less than `prune_thres` are pruned.
#' The process repeats up to `depth` times or until no further edges can be pruned.
#' Applicable only when `interaction_term = FALSE`.
#'
#' Supports pruning using multiple samples (cells) jointly, where edge contributions (\eqn{R^2}) are averaged across samples.
#'
#' @inheritParams kernelODE_step2
#' @param network_original The original network (adjacency matrix) to prune.
#' @param prune_thres Threshold for pruning. Edges explaining less than this
#'   proportion of variation are removed.
#' @param depth Maximum number of pruning steps per variable. Use `NULL` to
#'   continue pruning until all edges meet the threshold.
#' @param eval_edge_R2_pruned A logical value. If `TRUE`, evaluates the contribution
#'   (\eqn{R^2}) of each remaining edge in the pruned network.
#' @param Y_list A list of observed trajectories (`Y`), each corresponding to a sample (cell). See [kernelODE_step2()] for the format of `Y`.
#' @param yy_smth_list A list of smoothed trajectories (`yy_smth`).
#' @param obs_time_list A list of standardized observation time points (`obs_time`).
#' @param kernel_params_list A list of kernel parameter lists (`kernel_params`).
#' @param theta_initial_list A list of initial \eqn{\theta_j} matrices for iterative optimization (`theta_initial`).
#' @param parallel A logical value. If `TRUE`, use parallel computing (Linux/macOS only). Parallelization is applied over variables and samples.
#' @param mc.cores An integer specifying the number of cores for parallel computing. Defaults to all but one available cores.
#'
#' @returns A list with components:
#' \describe{
#'   \item{res_prune_path}{History of pruning steps.}
#'   \item{network_pruned}{Final pruned network.}
#'   \item{R2_avg_mat_pruned}{\eqn{R^2} contribution per edge in the pruned network, averaged across all samples.}
#'   \item{R2_multi_arr_pruned}{Per-sample edge \eqn{R^2} contributions.}
#'   \item{config}{List of pruning configuration and inputs.}
#' }
#'
#' @export
prune_network <- function(network_original,
                          prune_thres = 0.05,
                          depth = NULL,
                          eval_edge_R2_pruned = FALSE,
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
                          parallel = FALSE,
                          mc.cores = NULL,
                          verbose = 0){
  if (length(kernel_params_list) == 1) {kernel_params_list <- replicate(length(Y_list), kernel_params_list[[1]], simplify = FALSE)}  # if only one kernel parameter set is specified, it is used for all cells
  if (length(obs_time_list) == 1) {obs_time_list <- replicate(length(Y_list), obs_time_list[[1]], simplify = FALSE)}  # similarly for obs_time

  if (verbose > 0){
    cat("\n------ Pruning network ------\n")
    cat("--- Input network to be pruned ---\n")
    print(network_original)
    cat("\n")
  }

  # config parallel computing (Parallel computing will be employed on (1) pruning edges for each variable, and (2) evaluating the goodness of fit over the samples.)
  if (parallel) {
    max_num_cores <- parallel::detectCores()
    if (is.null(mc.cores)){
      mc.cores <- ifelse(max_num_cores > 1, max_num_cores-1, 1)
    }
    if (mc.cores > max_num_cores){
      mc.cores <- max_num_cores
    }

    mc.cores <- floor(mc.cores)  # in case a decimal number is given
  } else {  # sequential computation
    mc.cores <- 1  # This effectively disables parallel execution in mclapply
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

  # recyclable quantities
  B_list <- lapply(1:K, function(k){.construct_B(obs_time_list[[k]])})
  Sigma_k_kl_list <- lapply(1:K,
                            FUN = function(k) {
                              .construct_Sigma_k_kl(interaction_term = interaction_term,
                                                    kernel = kernel,
                                                    kernel_params = kernel_params_list[[k]],
                                                    obs_time = obs_time_list[[k]],
                                                    tt = tt,
                                                    yy_smth = yy_smth_list[[k]])
                            })  # note dim of Sigma_k_kl may vary across cells, since different cells have different #obs tps
  kk_array_list <- lapply(1:K,
                          FUN = function(k){
                            kk_array <- .construct_kk_array(Y1 = yy_smth_list[[k]],
                                                            Y2 = yy_smth_list[[k]],
                                                            interaction_term = interaction_term,
                                                            kernel = kernel,
                                                            kernel_params = kernel_params_list[[k]])
                          })  # analog of Sigma_k_kl but on `tt`. It has the same dimension (len, len, p) or (len, len, p^2) for all cells.


  ### Assess the recovered trajectories from the unpruned (original) network ###
  # a list of K lists each containing the MS and R2 metrics for each of the p variables
  # re-fit kernel ODE using the unpruned network
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
                                                                       max_iter = max_iter,
                                                                       tol = tol,
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
  metrics_unpruned <- lapply(1:K,
                             FUN = function(k){
                               assess_recov_traj(Y = Y_list[[k]],
                                                 Y_est = Y_est_list[[k]])
                             })


  ### Pruning process ###
  res_prune_path <- list()  # record the path of the pruning

  network_bef <- network_original
  var_to_prune_vec <- 1:p

  for (depth_idx in 1:depth){
    if (verbose > 0) {cat("\n------ Depth", depth_idx, "------\n")}

    res_prune_current <- .prune_network_step(network_bef = network_bef,
                                             prune_thres = prune_thres,
                                             var_to_prune_vec = var_to_prune_vec,
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
                                             parallel = parallel,
                                             mc.cores = mc.cores,
                                             verbose = verbose)
    R2_avg_mat <- res_prune_current$R2_avg_mat
    pruned_edge_mat <- res_prune_current$pruned_edge_mat
    network_aft <- res_prune_current$network_aft

    if (verbose > 0){
      cat("\n--- R2 of each selected edge ---\n")
      print(round(R2_avg_mat, 2))
      cat("\n")

      cat("The following edges are pruned:",
          ifelse(nrow(pruned_edge_mat) == 0,
                 yes = "None",
                 no = paste(pruned_edge_mat[,1], pruned_edge_mat[,2], sep = "-", collapse = ", ")))

      cat("\n--- network after pruning ---\n")
      print(network_aft)
      cat("\n")
    }

    # record into the path
    res_prune_path[[depth_idx]] <- list(depth_idx = depth_idx,
                                        network_bef = res_prune_current$network_bef,
                                        network_aft = res_prune_current$network_aft,
                                        pruned_edge_mat = res_prune_current$pruned_edge_mat,
                                        var_pruned = res_prune_current$var_pruned,
                                        R2_avg_mat = res_prune_current$R2_avg_mat,
                                        R2_multi_arr = res_prune_current$R2_multi_arr)

    # update for next iteration
    var_to_prune_vec <- unique(pruned_edge_mat[,2])  # if a variable is not pruned at current iteration (i.e. not in `pruned_edge_mat[,2]`), then it is already pruned to the simplest
    network_bef <- network_aft

    if (length(var_to_prune_vec) == 0) {
      if (verbose > 0){cat("\nPruning finished for all variables at depth", depth_idx, "\n")}
      break
    }
  }

  # save the results
  network_pruned <- network_aft


  ### re-assess importance (improvement in R2) of each edge in the final pruned network ###
  if (!eval_edge_R2_pruned){
    R2_multi_arr_pruned <- NULL
    R2_avg_mat_pruned <- NULL
  } else {
    if (verbose > 0) {cat("\n--- Assessing R2 of each selected edge in the pruned network ---\n")}

    R2_multi_arr_pruned <- array(NA, c(dim(network_pruned)[1], dim(network_pruned)[2], K))
    R2_avg_mat_pruned <- matrix(NA, nrow = dim(network_pruned)[1], ncol = dim(network_pruned)[2])

    res_edge_eval <- parallel::mclapply(1:p, function(j){
      # For each variable, prune each selected edge separately to assess its influence on trajectory reconstruction
      regulator_var_idx_vec <- which(network_pruned[,j] == 1)    # obtain indices of the selected edges for this variable
      if (length(regulator_var_idx_vec) == 0){
        return (list(NULL))
      }  # no variables affect variable j

      R2_j_multi_mat_pruned <- R2_multi_arr_pruned[,j,,drop = F]  # slice j, an array of dim (p, 1 ,K). This suppress the dimension dropping when K=1.
      R2_j_avg_vec_pruned <- R2_avg_mat_pruned[,j,drop = F]  # slice j, a matrix of dim (p, 1)
      for (var_idx in regulator_var_idx_vec){
        adj_col <- network_pruned[,j]
        adj_col[var_idx] <- 0  # remove this edge

        R2_temp_vec <- sapply(1:K,
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

                                # compute the R2 of the new recovered variable j trajectory with var_idx pruned
                                R2_est_j <- assess_recov_traj(Y = Yj, Y_est = Yj_est_new)$R2
                                # the R2 contribution of this edge (we expect it to be positive, but it can be negative due to overfitting)
                                R2_edge <- metrics_unpruned[[k]][["R2_per_var_vec"]][j] - R2_est_j

                                return (R2_edge)
                              })  # of length K

        # update the R2 results
        R2_j_multi_mat_pruned[var_idx, 1, ] <- R2_temp_vec
        R2_j_avg_vec_pruned[var_idx, 1] <- mean(R2_temp_vec)
      }  # end of iterating through all edges/regulators for this variable

      return (list(R2_j_multi_mat_pruned = R2_j_multi_mat_pruned,
                   R2_j_avg_vec_pruned = R2_j_avg_vec_pruned))
    }, mc.cores = mc.core)

    for (j in 1:p){
      res <- res_edge_eval[[j]]
      if (is.null(res)){next}

      R2_multi_arr_pruned[,j,] <- res$R2_j_multi_mat_pruned
      R2_avg_mat_pruned[,j] <- res$R2_j_avg_vec_pruned
    }
  }  # end of evaluating each selected edge for all varibales


  # print results
  if (verbose > 0){
    cat("\n------ Pruning finished ------\n")
    cat("\n------ Pruned network ------\n")
    print(network_pruned)
    cat("\n--- R2 of each selected edge in the pruned network ---\n")
    print(round(R2_avg_mat_pruned, 2))
  }

  res_prune_network <- list(res_prune_path = res_prune_path,
                            network_pruned = network_pruned,
                            R2_avg_mat_pruned = R2_avg_mat_pruned,
                            R2_multi_arr_pruned = R2_multi_arr_pruned,
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
                              tol = tol,
                              parallel = parallel,
                              mc.cores = mc.cores))

  return (res_prune_network)
}






## ---- Internal helpers (private) ----
.kernelODE_step2_single_variable <- function(adj_col,
                                             B,
                                             interaction_term,
                                             Sigma_k_kl,
                                             Yj,
                                             theta_j_initial = NULL,
                                             max_iter = 5,
                                             tol = 0.001){
  # A dedicated function for running the single-sample kernel ODE on a single variable j during the network pruning.
  # `adj_col` is the jth column of `network_est` with a specific edge removed.
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
    res_fun_est <- .function_estimation(Yj = Yj,
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
    res_theta_j_est <- .theta_j_estimation(bj = bj,
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
      improvement <- 2*tol  # ignore and continue
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
                                var_to_prune_vec = NULL,
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
                                parallel = FALSE,
                                mc.cores = NULL,
                                verbose = 0){
  # One step (depth) of `prune_network`.
  # Prune the current network by deleting the least important (only one) edge for each variable.
  # `network_bef`: current network to prune. Variables in the 1st column affect variables in the 2nd column. If (i,j)th entry is 1, then variable i affects variable j and one row of `all_edges_idx` is (i, j).
  # `var_to_prune_vec`: a numeric vector that specifies indices for which variables their selected edges need pruning (we skip variable j if its affecting variables are all important). Set to NULL to prune for all variables.

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
  if (!(is.null(var_to_prune_vec) | all(var_to_prune_vec %in% (1:p)))){stop("var_to_prune_vec should be a vector containing values in 1,...,p")}
  if (is.null(var_to_prune_vec)){var_to_prune_vec <- 1:p}

  if (verbose > 0) {cat("variables to prune:", paste(var_to_prune_vec, collapse = ", "), "\n")}


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
                                                                                 max_iter = max_iter,
                                                                                 tol = tol,
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
  metrics_bef_prune <- lapply(1:K,
                              FUN = function(k){
                                assess_recov_traj(Y = Y_list[[k]],
                                                  Y_est = Y_est_bef_prune_list[[k]])
                                })

  ### single step pruning (prune one edge/regulator for each variable) ###
  network_aft <- network_bef
  pruned_edge_mat <- matrix(NA, nrow = 0, ncol = 2)
  colnames(pruned_edge_mat) <- c("from", "to")

  mat_temp <- network_bef
  mat_temp[,-var_to_prune_vec] <- NA  # NA indicates either that variable column is not pruned or an edge is not selected in the input network to prune
  mat_temp[network_bef == 0] <- NA
  R2_multi_arr <- lapply(1:K, function(k){mat_temp})  # R2 (% variation explained) by each selected edge (in the input network) for each cell. Set to be NA for non-selected edges.
  R2_multi_arr <- simplify2array(R2_multi_arr)  # (p, p, K)
  R2_avg_mat <- mat_temp  # (p, p), R2 of each selected edge (in the input network) averaged over all cells.
  # To see how much variable j (rather than an edge) explains the variation, simply take the jth row of `R2_avg_mat` (i.e. how much does variable j explains each variable)

  # evaluate each edge by assessing change in R2 due to their removal
  res_edge_eval_single_step <- parallel::mclapply(var_to_prune_vec, function(j){
    # For each variable, prune each selected edge separately to assess its influence on trajectory reconstruction
    regulator_var_idx_vec <- which(network_bef[,j] == 1)    # obtain indices of the selected edges for this variable

    # If no variables affect variable j, continue to the next variable to prune
    if (length(regulator_var_idx_vec) == 0){
      return (list(NULL))
    }

    R2_j_multi_mat <- R2_multi_arr[,j,,drop = F]  # slice j, an array of dim (p, 1 ,K). This suppress the dimension dropping when K=1.
    R2_j_avg_vec <- R2_avg_mat[,j,drop = F]  # slice j, a matrix of dim (p, 1)
    for (var_idx in regulator_var_idx_vec){
      adj_col <- network_bef[,j]
      adj_col[var_idx] <- 0  # remove this edge

      # compute the R2 for each sample after removing this edge
      R2_temp_vec <- sapply(1:K,
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

                              # compute the R2 of the new recovered variable j trajectory with var_idx pruned
                              R2_est_j <- assess_recov_traj(Y = Yj, Y_est = Yj_est_new)$R2
                              # the R2 contribution of this edge (we expect it to be positive, but it can be negative due to overfitting)
                              R2_edge <- metrics_bef_prune[[k]][["R2_per_var_vec"]][j] - R2_est_j

                              return (R2_edge)
                            })  # of length K

      # update the R2 results
      R2_j_multi_mat[var_idx, 1, ] <- R2_temp_vec
      R2_j_avg_vec[var_idx, 1] <- mean(R2_temp_vec)
    }

    # prune the least important edge for this variable
    R2_j_min <- min(c(R2_j_avg_vec), na.rm = T)  # the least important selected regulator

    if (R2_j_min < prune_thres) {
      # at least one edge can be pruned. Prune the least important regulator.
      var_idx_pruned <- which(c(R2_j_avg_vec) == R2_j_min)[1]
      pruned_edge <- c(from = var_idx_pruned, to = j)
    } else {
      # all affecting variables are important. No pruning is done.
      pruned_edge <- c(NA, j)
    }

    return (list(var_to_prune = j,
                 pruned_edge = pruned_edge,
                 R2_j_multi_mat = R2_j_multi_mat,
                 R2_j_avg_vec = R2_j_avg_vec
    ))
  }, mc.cores = mc.cores)

  # prune the network
  for (res in res_edge_eval_single_step){
    if (is.null(res)) {next}

    j <- res$var_to_prune
    pruned_edge <- res$pruned_edge
    if (any(is.na(pruned_edge))){next}  # no edge to prune

    network_aft[pruned_edge[1], pruned_edge[2]] <- 0
    pruned_edge_mat <- rbind(pruned_edge_mat, pruned_edge)
    R2_multi_arr[,j,] <- res$R2_j_multi_mat
    R2_avg_mat[,j] <- res$R2_j_avg_vec
  }

  # identify edges where removing them increases R2 (possible overfitting)
  neg_edge_mat <- which(R2_avg_mat < 0, arr.ind = TRUE)
  if (nrow(neg_edge_mat) > 0) {
    cat("Note: R2 increases when removing these edges:: ",
        paste(neg_edge_mat[,1], neg_edge_mat[,2], sep = "->", collapse = ", "), "\n")
  }

  # If no pruning is done, pruned_edge_mat is an empty matrix with dimension (0,2).
  if (nrow(pruned_edge_mat) == 0) {cat("No pruning is done. Network is already pruned to the simplest.\n")}

  var_pruned <- unique(pruned_edge_mat[,2])

  return (list(network_bef = network_bef,
               network_aft = network_aft,
               pruned_edge_mat = pruned_edge_mat,
               var_pruned = var_pruned,
               R2_avg_mat = R2_avg_mat,
               R2_multi_arr = R2_multi_arr))
}



