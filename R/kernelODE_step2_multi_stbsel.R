kernelODE_step2_multi <- function(Y,
                                  obs_time,
                                  yy_smth,
                                  tt,
                                  kernel,
                                  kernel_params,
                                  interaction_term = FALSE,
                                  theta_initial = NULL,
                                  # adj_matrix = NULL,
                                  nzero_thres = NULL,
                                  stability_thres = 0.6,
                                  tol = 0.001,
                                  max_iter = 10,
                                  verbose = 0){
  ## kernel ODE Step 2 (multi-sample version): Iterative optimization algorithm with stability selection in network inference (SUPPORTED ONLY WHEN `interaction_term` is FALSE)
  ## `Y` with shape (n, p, K) is the an array of p variables observed at n time points for K cells. The time points of observations are assumed to be [0,1]. All series are assumed to be observed on a common grid of time points.
  ## `yy_smth` with shape (len, p, K) is the an array of p variables observed at n time points for K cells. The time points of observations are assumed to be [0,1]. All series are assumed to be observed on a common grid of time points.
  ## all cell series in `Y` and `yy_smth` should be centered.

  # if only one kernel parameter set is specified, it is used for all variables
  if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y), kernel_params[[1]], simplify = FALSE)}

  # just in case `Y` and `yy_smth` are not centered at 0
  Y <- apply(Y, c(2,3), scale, center = TRUE, scale = FALSE)
  # DO NOT CENTER yy_smth!

  if (interaction_term) {warning("Stability selection NOT supported when `interaction_term` is TRUE.")}

  # parallel computing
  parallel <- TRUE
  if (parallel) {
    detectCores <- parallel::detectCores()
    mc.cores <- ifelse(detectCores > 1, detectCores - 1, 1)
  } else {  # sequential computation
    mc.cores <- 1  # This effectively disables parallel execution in mclapply
  }


  n <- dim(Y)[1]
  p <- dim(Y)[2]
  K <- dim(Y)[3]


  # recyclable quantities
  B <- .construct_B(obs_time)
  qr_decomp <- qr(B)
  Q <- qr.Q(qr_decomp, complete = TRUE)
  Q1 <- matrix(Q[,1], nrow = n, ncol = 1)
  Q2 <- Q[,2:n]
  R <- qr.R(qr_decomp)[1,1]


  # Construct Sigma's for each cell
  if (verbose > 0) {cat("-------- constructing Sigma_k_kl for each cell --------\n")}

  ## without parallel computing
  # Sigma_k_kl <- array(NA, c(n, n, p, K))
  # for (k in 1:K){
  #   Sigma_k_kl[,,,k] <- .construct_Sigma_k_kl(interaction_term = interaction_term,
  #                                             kernel = kernel,
  #                                             kernel_params = kernel_params,
  #                                             obs_time = obs_time,
  #                                             tt = tt,
  #                                             yy_smth = yy_smth[,,k])
  # }

  ## with parallel computing
  Sigma_k_kl <- parallel::mclapply(1:K,
                                   FUN = function(k) {
                                     .construct_Sigma_k_kl(interaction_term = interaction_term,
                                                           kernel = kernel,
                                                           kernel_params = kernel_params,
                                                           obs_time = obs_time,
                                                           tt = tt,
                                                           yy_smth = yy_smth[,,k])
                                   },
                                   mc.cores = mc.cores)
  Sigma_k_kl <- simplify2array(Sigma_k_kl)  # of dim (n, n, p, K)

  if ((!interaction_term) & any(dim(Sigma_k_kl) != c(n, n, p, K))) stop("Sigma_k_kl should have dimension (n, n, p, K): `interaction_term` is FALSE")
  if (interaction_term & any(dim(Sigma_k_kl) != c(n, n, p^2, K))) stop("Sigma_k_kl should have dimension (n, n, p^2, K): `interaction_term` is TRUE.")
  if (any(is.na(Sigma_k_kl))) stop("NA occurs in Sigma_k_kl")


  # Iterative optimization
  # initialization: if no specified theta_initial, then initialize to 1 for all entries.
  if (is.null(theta_initial)){
    num_row <- ifelse(interaction_term, yes = p^2, no = p)
    theta_initial <- array(1, dim = c(num_row, p, K))
  }
  res_theta_multi <- theta_initial  # of shape (p^2, p, K) if `interaction_term` is TRUE or (p, p, K) if `interaction_term` is FALSE


  # iteration begins
  for (num_iter in 1:max_iter){
    if (verbose > 0) {cat("\n-------- iteration", num_iter, " --------\n")}

    res_theta_multi_prev <- res_theta_multi

    # given theta_j, estimate F_j
    # infer network using stability selection by pooling the adjacency matrices from all cells
    if (verbose > 0) {cat("-------- estimating Fj's and network by stabilty selection --------\n")}
    res_multi <- .stability_selection(B = B,
                                      interaction_term = interaction_term,
                                      mc.cores = mc.cores,
                                      nzero_thres = nzero_thres,
                                      Q1 = Q1,
                                      Q2 = Q2,
                                      res_theta_multi = res_theta_multi,
                                      R = R,
                                      stability_thres = stability_thres,
                                      Sigma_k_kl = Sigma_k_kl,
                                      verbose = verbose,
                                      Y = Y)
    adj_matrix <- res_multi$adj_matrix  # (p, p)
    adj_matrix_multi <- res_multi$adj_matrix_multi  # (p, p, K)
    # res_bj_multi <- res_multi$res_bj_multi  # (p, K)
    # res_cj_multi <- res_multi$res_cj_multi  # (n, p, K)
    # res_best_eta_multi <- res_multi$res_best_eta_multi  # (p, K)

    # # fit theta_j for each cell k using the network from stability selection
    # if (verbose > 0) {cat("-------- estimating theta_j's for each cell --------\n")}
    # res_theta_est <- parallel::mclapply(1:K,
    #                                       FUN = function(k) {
    #                                         res_theta_single <- res_theta_multi[,,k]  #  (p^2, p) or  (p, p), initialized to store the estimated theta_j's for this cell
    #                                         res_best_kappa_single <- rep(NA, p)
    #                                         for (j in 1:p){
    #                                           Yj <- Y[,j,k]
    #                                           bj <- res_bj_multi[j,k]
    #                                           cj <- res_cj_multi[,j,k]
    #                                           eta_j <- res_best_eta_multi[j,k]
    #                                           G <- .construct_G(cj = cj,
    #                                                             interaction_term = interaction_term,
    #                                                             Sigma_k_kl = Sigma_k_kl[,,,k])
    #
    #                                           # lasso fit
    #                                           adj_col <- adj_matrix[,j]  # variable j corresponds to column j of the adjacency matrix. Note that `adj_col` is NULL if `adj_matrix` is NULL.
    #                                           res_theta_j_est <- theta_j_estimation(bj = bj,
    #                                                                                 B = B,
    #                                                                                 cj = cj,
    #                                                                                 eta_j = eta_j,
    #                                                                                 G = G,
    #                                                                                 interaction_term = interaction_term,
    #                                                                                 Yj = Yj,
    #                                                                                 adj_col = adj_col,
    #                                                                                 nzero_thres = nzero_thres)
    #                                           res_theta_single[,j] <- res_theta_j_est$theta_j  # update theta_j
    #                                           res_best_kappa_single[j] <- res_theta_j_est$best_kappa
    #                                         }
    #
    #                                         return (list(res_theta_single = res_theta_single,
    #                                                      res_best_kappa_single = res_best_kappa_single))
    #                                       },
    #                                       mc.cores = mc.cores)
    res_theta_multi <- simplify2array(lapply(1:K, function(k){res_theta_est[[k]][["res_theta_single"]]}))  # (p^2, p, K) or (p, p, K)
    res_best_kappa_multi <- simplify2array(lapply(1:K, function(k){res_theta_est[[k]][["res_best_kappa_single"]]}))  # (p, K)

    if (any(dim(res_theta_multi) != c(ifelse(interaction_term, p^2, p), p, K))) {stop("res_theta_multi should have dimension (p, p, K) or (p^2, p, K)")}

    # At each iteration, truncate small components of theta_j to 0
    res_theta_multi[res_theta_multi < 0.01] <- 0

    # check improvement in this iteration
    improvement <- sum((res_theta_multi - res_theta_multi_prev)^2) / sum(res_theta_multi_prev^2)

    if (verbose > 0) {
      cat("theta (for all cells) improvement:", round(improvement, 4), "\n")

      avg_Fnorm_res_theta <- mean(apply(res_theta_multi, 3, function(res_theta){norm(res_theta, type = "F")}))
      cat("average F-norm of current theta for all cells:", round(avg_Fnorm_res_theta, 4), "\n")
    }

    if (verbose > 0) {
      cat("--- inferred network ---\n")
      print(adj_matrix)
      cat("\n")
      cat("--- selection rate ---\n")
      print(round(apply(adj_matrix_multi, c(1,2), mean), 2))
      cat("\n")
    }

    if (is.na(improvement)){  # usually because `norm(res_theta_prev, type = "F")` is 0 (i.e. `res_theta_prev` is all 0) at the first few iterations
      warning("improvement is NA")
      improvement <- 2*tol # continue
    }

    # check stopping criteria
    if (improvement < tol) {break}
  }  # end of iteration loop

  ### CAVEAT ###
  ### DO NOT TAKE THE MEAN OF theta_j, bj, cj across all cells. The underlying math is more complicated than a simple average and may involve non-linearity.
  ### YOU WILL NOT OBTAIN A HIGHER PRECISION ESTIMATION BY SIMPLY TAKING MEAN.
  ### INSTEAD, USE `evaluate_Fj` to obtain the estimated derivative functions Fj for each cell, then take the mean of these Fj's (over all cells) for each variable j.

  if (verbose > 0) {cat("\n -------- finished -------- \n")}


  res_kernelODE <- list(adj_matrix = adj_matrix,   # (p, p)
                        adj_matrix_multi = adj_matrix_multi,  # (p, p, K)
                        res_theta_multi = res_theta_multi, # (p^2, p, K) or (p, p, K)
                        res_best_kappa_multi = res_best_kappa_multi,  # (p, K)
                        res_bj_multi = res_bj_multi,  # (p, K)
                        res_cj_multi = res_cj_multi,  # (n, p, K)
                        res_best_eta_multi = res_best_eta_multi,  # (p, K)
                        num_iter = num_iter,
                        last_improvement = improvement,
                        config = list(
                          interaction_term = interaction_term,
                          kernel = kernel,
                          kernel_params = kernel_params,
                          max_iter = max_iter,
                          nzero_thres = nzero_thres,
                          obs_time = obs_time,
                          stability_thres = stability_thres,
                          theta_initial = theta_initial,
                          tol = tol,
                          tt = tt,
                          Y = Y,
                          yy_smth = yy_smth))


  return(res_kernelODE)
}
