kernelODE_step2_multi_final <- function(Y,
                                        obs_time,
                                        yy_smth,
                                        tt,
                                        kernel,
                                        kernel_params,
                                        interaction_term = FALSE,
                                        theta_initial = NULL,
                                        adj_matrix = NULL,
                                        nzero_thres = NULL,
                                        tol = 0.001,
                                        max_iter = 10,
                                        verbose = 0){
  ## kernel ODE Step 2 (multi-sample version): Iterative optimization algorithm for estimating the (shared) Fj and theta_j using all cells without cell interactions.
  ## Specifically, no cell interactions implies the big Sigma (Kn x Kn) is block-diagonal, with each block being nxn as the big Sigma for the corresponding cell.
  ## (1) In function estimation step, solve bj and cj using all cells assuming no cell interactions (i.e. ).
  ## (2) In theta estimation step, solve theta_j by concatenating row-wise the G's and zj's from all cells and then performing Lasso on this concatenated system.
  ## `Y` with shape (n, p, K) is the an array of p variables observed at n time points for K cells. The time points of observations are assumed to be [0,1]. All series are assumed to be observed on a common grid of time points.
  ## `yy_smth` with shape (len, p, K) is the an array of p variables observed at len time points for K cells. The time points of observations are assumed to be [0,1]. All series are assumed to be observed on a common grid of time points.
  ## all cell series in `Y` and `yy_smth` should be centered.
  ## TODO: verify whether it works when `interaction_term` is TRUE.


  # if only one kernel parameter set is specified, it is used for all variables
  if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y), kernel_params[[1]], simplify = FALSE)}

  # just in case `Y` and `yy_smth` are not centered at 0
  Y <- apply(Y, c(2,3), scale, center = TRUE, scale = FALSE)
  # DO NOT CENTER yy_smth!

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


  # Construct Sigma's for each cell
  if (verbose > 0) {cat("-------- constructing Sigma_k_kl for each cell --------\n")}

  ## with parallel computing
  Sigma_k_kl_multi <- parallel::mclapply(1:K,
                                   FUN = function(k) {
                                     .construct_Sigma_k_kl(interaction_term = interaction_term,
                                                           kernel = kernel,
                                                           kernel_params = kernel_params,
                                                           obs_time = obs_time,
                                                           tt = tt,
                                                           yy_smth = yy_smth[,,k])
                                   },
                                   mc.cores = mc.cores)
  Sigma_k_kl_multi <- simplify2array(Sigma_k_kl_multi)  # (n, n, p, K) or (n, n, p^2, K)

  if ((!interaction_term) & any(dim(Sigma_k_kl_multi) != c(n, n, p, K))) stop("Sigma_k_kl_multi has incorrect dimensions: `interaction_term` is FALSE.")
  if (interaction_term & any(dim(Sigma_k_kl_multi) != c(n, n, p^2, K))) stop("Sigma_k_kl_multi has incorrect dimensions: `interaction_term` is TRUE.")
  if (any(is.na(Sigma_k_kl_multi))) stop("NA occurs in Sigma_k_kl_multi.")


  # Iterative optimization
  res_bj <- rep(NA, p)  # length p
  res_cj <- matrix(NA, nrow = n, ncol = p)  # (n,p)
  res_best_eta <- rep(NA, p)  # length p
  res_best_kappa <- rep(NA, p)  # in theta_j estimation step

  # initialization: if no specified theta_initial, then initialize to 1 for all entries.
  if (is.null(theta_initial)){
    num_row <- ifelse(interaction_term, yes = p^2, no = p)
    theta_initial <- matrix(1, nrow = num_row, ncol = p)  # p^2 x p or p x p
  }
  res_theta <- theta_initial  # (p^2, p) or (p, p)

  # iteration begins
  for (num_iter in 1:max_iter){
    if (verbose > 0) {cat("\n-------- iteration", num_iter, " --------\n")}

    res_theta_prev <- res_theta

    # given theta_j, estimate F_j
    if (verbose > 0) {cat("-------- estimating Fj's for each cell --------\n")}

    for (j in 1:p) {
      if (verbose > 0) {cat("Working on variable", j, "\n")}

      yj_multi <- Y[,j,]
      theta_j <- res_theta[,j]

      # update Sigma using new theta_j
      Sigma_multi <- parallel::mclapply(1:K,
                                        FUN = function(k) {
                                          .construct_Sigma(Sigma_k_kl = Sigma_k_kl_multi[,,,k],
                                                           theta_j = theta_j)
                                        },
                                        mc.cores = mc.cores)
      Sigma_multi <- simplify2array(Sigma_multi)  # (n, n, K)

      res_fun_est <- function_estimation_multi(B = B,
                                               yj_multi = yj_multi,
                                               Sigma_multi = Sigma_multi,
                                               verbose = verbose)

      res_best_eta[j] <- res_fun_est$best_eta
      res_bj[j] <- res_fun_est$bj
      res_cj[,j] <- res_fun_est$cj

      if (verbose > 1) {cat("best eta =", res_fun_est$best_eta, "\n")}
    }

    # given F_j, estimate theta_j
    if (verbose > 0) {cat("-------- estimating theta_j's --------\n")}
    for (j in 1:p){
      G_multi <- parallel::mclapply(1:K,
                                    FUN = function(k) {
                                      .construct_G(cj = res_cj[,j],
                                                   interaction_term = interaction_term,
                                                   Sigma_k_kl = Sigma_k_kl_multi[,,,k])
                                    },
                                    mc.cores = mc.cores)
      G_multi <- do.call(rbind, G_multi)  # rbind G's for all cells, resulting in a matrix of dimension (Kn, p) or (Kn, p^2)

      if (any(dim(G_multi) != c(K*n, ifelse(interaction_term, p^2, p)))) {stop("G_multi should be of dimension (Kn, p) or (Kn, p^2).")}

      # lasso fit (optionally non-negative lasso)
      adj_col <- adj_matrix[,j]  # variable j corresponds to column j of the adjacency matrix. Note that `adj_col` is NULL if `adj_matrix` is NULL.

      # bj, cj and eta_j are shared across all cells
      res_bj_multi <- matrix(res_bj, nrow = p, ncol = K)  # (p,K)
      res_best_eta_multi <- matrix(res_best_eta, nrow = p, ncol = K)  # (p,K)
      res_cj_multi <- array(c(res_cj), dim = c(nrow(res_cj), ncol(res_cj), K))  # (n,p,K)

      res_theta_j_est <- theta_j_estimation_multi(bj_multi = res_bj_multi[j,],
                                                  B = B,
                                                  cj_multi = res_cj_multi[,j,],
                                                  eta_j_multi = res_best_eta_multi[j,],
                                                  G_multi = G_multi,
                                                  interaction_term = interaction_term,
                                                  yj_multi = Y[,j,],
                                                  adj_col = adj_col,
                                                  nzero_thres = nzero_thres)
      res_theta[,j] <- res_theta_j_est$theta_j
      res_best_kappa[j] <- res_theta_j_est$best_kappa
    }

    if (any(dim(theta_j) != c(ifelse(interaction_term, p^2, p), p))) {stop("Incorrect dimension of `theta_j`.")}

    # At each iteration, truncate small components of theta_j to 0
    res_theta[res_theta < 0.01] <- 0

    # check improvement in this iteration
    improvement <- norm(res_theta - res_theta_prev, type = "F") / norm(res_theta_prev, type = "F")

    if (verbose > 0) {
      cat("theta improvement:", round(improvement, 4), "\n")
      cat("F-norm of current theta:", round(norm(res_theta, type = "F"), 4), "\n")
    }

    # print the inferred network
    # only supported when `interaction_term` is FALSE.
    network_est <- theta_to_adj_matrix(interaction_term = interaction_term,
                                            res_theta = res_theta)
    if (verbose > 0) {
      cat("--- inferred network ---\n")
      print(network_est)
      cat("\n")
    }

    if (is.na(improvement)){  # usually because `norm(res_theta_prev, type = "F")` is 0 (i.e. `res_theta_prev` is all 0) at the first few iterations
      # warning("improvement is NA, iteration terminated")
      # improvement <- 0  # break the loop by the if statement below
      warning("improvement is NA")
      improvement <- 2*tol  # continue
    }

    # check stopping criteria
    if (improvement < tol) {break}
  }  # end of iteration loop

  ### CAVEAT ###
  ### DO NOT TAKE THE MEAN OF theta_j, bj, cj across all cells. The underlying math is more complicated than a simple average and may involve non-linearity.
  ### YOU WILL NOT OBTAIN A HIGHER PRECISION ESTIMATION BY SIMPLY TAKING MEAN.
  ### INSTEAD, USE `evaluate_Fj` to obtain the estimated derivative functions Fj for each cell, then take the mean of these Fj's (over all cells) for each variable j.

  if (verbose > 0) {cat("\n -------- finished -------- \n")}

  res_kernelODE <- list(network_est = network_est,
                        res_theta = res_theta,  # (p^2, p) or (p, p)
                        res_best_kappa = res_best_kappa,
                        res_bj = res_bj,  # length p
                        res_cj = res_cj,  # (n,p)
                        res_best_eta = res_best_eta,  # length p
                        num_iter = num_iter,
                        last_improvement = improvement,
                        config = list(
                          adj_matrix = adj_matrix,
                          interaction_term = interaction_term,
                          kernel = kernel,
                          kernel_params = kernel_params,
                          max_iter = max_iter,
                          nzero_thres = nzero_thres,
                          obs_time = obs_time,
                          theta_initial = theta_initial,
                          tol = tol,
                          tt = tt,
                          Y = Y,
                          yy_smth = yy_smth))

  return(res_kernelODE)
}
