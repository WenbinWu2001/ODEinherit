.stability_selection <- function(B,
                                 interaction_term,
                                 mc.cores = 1,
                                 nzero_thres = NULL,
                                 Q1,
                                 Q2,
                                 res_theta_multi,
                                 R,
                                 stability_thres = 0.6,
                                 Sigma_k_kl,
                                 verbose = 0,
                                 Y){
  ## `res_theta_multi` is of shape (p^2, p, K) or (p, p, K). 'Y' is of shape (n, p, K).
  # stability selection by estimating the network `adj_matrix_multi` by all cells (supported only when 'interaction_term' is FALSE)
  if (verbose > 0) {cat("-------- stability selection --------\n")}

  n <- dim(Y)[1]
  p <- dim(Y)[2]
  K <- dim(Y)[3]

  res_multi_temp <- parallel::mclapply(1:K,
                                       FUN = function(k) {
                                         # infer the network using each cell
                                         # if (verbose > 1) {cat("\nWorking on cell ", k, "\n")}  # I suggest only output this message when no parallel computing is used

                                         res_theta_single <- res_theta_multi[,,k]  # initialize to store the estimated theta_j's for this cell
                                         res_best_kappa_single <- rep(NA, p)
                                         # res_bj_single <- rep(NA, p)
                                         # res_cj_single <- matrix(NA, nrow = n, ncol = p)
                                         # res_best_eta_single <- rep(NA, p)

                                         for (j in 1:p) {
                                           # function estimation: given theta_j, estimate F_j
                                           if (verbose > 1) {cat("Working on variable ", j, "\n")}

                                           Yj <- Y[,j,k]
                                           theta_j <- res_theta_single[,j]  # current theta_j
                                           Sigma <- .construct_Sigma(Sigma_k_kl = Sigma_k_kl[,,,k],
                                                                     theta_j = theta_j)  # update Sigma using new theta_j
                                           res_fun_est <- function_estimation(Yj = Yj,
                                                                              Sigma = Sigma,
                                                                              Q1 = Q1,
                                                                              Q2 = Q2,
                                                                              R = R,
                                                                              verbose = verbose)
                                           # res_best_eta_single[j] <- res_fun_est$best_eta
                                           # res_bj_single[j] <- res_fun_est$bj
                                           # res_cj_single[,j] <- res_fun_est$cj

                                           # theta estimation: given F_j, estimate theta_j
                                           bj <- res_fun_est$bj
                                           cj <- res_fun_est$cj
                                           eta_j <- res_fun_est$best_eta
                                           G <- .construct_G(cj = cj,
                                                             interaction_term = interaction_term,
                                                             Sigma_k_kl = Sigma_k_kl[,,,k])

                                           # non-negative lasso
                                           res_theta_j_est <- theta_j_estimation(bj = bj,
                                                                                 B = B,
                                                                                 cj = cj,
                                                                                 eta_j = eta_j,
                                                                                 G = G,
                                                                                 interaction_term = interaction_term,
                                                                                 Yj = Yj,
                                                                                 adj_col = NULL,
                                                                                 nzero_thres = nzero_thres)
                                           res_theta_single[,j] <- res_theta_j_est$theta_j  # update theta_j
                                           res_best_kappa_single[j] <- res_theta_j_est$best_kappa
                                         }

                                         # obtain inferred network by cell k
                                         res_theta_single[res_theta_single < 0.01] <- 0  # coerce small values to 0
                                         adj_matrix_single <- theta_to_adj_matrix(interaction_term = interaction_term,
                                                                                  res_theta = res_theta_single)

                                         res_single <- list(adj_matrix_single = adj_matrix_single)
                                         # res_single <- list(adj_matrix_single = adj_matrix_single,
                                         #                    res_bj_single = res_bj_single,
                                         #                    res_cj_single = res_cj_single,
                                         #                    res_best_eta_single = res_best_eta_single)
                                         return (res_single)
                                       },
                                       mc.cores = mc.cores)

  # convert to array
  adj_matrix_multi <- simplify2array(lapply(res_multi_temp, function(res_single){res_single$adj_matrix_single}))  # shape is (p^2, p, K) or (p, p, K). [,,k] slice is the inferred network by cell k
  # res_bj_multi <- simplify2array(lapply(res_multi_temp, function(res_single){res_single$res_bj_single}))  # (p, K)
  # res_cj_multi <- simplify2array(lapply(res_multi_temp, function(res_single){res_single$res_cj_single}))  # (n, p, K)
  # res_best_eta_multi <- simplify2array(lapply(res_multi_temp, function(res_single){res_single$res_best_eta_single}))  # (p, K)

  # stability selection step: form the network (adjacency matrix) where an edge is selected only if it is selected more than `stability_thres` among all cells
  sel_rate <- apply(adj_matrix_multi, c(1,2), mean)
  adj_matrix <- ifelse(sel_rate >= stability_thres, yes = 1, no = 0)  # of shape (p, p) after pooling the networks from all cells

  if (any(dim(adj_matrix) != c(p, p))) {stop("`adj_matrix` should be pxp.")}
  if (!all(adj_matrix %in% c(0, 1))) {stop("`adj_matrix` should be binary.")}
  # if (any(dim(res_cj_multi) != c(n, p, K))) {stop("`res_cj_multi` should be (n, p, K).")}
  # if (any(dim(res_bj_multi) != c(p, K))) {stop("`res_bj_multi` should be (p, K).")}
  # if (any(dim(res_best_eta_multi) != c(p, K))) {stop("`res_best_eta_multi` should be (p, K).")}

  # res_multi <- list(adj_matrix = adj_matrix,
  #                   adj_matrix_multi = adj_matrix_multi,
  #                   res_bj_multi = res_bj_multi,
  #                   res_cj_multi = res_cj_multi,
  #                   res_best_eta_multi = res_best_eta_multi)

  res_multi <- list(adj_matrix = adj_matrix,
                    adj_matrix_multi = adj_matrix_multi)

  return (res_multi)
}
