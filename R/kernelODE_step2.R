kernelODE_step2 <- function(Y,
                            obs_time,
                            yy_smth,
                            tt,
                            kernel,
                            kernel_params,
                            interaction_term = FALSE,
                            theta_initial = NULL,
                            adj_matrix = NULL,
                            nzero_thres = NULL,
                            eval_edge_R2 = FALSE,
                            eval_loss = FALSE,
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

  # # ensure `Y` is centered
  # Y <- scale(Y, center = TRUE, scale = FALSE)
  # # DO NOT CENTER yy_smth!

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
                                      yy_smth = yy_smth)  # (n, n, p^2) or (n, n, p)
  kk_array <- .construct_kk_array(Y1 = yy_smth,
                                  Y2 = yy_smth,
                                  interaction_term = interaction_term,
                                  kernel = kernel,
                                  kernel_params = kernel_params)  # (len, len, p) or (len, len, p^2)

  # Iterative optimization
  res_bj <- rep(NA, p)
  res_cj <- matrix(NA, nrow = n, ncol = p)
  res_best_eta <- rep(NA, p)
  res_best_kappa <- rep(NA, p)  # in theta_j estimation step

  # initialization: if no specified theta_initial, then initialize to 1 for all entries.
  if (is.null(theta_initial)){
    num_row <- ifelse(interaction_term, yes = p^2, no = p)
    theta_initial <- matrix(1, nrow = num_row, ncol = p)  # p^2 x p or p x p
  }
  res_theta <- theta_initial
  res_loss_path <- list()

  # iteration begins
  for (num_iter in 1:max_iter){
    if (verbose > 0) {cat("\n-------- iteration", num_iter, " --------\n")}

    res_theta_prev <- res_theta

    # given theta_j, estimate F_j
    if (verbose > 0) {cat("-------- estimating Fj's --------\n")}

    for (j in 1:p) {
      if (verbose > 1) {cat("Working on variable", j, "\n")}

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

    if (eval_loss) {
      # compute the loss function
      res_loss_iter <- lapply(1:p, function(j){compute_loss(bj = res_bj[j],
                                                            cj = res_cj[,j],
                                                            eta_j = res_best_eta[j],
                                                            interaction_term = interaction_term,
                                                            kappa_j = res_best_kappa[j],
                                                            kernel = kernel,
                                                            kernel_params = kernel_params,
                                                            kk_array = kk_array,
                                                            obs_time = obs_time,
                                                            Sigma_k_kl = Sigma_k_kl,
                                                            theta_j = res_theta[,j],
                                                            tt = tt,
                                                            Yj = Y[,j],
                                                            yy_smth = yy_smth)})
      res_loss_path[[num_iter]] <- res_loss_iter
    }

    # check improvement in this iteration
    improvement <- norm(res_theta - res_theta_prev, type = "F") / norm(res_theta_prev, type = "F")

    if (verbose > 0) {
      cat("theta improvement:", round(improvement, 4), "\n")
      cat("F-norm of current theta:", round(norm(res_theta, type = "F"), 4), "\n")
      if (eval_loss) {
        cat("obj fun value of each variable:",
            paste(round(sapply(res_loss_iter, function(obj){obj$loss}), 2),
                  sep = ", "))
      }
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


  res_step2 <- list(res_theta = res_theta,
                    res_best_kappa = res_best_kappa,
                    res_bj = res_bj,
                    res_cj = res_cj,
                    res_best_eta = res_best_eta,
                    res_loss_path = res_loss_path,
                    num_iter = num_iter,
                    last_improvement = improvement,
                    config = list(adj_matrix = adj_matrix,
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

  return (res_step2)
}
