#' Title
#'
#' This is an iterative process.
#' @param Y
#' @param obs_time
#' @param yy_smth
#' @param tt
#' @param kernel
#' @param kernel_params
#' @param interaction_term
#' @param theta_initial
#' @param adj_matrix
#' @param nzero_thres
#' @param eval_edge_R2
#' @param eval_loss
#' @param tol
#' @param max_iter
#' @param verbose
#'
#' @returns
#' @export
#'
#' @examples
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
                            verbose = 0){

  ## kernel ODE Step 2: Iterative optimization algorithm (single sample version)
  # structure check of kernel_params

  if (!(is.list(kernel_params) & all(sapply(kernel_params, is.list)))) {stop("kernel_params should be of a list of lists.")}  # must be a list of lists
  if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y), kernel_params[[1]], simplify = FALSE)}  # if only one parameter set is specified, it is used for all variables
  if (length(kernel_params) != ncol(Y)) {stop("kernel_params should be of length p.")}

  n <- nrow(Y)
  p <- ncol(Y)

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

  # initialization: if no specified theta_initial, then initialize to 1 for all entries.
  res_bj <- rep(NA, p)
  res_cj <- matrix(NA, nrow = n, ncol = p)
  res_best_eta <- rep(NA, p)
  res_best_kappa <- rep(NA, p)  # in theta_j estimation step
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
    res_Fj_est_list <- lapply(1:p, FUN = function(j){
      Yj <- Y[,j]
      theta_j <- res_theta[,j]
      Sigma <- .construct_Sigma(Sigma_k_kl = Sigma_k_kl,
                                theta_j = theta_j)  # update Sigma using new theta_j
      res_fun_est <- .function_estimation(Yj = Yj,
                                          Sigma = Sigma,
                                          Q1 = Q1,
                                          Q2 = Q2,
                                          R = R,
                                          verbose = verbose)
      res_fun_est
    })

    res_best_eta <- sapply(res_Fj_est_list, function(lst){lst$best_eta})
    res_bj <- sapply(res_Fj_est_list, function(lst){lst$bj})
    res_cj <- sapply(res_Fj_est_list, function(lst){lst$cj})


    # given F_j, estimate theta_j
    if (verbose > 0) {cat("-------- estimating theta_j's --------\n")}
    res_theta_j_est_list <- lapply(1:p, FUN = function(j){
      Yj <- Y[,j]
      bj <- res_bj[j]
      cj <- res_cj[,j]
      eta_j <- res_best_eta[j]
      G <- .construct_G(cj = cj,
                        interaction_term = interaction_term,
                        Sigma_k_kl = Sigma_k_kl)

      # non-negative lasso fit
      adj_col <- adj_matrix[,j]  # variable j corresponds to column j of the adjacency matrix. Note that `adj_col` is NULL if `adj_matrix` is NULL.
      res_theta_j_est <- .theta_j_estimation(bj = bj,
                                             B = B,
                                             cj = cj,
                                             eta_j = eta_j,
                                             G = G,
                                             interaction_term = interaction_term,
                                             Yj = Yj,
                                             adj_col = adj_col,
                                             nzero_thres = nzero_thres)
      res_theta_j_est
    })

    res_theta <- sapply(res_theta_j_est_list, function(lst){lst$theta_j})
    res_best_kappa <- sapply(res_theta_j_est_list, function(lst){lst$best_kappa})


    # (optional) evaluate the loss function
    if (eval_loss) {
      res_loss_iter <- lapply(1:p, FUN = function(j){
        .compute_loss(bj = res_bj[j],
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
                      yy_smth = yy_smth)
        })

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
      # This usually happens for one of two reasons:
      # (1) `norm(res_theta_prev, type = "F")` is 0 (i.e. `res_theta_prev` is all 0) at the first few iterations, or
      # (2) the jth column of the input network `adj_matrix` is all zero (i.e. no variable regulates variable j) so that theta_j is all zero.

      # warning("improvement is NA")  # suppressed to make output cleaner
      improvement <- 2*tol  # ignore and continue
    }

    if (improvement < tol) {break}  # check stopping criteria
  }  # All iteration finish.

  # after the iterative process, truncate small components of theta_j to 0
  res_theta[res_theta < 0.01] <- 0

  if (verbose > 0) {cat("\n -------- finished -------- \n")}

  # obtain the regulatory network
  if (!interaction_term) {  # without interaction
    network_est <- ifelse(res_theta > 0, yes = 1, no = 0)
  } else { # with interaction
    # network for model with interaction, currently NOT supported
    warning("network construction with interaction_term not supported. Returned a fully connected network.")
    network_est <- matrix(1, nrow = p, ncol = p)
  }

  res_step2 <- list(res_theta = res_theta,
                    res_best_kappa = res_best_kappa,
                    res_bj = res_bj,
                    res_cj = res_cj,
                    res_best_eta = res_best_eta,
                    res_loss_path = res_loss_path,
                    network_est = network_est,
                    num_iter = num_iter,
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
