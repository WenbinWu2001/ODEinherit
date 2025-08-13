.compute_loss <- function(bj,
                          cj,
                          eta_j,
                          interaction_term,
                          kappa_j,
                          kernel,
                          kernel_params,
                          kk_array = NULL,
                          obs_time,
                          Sigma_k_kl,
                          theta_j,
                          tt,
                          Yj,
                          yy_smth){
  # compute the objective function value for a specific variable j
  # Eq (11) in the kernel ODE paper

  obs_idx <- .map_to(obs_time, tt)

  res_recov <- evaluate_Fj(bj = bj,
                           cj = cj,
                           interaction_term = interaction_term,
                           kernel = kernel,
                           kernel_params = kernel_params,
                           kk_array = kk_array,
                           obs_time = obs_time,
                           theta_j = theta_j,
                           tt = tt,
                           Yj = Yj,
                           yy_smth = yy_smth)
  Yj_est <- res_recov$yy_est[obs_idx]

  # MSE
  MSE <- mean((Yj - Yj_est)^2)

  # penalty term
  PkFj_norm_sq <- sapply(1:length(theta_j), function(k_kl){t(cj) %*% Sigma_k_kl[,,k_kl] %*% cj * theta_j[k_kl]^2})
  term1 <- sum(PkFj_norm_sq / theta_j, na.rm = TRUE) * eta_j  # theta_j could be 0
  term2 <- sum(theta_j) * kappa_j
  penalty <- term1 + term2

  # overall loss
  loss <- MSE + penalty

  return (list(loss = loss,
               MSE = MSE,
               penalty = penalty))
}
