function_estimation <- function(Yj,
                                Sigma,
                                Q1,
                                Q2,
                                R,
                                eta_vector = 10^seq(-10, 0, 0.25),
                                verbose = 0){
  # `Yj` is a vector of length n. It is the time series of variable j.
  # `Sigma` is of dim (n, n).
  # `Q1`, `Q2` are nx1. `R` is a scalar.

  n <- nrow(Sigma)
  Yj <- matrix(Yj, nrow = n, ncol = 1)  # convert to matrix form for matrix operation

  # GCV tuning for eta
  num_eta <- length(eta_vector)
  GCVs <- rep(NA, num_eta)
  for (i in 1:num_eta){
    if(verbose > 1 && num_eta > 10 && i %% floor(num_eta/10) == 0) cat('*')  # update per 10 trials

    eta <- eta_vector[i]
    tryCatch({  # t(Q2) %*% W %*% Q2 may be singular for small eta
      W <- Sigma
      diag(W) <- diag(Sigma) + n*eta

      A <- -n*eta*Q2 %*% tcrossprod(solve(crossprod(Q2, W) %*% Q2), Q2)
      diag(A) <- diag(A) + 1

      A_temp <- A
      diag(A_temp) <- diag(A_temp) - 1  # A_temp = A - In
      GCVs[i] <- sum((A_temp %*% (Yj - mean(Yj)))^2) / (sum(diag(A_temp))/n)^2
    }, error = function(e){})  # if matrix inversion fails, simply use next (larger) eta
  }
  if(verbose > 1) cat("\n")

  best_eta <- eta_vector[which.min(GCVs)]

  W <- Sigma
  diag(W) <- diag(Sigma) + n*best_eta

  cj <- Q2 %*% tcrossprod(solve(crossprod(Q2, W) %*% Q2), Q2) %*% (Yj - mean(Yj))
  bj <- crossprod(Q1, Yj - mean(Yj) - W%*%cj) / R

  return (list(bj = bj, cj = cj, best_eta = best_eta))
}

# function_estimation_alt2 <- function(Yj,
#                                      Sigma,
#                                      Q1,
#                                      Q2,
#                                      R,
#                                      K,
#                                      eta_vector = 10^seq(-10, 0, 0.25),
#                                      verbose = 0){
#   # `Yj` is a vector of length n. It is the time series of variable j.
#   # `Sigma` is of dim (n, n).
#   # `Q1`, `Q2` are nx1. `R` is a scalar.
#   # `K` is the number of cells, though this function only applies to one cell series.
#
#   n <- nrow(Sigma)
#   Yj <- matrix(Yj, nrow = n, ncol = 1)  # convert to matrix form for matrix operation
#
#   # GCV tuning for eta
#   num_eta <- length(eta_vector)
#   GCVs <- rep(NA, num_eta)
#   for (i in 1:num_eta){
#     if(verbose > 1 && num_eta > 10 && i %% floor(num_eta/10) == 0) cat('*')  # update per 10 trials
#
#     eta <- eta_vector[i]
#     tryCatch({  # t(Q2) %*% W %*% Q2 may be singular for small eta
#       W <- Sigma
#       diag(W) <- diag(Sigma) + n*eta * K  ### multi-sample version ###
#
#       A <- -n*eta*Q2 %*% tcrossprod(solve(crossprod(Q2, W) %*% Q2), Q2)
#       diag(A) <- diag(A) + 1
#
#       A_temp <- A
#       diag(A_temp) <- diag(A_temp) - 1  # A_temp = A - In
#       GCVs[i] <- sum((A_temp %*% (Yj - mean(Yj)))^2) / (sum(diag(A_temp))/n)^2
#     }, error = function(e){})  # if matrix inversion fails, simply use next (larger) eta
#   }
#   if(verbose > 1) cat("\n")
#
#   best_eta <- eta_vector[which.min(GCVs)]
#
#   W <- Sigma
#   diag(W) <- diag(Sigma) + n*best_eta * K  ### multi-sample version ###
#
#   cj <- Q2 %*% tcrossprod(solve(crossprod(Q2, W) %*% Q2), Q2) %*% (Yj - mean(Yj))
#   bj <- crossprod(Q1, Yj - mean(Yj) - W%*%cj) / R
#
#   return (list(bj = bj, cj = cj, best_eta = best_eta))
# }


function_estimation_multi <- function(B,
                                      yj_multi,
                                      Sigma_multi,
                                      eta_vector = 10^seq(-10, 0, 0.25),
                                      verbose = 0){
  # `yj_multi` is a matrix of dim (n, K) containing the time series data of variable j for all experiments.
  # `Sigma_multi` is an array with dim (n, n, K) containing the nxn Sigma's for all experiments.

  n <- dim(Sigma_multi)[1]
  K <- dim(Sigma_multi)[3]  # number of experiments

  # shared quantities across tuning
  yj_tilde_multi <- apply(yj_multi, 2, scale, center = TRUE, scale = FALSE)  # (n,K), centered Yj (Yj should already be centered before being fed into kernel ODE, but we follow the math for clarity)
  sum_yj_tilde <- matrix(apply(yj_tilde_multi, 1, sum), ncol = 1)  # (n,1)
  Lambda <- apply(Sigma_multi, c(1,2), sum)  # (n,n)  # sum of Sigma's from all cells
  q <- matrix(B / sqrt(sum(B^2)), ncol = 1)  # (n,1)

  zj <- matrix(
    c(apply(yj_tilde_multi, 2, function(yj_tilde_r){yj_tilde_r - q %*% crossprod(q, sum_yj_tilde) / K})),
    ncol = 1)  # (Kn,1)
  G <- lapply(1:K, function(k){Sigma_multi[,,k] - q %*% crossprod(q, Lambda) / K})
  G <- do.call(rbind, G)  # (Kn,n)
  GTG <- crossprod(G)  # (n,n)

  # GCV tuning for eta
  num_eta <- length(eta_vector)
  GCVs <- rep(NA, num_eta)
  for (i in 1:num_eta){
    if(verbose > 1 && num_eta > 10 && i %% floor(num_eta/10) == 0) cat('*')  # update per 10 trials

    eta <- eta_vector[i]
    tryCatch({
      W <- GTG + K*n * eta * Lambda
      W_inv <- MASS::ginv(W)  # use pseudo inverse to obtain one valid solution, since Sigma (thus GTG and Lambda) is rank deficient and gives infinitely many solutions
      A <- G %*% tcrossprod(W_inv, G)

      A_temp <- A
      diag(A_temp) <- diag(A_temp) - 1  # A_temp = A - I_{Kn}
      GCVs[i] <- K*n * sum((A_temp %*% zj)^2) / (sum(diag(A_temp)))^2
    }, error = function(e){})  # if matrix inversion fails, simply use next (larger) eta
  }

  best_eta <- eta_vector[which.min(GCVs)]
  W <- GTG + K*n * best_eta * Lambda
  W_inv <- MASS::ginv(W)
  cj <- W_inv %*% crossprod(G, zj)
  bj <- crossprod(B, sum_yj_tilde - Lambda %*% cj) / (K * sum(B^2))

  return (list(bj = bj, cj = cj, best_eta = best_eta))
}



