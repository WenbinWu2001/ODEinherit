.function_estimation <- function(Yj,
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

  best_eta <- eta_vector[which.min(GCVs)]

  W <- Sigma
  diag(W) <- diag(Sigma) + n*best_eta

  cj <- Q2 %*% tcrossprod(solve(crossprod(Q2, W) %*% Q2), Q2) %*% (Yj - mean(Yj))
  bj <- crossprod(Q1, Yj - mean(Yj) - W%*%cj) / R

  return (list(bj = bj, cj = cj, best_eta = best_eta))
}
