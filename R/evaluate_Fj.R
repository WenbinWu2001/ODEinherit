evaluate_Fj <- function(bj,
                        cj,
                        interaction_term,
                        kernel,
                        kernel_params,
                        kk_array = NULL,
                        obs_time,
                        theta_j,
                        tt,
                        Yj,
                        yy_smth){
  ## evaluate Fj at the integration grid tt. Returns the estimated initial condition (a scalar), the estimated Fj on tt, the reconstructed trajectory on tt, and the total variation of the reconstructed trajectory.
  ## Use Eq (13) on page 10.
  ## `Yj`: a vector of length n. It is the observed time series for variable j.
  ## `kk_array`: The analog of Sigma_k_kl on the integration grid `tt`. It is the output from `.construct_kk_array`, of dimension (len, len, p) for `interaction_term` == FALSE or (len, len, p^2) for `interaction_term` == TRUE, where `len` is `length(tt)`.
  ## It is shared across all variables. Pre-computing it and passing it into this function for reusing can largely save computations when iterating over all variables and evaluating their Fj's.

  n <- length(obs_time)
  p <- ncol(yy_smth)
  len <- length(tt)
  delta <- 1/len
  tt_mean <- .construct_tt_mean(obs_time, tt)

  if ((length(bj) != 1)) {stop("bj should be a scalar.")}
  if (length(cj) != n) {stop("cj should be a vector of length n.")}
  if (length(theta_j) != ifelse(interaction_term, p^2, p)) {stop("theta_j should be a vector of length p or p^2.")}
  if ((!is.null(kk_array)) & any(dim(kk_array) != c(len, len, length(theta_j)))) {stop("kk_array should be an array of dimension (len, len, p) or (len, len, p^2).")}

  cj <- matrix(cj, ncol = 1)  # nx1

  # construct kk_theta
  if (is.null(kk_array)){  # without given kk_array
    kk_theta <- .construct_kk_theta(theta_j,
                                    Y1 = yy_smth,
                                    Y2 = yy_smth,
                                    interaction_term = interaction_term,
                                    kernel = kernel,
                                    kernel_params = kernel_params)
  } else {  # given kk_array
    kk_theta <- .kk_array_to_kk_theta(kk_array,
                                      theta_j)
  }

  # evaluate Fj
  V <- matrix(NA, nrow = len, ncol = n)
  for (i in 1:n){
    V[ ,i] <- matrix((tt <= obs_time[i]) - tt_mean, nrow = length(tt), ncol = 1)
  }
  Fj_est <- (kk_theta %*% V * delta) %*% cj + c(bj)  # note bj is a scalar
  Fj_est <- c(Fj_est)

  # estimate the initial condition theta_j0 (above Eq. 12)
  theta_j0 <- mean(Yj) - sum(Fj_est * tt_mean) * delta

  # recover the trajectory xj
  yy_est <- theta_j0 + cumsum(Fj_est) * delta

  # calculate the total variation of the trajectory (= integration of |Fj| over time)
  TV_est <- sum(abs(Fj_est)) * delta

  return (list(theta_j0 = theta_j0,  # scalar
               Fj_est = Fj_est,  # vector of length `len`
               yy_est = yy_est,  # vector of length `len`
               TV_est = TV_est))  # scalar
}
