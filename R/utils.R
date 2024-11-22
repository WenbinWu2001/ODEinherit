auto_select_kernel_params <- function(kernel, Y){
  # auto-selecting kernel parameters for each variable
  # `Y` is a matrix of shape (n, p).
  # returns `kernel_params` as a list of p sublists, each sublist containing the kernel params for a variable.
  # if `kernel_params` only contains one list, then that kernel parameter is used for all variables.
  if (kernel == "gaussian"){
    kernel_params <- compute_bw(Y)
  } else if (kernel == "linear"){
    kernel_params <- list(list())
  } else if (kernel == "polynomial") {
    kernel_params <- list(list(intercept=1, degree=3))
  } else if (kernel == "matern") {
    kernel_params <- list(list(lengthscale=1))
  } else {
    kernel_params <- NULL
    stop("Invalid kernel.")
  }

  return (kernel_params)
}


.construct_tt_mean <- function(obs_time, tt){
  ## computes \bar{T}(t) evaluated at each time point t in `tt`.
  ## `obs_time` are ti's on page (10) in the paper.
  ## tt_mean is a numeric vector of length len
  n <- length(obs_time)
  len <- length(tt)

  tt_mean <- rep(0, len)
  for (s in 1:n){tt_mean <- tt_mean + (tt <= obs_time[s])}
  tt_mean <- tt_mean / n

  return (tt_mean)
}


.construct_B <- function(obs_time){
  n <- length(obs_time)
  B <- rep(0, n)
  obs_time_mean <- mean(obs_time)
  for (i in 1:n){B[i] <- obs_time[i] - obs_time_mean}
  return (B)
}


.construct_V <- function(obs_time,
                         tt){
  ## V is a matrix with shape (len, n), where len = length(tt), n = length(obs_time)
  ## each column i of V corresponds to T_i(t) - \bar{T}(t) evaluated at each time point t in `tt`. i is the index for `obs_time` (i.e. ti's in the paper)
  tt_mean <- .construct_tt_mean(obs_time, tt)
  V <- sapply(1:length(obs_time), function(i) {(tt <= obs_time[i]) - tt_mean})

  return (V)
}


.construct_kk_main <- function(y1,
                               y2,
                               kernel = c("gaussian", "linear", "polynomial", "matern"),
                               param = NULL){
  ## Computes kernel values for main effects of ONE variable.
  ## This function is a vectorized version of `kernel_main`. Mainly used for efficient construction of Sigma's.
  ## `y1` and `y2` are two **vector** samples of the same variable of length n1 and n2 respectively.
  ## These kernel values are computed on the grid y1 x y2, thus returning a **matrix** of dimension (n1, n2). In other words, y1 is measured at n1 evenly spaced time points and y2 is measured at n2 evenly spaced time points.
  ## `param` is a single list containing the parameters for the specified kernel for this particular variable.
  ## For kernel="gaussian", it should include `bandwidth` as the bandwidth.
  ## For kernel="linear", no parameter is needed,
  ## For kernel="polynomial", it should include `intercept` and `degree`.
  ## For kernel="matern", it should include `lengthscale`.
  kernel <- match.arg(kernel)  # if no kernel is specified, use gaussian kernel

  grid <- expand.grid(y1 = y1, y2 = y2)  # expand the grid y1 x y2 for vectorized computation. Each row is a particular (y1, y2) element combination.
  if (kernel == "gaussian"){
    if (is.null(param$bandwidth)){stop("Please specify bandwidth for gaussian kernel.")}
    kk_main <- gaussian_kernel(y1 = grid$y1,
                               y2 = grid$y2,
                               bandwidth = param$bandwidth)
  } else if (kernel == "linear"){
    kk_main <- linear_kernel(y1 = grid$y1,
                             y2 = grid$y2)
  } else if (kernel == "polynomial"){
    if (is.null(param$intercept) | is.null(param$degree)){stop("Please specify intercept and degree for polynomial kernel.")}
    kk_main <- polynomial_kernel(y1 = grid$y1,
                                 y2 = grid$y2,
                                 intercept = param$intercept,
                                 degree = param$degree)
  } else if (kernel == "matern"){
    if (is.null(param$lengthscale)){stop("Please specify lengthscale for matern kernel.")}
    kk_main <- matern_kernel(y1 = grid$y1,
                             y2 = grid$y2,
                             lengthscale = param$lengthscale)
  } else {
    stop("Does not support the specified kernel.")
    kk_main <- NA
  }

  n1 <- length(y1)
  n2 <- length(y2)
  kk_main <- matrix(kk_main, nrow = n1, ncol = n2, byrow = FALSE)  # reshape to matrix

  if (any(is.na(kk_main) | is.nan(kk_main))) {stop("NaN produced in kk_main. Please check kernel computation.")}

  return (kk_main)
}


.construct_kk_inter <- function(Y1,
                                Y2,
                                k,
                                l,
                                kernel = c("gaussian", "linear", "polynomial", "matern"),
                                kernel_params = NULL){
  ## Computes kernel values for interaction effect between variables k and l (k is not equal to l).
  ## This function is a vectorized version of `kernel_inter`. Mainly used for efficient construction of Sigma's.
  ## `Y1` and `Y2` are two **matrix** samples of dimension (n1, p) and (n2, p) respectively. In other words, y1 is measured at n1 evenly spaced time points and y2 is measured at n2 evenly spaced time points.
  ## These kernel values are computed on the grid specified by rows of `Y1` and `Y2`, thus returning a **matrix** of dimension (n1, n2).
  ## `kernel_params` is a list of p lists. kernel_params[[k]] contains the parameters of variable k for the specified kernel. Check `kernel_main` to see the parameters needed for each kernel.
  if (ncol(Y1) != ncol(Y2)) {stop("Y1 and Y2 should contain the same variables.")}
  if (k == l) {stop("Interaction terms are only computed between different variables.")}
  kernel <- match.arg(kernel)  # if no kernel is specified, use gaussian kernel

  res_main_k <- .construct_kk_main(y1 = Y1[,k],
                                   y2 = Y2[,k],
                                   kernel = kernel,
                                   param = kernel_params[[k]])  # n1 x n2
  res_main_l <- .construct_kk_main(y1 = Y1[,l],
                                   y2 = Y2[,l],
                                   kernel = kernel,
                                   param = kernel_params[[l]])  # n1 x n2

  # K_{kl} = K_k \times K_l for interaction between variable k and l
  kk_inter <- res_main_k * res_main_l

  if (any(is.na(kk_inter) | is.nan(kk_inter))) {stop("NaN produced in kk_inter. Please check kernel computation.")}

  return(kk_inter)
}


.construct_kk_theta <- function(theta_j,
                                Y1,
                                Y2,
                                interaction_term,
                                kernel = c("gaussian", "linear", "polynomial", "matern"),
                                kernel_params = NULL){
  ## Computes values of kernel theta in the paper that generates the RKHS of F_j (which we estimate), specified by paragraph below Eq. (12).
  ## This function is a vectorized version of `kernel_theta`. Mainly used for efficient construction of Sigma's.
  ## `theta_j` is a vector with length p^2 if `interaction_term` == TRUE OR with length p if `interaction_term` == FALSE.
  ## When `interaction_term` == TRUE, the order of `theta_j` components (k and kl) differs from that in the paper, with the theta's of main effects on the diagonal and those of interaction effects off-diagonal, if `theta_j` is reshaped to a pxp square matrix.
  ## `Y1` and `Y2` are two **matrix** samples of dimension (n1, p) and (n2, p) respectively. In other words, `Y1` is measured at n1 evenly spaced time points and `Y2` is measured at n2 evenly spaced time points.
  ## These kernel values are computed on the grid specified by rows of `Y1` and `Y2`, thus returning a **matrix** of dimension (n1, n2).
  ## `kernel_params` is a list of p lists. kernel_params[[k]] contains the parameters of variable k for the specified kernel. Check `kernel_main` to see the parameters needed for each kernel.
  if (ncol(Y1) != ncol(Y2)) {stop("Y1 and Y2 should contain the same variables.")}
  if (interaction_term & (length(theta_j) != ncol(Y1)^2)) {stop("interaction_term = TRUE: theta_j should have length p^2.")}
  if ((!interaction_term) & (length(theta_j) != ncol(Y1))) {stop("interaction_term = FALSE: theta_j should have length p.")}
  if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y1), kernel_params[[1]], simplify = FALSE)}  # if only one parameter set is specified, it is used for all variables
  if (length(kernel_params) != ncol(Y1)) {stop("kernel_params should be of length p.")}
  kernel <- match.arg(kernel)  # if no kernel is specified, use gaussian kernel

  n1 <- nrow(Y1)
  n2 <- nrow(Y2)
  p <- ncol(Y1)

  # compute kernel values of main effects of all variables
  kk_main <- mapply(FUN = .construct_kk_main,
                    as.data.frame(Y1), as.data.frame(Y2), kernel, kernel_params,
                    SIMPLIFY = "array")  # apply to each column (variable) of Y1 and Y2 with element-wise argument specification. Simplified to 3D array.
  if (any(dim(kk_main) != c(n1, n2, p))) stop("Error in computing main effects in kk_main. Dimension doesn't match.")


  if (!interaction_term){
    kk_theta <- apply(sweep(kk_main,
                            MARGIN = 3,
                            STATS = theta_j,
                            FUN = "*"),  # weight each main effect (n1 x n2) matrix by theta_j
                      MARGIN = c(1,2),
                      sum)
  } else {  # compute kernel for interaction effect
    expanded_array <- array(rep(kk_main, p), dim = c(n1, n2, p, p))

    # specify diagonal indices, i.e. [i, j, k, k] for each i, j, k
    diag_indices <- expand.grid(1:n1, 1:n2, 1:p)
    diag_indices <- cbind(diag_indices, diag_indices[,3])
    diag_indices <- as.matrix(diag_indices)

    multiplied_array <- expanded_array * aperm(expanded_array, c(1, 2, 4, 3))  # On each [i,j,,] slice (pxp matrix), non-diagonal entries are interaction effects while diagonal entries are squares of main effects.
    multiplied_array[diag_indices] <- sqrt(multiplied_array[diag_indices])  # recover main effects on the diagonal by taking square root
    kk_array <- aperm(apply(multiplied_array, c(1, 2), c), c(2, 3, 1))  # reshape to 3D array with dimension (n1, n2, p^2) giving kernel values of both main effects (diagonal) and interaction effects (non-diagonal) on each [i,j,] slice
    kk_theta <- apply(sweep(kk_array,
                            MARGIN = 3,
                            STATS = theta_j,
                            FUN = "*"),   #  weight each main/interaction effect (n1 x n2) matrix by theta_j
                      MARGIN = c(1,2),
                      sum)
  }

  if (any(is.na(kk_theta) | is.nan(kk_theta))) {stop("NaN produced in kk_theta. Please check kernel computation.")}

  return (kk_theta)
}


.construct_kk_array <- function(Y1,
                                Y2,
                                interaction_term,
                                kernel = c("gaussian", "linear", "polynomial", "matern"),
                                kernel_params = NULL){
  ## Computes the kernel values of ALL main effects (and interaction effects if `interaction_term` == TRUE), returned as a 3D array.
  ## This function is vectorized. Mainly used for efficient construction of Sigma's.
  ## `Y1` and `Y2` are two **matrix** samples of dimension (n1, p) and (n2, p) respectively. In other words, `Y1` is measured at n1 evenly spaced time points and `Y2` is measured at n2 evenly spaced time points.
  ## These kernel values are computed on the grid specified by rows of `Y1` and `Y2`.
  ## Returns an **array** of dimension (n1, n2, p^2) if `interaction_term` == TRUE or (n1, n2, p) if `interaction_term` == FALSE. The order of components (k and kl, third dimension of the array) aligns with that of `theta_j`.
  ## The third dimension of the returned array is named to indicate the effect each (n1, n2) component matrix corresponds to.
  ## `kernel_params` is a list of p lists. kernel_params[[k]] contains the parameters of variable k for the specified kernel. Check `kernel_main` to see the parameters needed for each kernel.
  if (ncol(Y1) != ncol(Y2)) {stop("Y1 and Y2 should contain the same variables.")}
  if (length(kernel_params) == 1) {kernel_params <- replicate(ncol(Y1), kernel_params[[1]], simplify = FALSE)}  # if only one parameter set is specified, it is used for all variables
  if (length(kernel_params) != ncol(Y1)) {stop("kernel_params should be of length p.")}
  kernel <- match.arg(kernel)  # if no kernel is specified, use gaussian kernel

  n1 <- nrow(Y1)
  n2 <- nrow(Y2)
  p <- ncol(Y1)

  # compute kernel values of main effects of all variables
  kk_main <- mapply(FUN = .construct_kk_main,
                    as.data.frame(Y1), as.data.frame(Y2), kernel, kernel_params,
                    SIMPLIFY = "array")  # apply to each column (variable) of Y1 and Y2 with element-wise argument specification. Simplified to 3D array.
  if (any(dim(kk_main) != c(n1, n2, p))) stop("Error in computing main effects: kk_main should be n1 x n2 x p.")

  if (!interaction_term){
    # compute kernel gram matrices only for main effects
    kk_array <- kk_main
  } else {
    # compute kernel gram matrices for both main effects and interaction effects
    expanded_array <- array(rep(kk_main, p), dim = c(n1, n2, p, p))

    # specify diagonal indices, i.e. [i1, i1, j, j] for each (i1, i2, j) combination
    diag_indices <- expand.grid(1:n1, 1:n2, 1:p)
    diag_indices <- cbind(diag_indices, diag_indices[,3])
    diag_indices <- as.matrix(diag_indices)

    multiplied_array <- expanded_array * aperm(expanded_array, c(1, 2, 4, 3))  # On each [i,j,,] slice (pxp matrix), non-diagonal entries are interaction effects while diagonal entries are squares of main effects.
    multiplied_array[diag_indices] <- sqrt(multiplied_array[diag_indices])  # recover main effects on the diagonal by taking square root
    kk_array <- aperm(apply(multiplied_array, c(1, 2), c), c(2, 3, 1))  # reshape to 3D array with dimension (n1, n2, p^2) giving kernel values of both main effects (diagonal) and interaction effects (non-diagonal) on each [i,j,] slice
  }

  # name the 3rd dimension
  if (!interaction_term){
    dimnames(kk_array)[[3]] <- paste0("main_", 1:p)
  } else {
    names_3rdDim <- outer(1:p, 1:p, paste, sep = "_")
    names_3rdDim <- apply(names_3rdDim, c(1, 2), function(name) paste0("inter_", name))
    diag(names_3rdDim) <- paste0("main_", 1:p)
    dimnames(kk_array)[[3]] <- c(t(names_3rdDim))
  }

  if (any(is.na(kk_array) | is.nan(kk_array))) {stop("NaN produced in kk_array. Please check kernel computation.")}

  return (kk_array)
}


.kk_array_to_kk_theta <- function(kk_array,
                                  theta_j){
  ## Computes values of kernel theta as the sum of `kk_array` weighted by `theta_j`.
  ## `kk_array` is the output from `.construct_kk_array`, of dimension (n1, n2, p) for `interaction_term` == FALSE or (n1, n2, p^2) for `interaction_term` == TRUE.
  ## `theta_j` should have the same components order (k and kl) as the third dimension of `kk_array`.
  kk_theta <- apply(sweep(kk_array,
                          MARGIN = 3,
                          STATS = theta_j,
                          FUN = "*"),  # weight each main/interaction effect (n1 x n2) matrix by theta_j
                    MARGIN = c(1,2),
                    sum)

  if (any(is.na(kk_theta) | is.nan(kk_theta))) {stop("NaN produced in kk_theta. Please check kernel computation.")}

  return (kk_theta)
}


.construct_Sigma_k_kl <- function(interaction_term,
                                  kernel,
                                  kernel_params,
                                  obs_time,
                                  tt,
                                  yy_smth){
  ## returns an array of dimension (n, n, p^2) if `interaction_term` == TRUE or (n, n, p) if `interaction_term` == FALSE, where n = length(obs_time).
  ## return value is ALWAYS on original (non-log) scale, no matter `kk_array` is computed on log scale or not.

  p <- ncol(yy_smth)
  len <- nrow(yy_smth)
  n <- length(obs_time)

  delta <- 1 / nrow(yy_smth)

  V <- .construct_V(obs_time = obs_time, tt = tt)
  kk_array <- .construct_kk_array(Y1 = yy_smth,
                                  Y2 = yy_smth,
                                  interaction_term = interaction_term,
                                  kernel = kernel,
                                  kernel_params = kernel_params)

  Sigma_k_kl <- sapply(1:dim(kk_array)[3],
                             function(k) {crossprod(V, kk_array[,,k]) %*% V * delta^2},
                             simplify = "array")  # V^T kk_array[,,k] V
  dimnames(Sigma_k_kl)[[3]] <- dimnames(kk_array)[[3]]

  return (Sigma_k_kl)
}


.construct_Sigma <- function(Sigma_k_kl,
                            theta_j){
  ## construct Sigma given in top of page 11 in the paper.
  ## returns a matrix of dimension (n, n).
  ## `Sigma_k_kl` should always be on original (non-log) scale.
  ## returned value is always on original (non-log) scale.
  if (length(c(theta_j)) != dim(Sigma_k_kl)[3]) {stop("length(theta_j) should be the same as dim(Sigma_k_kl)[3].")}

  theta_j <- c(theta_j)  # convert to vector form

  Sigma <- apply(sweep(Sigma_k_kl,
                       MARGIN=3,
                       STATS=theta_j,
                       FUN="*"),  #  weight each Sigma_k or Sigma_kl matrix by theta_j
                 MARGIN = c(1,2),
                 sum)
  return (Sigma)
}


.construct_G <- function(cj,
                         interaction_term,
                         Sigma_k_kl){
  ## columns of G has different order as in the paper. The order aligns with that of `theta_j` (main effects on the diagonal and interaction effects off-diagonal, if `theta_j` is reshaped to a pxp square matrix.)
  if (dim(Sigma_k_kl)[1] != length(c(cj))) {stop("Sigma_k_kl[,,j] should be n x n. Here n can be obtained by length(cj).")}

  n <- dim(Sigma_k_kl)[1]
  cj <-  matrix(cj, nrow = n, ncol = 1)
  G <- sapply(1:dim(Sigma_k_kl)[3],  # 1:p without interaction or 1:p^2 with interaction
              function(kl) {Sigma_k_kl[,,kl] %*% cj},
              simplify = TRUE)

  if (any(is.na(G) | is.nan(G))) {stop("NA occurs in G.")}

  return (G)
}


.map_to <- function(vec_from, vec_to){
  # map each element in `vec_from` to the closest on in `vec_to` and return its index in `vec_to`.
  # alternative for `match` avoiding floating point problems
  n <- length(vec_from)
  res <- sapply(1:n, function(i){
    which.min(abs(vec_to - vec_from[i]))
  })
  res
}


theta_to_adj_matrix <- function(interaction_term, res_theta){
  # `res_theta` is (p^2, p) if `interaction_term` is TRUE or (p, p) if `interaction_term` is FALSE.
  p <- ncol(res_theta)

  if (interaction_term & (nrow(res_theta) != p^2)) {stop("res_theta should be (p^2, p) when `interaction_term` is TRUE.")}
  if ((!interaction_term) & (nrow(res_theta) != p)) {stop("res_theta should be (p, p) when `interaction_term` is FALSE.")}

  if (!interaction_term) {  # without interaction
    adj_matrix <- ifelse(res_theta > 0, yes = 1, no = 0)
  } else {  # with interaction
    # TODO: how to generate a network for model with interaction
    warning("network construction with interaction_term not supported. Returned a fully connected network.")
    adj_matrix <- matrix(1, nrow = p, ncol = p)
  }

  return (adj_matrix)
}

