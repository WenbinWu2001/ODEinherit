linear_kernel <- function(y1, y2){
  ## K(y1, y2) = y1*y2, equivalent to polynomial kernel with intercept=0, degree=1.
  ## Note this function is vectorized.
  if (length(y1) != length(y2)) {stop("y1 and y2 should have same length.")}

  res <- y1*y2

  if (any(is.na(res) | is.nan(res))) {stop("Computing linear kernel: non-positive value occurs in log(). NaN produced.")}

  return (res)
}


gaussian_kernel <- function(y1, y2, bandwidth){
  ## K(y1, y2) = exp(- (y1-y2)^2 / (2*bandwidth^2))
  ## Note this function is vectorized.
  if (length(y1) != length(y2)) {stop("y1 and y2 should have same length.")}
  if (is.na(bandwidth) | is.null(bandwidth)) {stop("bandwidth is NA or NULL.")}
  if (bandwidth <= 0) {stop("bandwidth should be positive.")}

  res <- exp(- (y1-y2)^2 / (2 * bandwidth^2))

  if (any(is.na(res) | is.nan(res))) {stop("Computing gaussian kernel: non-positive value occurs in log(). NaN produced.")}

  return (res)
}


polynomial_kernel <- function(y1, y2, intercept, degree){
  ## K(y1, y2) = (y1*y2 + intercept)^degree
  ## Note this function is vectorized.
  if (length(y1) != length(y2)) {stop("y1 and y2 should have same length.")}
  if (is.na(intercept) | is.null(intercept)) {stop("intercept is NA or NULL.")}
  if (is.na(degree) | is.null(degree)) {stop("degree is NA or NULL.")}
  if ((degree < 0) | (degree%%1 != 0)) {stop("degree should be a non-negative integer.")}

  res <- (y1*y2 + intercept)^degree

  if (any(is.na(res) | is.nan(res))) {stop("Computing polynomial kernel: non-positive value occurs in log(). NaN produced.")}

  return (res)
}


matern_kernel <- function(y1, y2, lengthscale){
  ## Matern kernel function with nu = 3/2
  ## K(y1, y2) = (1 + \sqrt(3)*|y1-y2|/lengthscale) * exp(-\sqrt(3)*|y1-y2|/lengthscale)
  ## Note this function is vectorized.
  if (length(y1) != length(y2)) {stop("y1 and y2 should have same length.")}
  if (is.na(lengthscale) | is.null(lengthscale)) {stop("lengthscale is NA or NULL.")}
  if (lengthscale <= 0) {stop("lengthscale should be positive.")}

  res <- (1 + sqrt(3)*abs(y1-y2)/lengthscale) * exp(- sqrt(3)*abs(y1-y2)/lengthscale)

  if (any(is.na(res) | is.nan(res))) {stop("Computing matern kernel: non-positive value occurs in log(). NaN produced.")}

  return (res)
}


kernel_main <- function(y1, y2, kernel = c("gaussian", "linear", "polynomial", "matern"), param = NULL){
  ## Computes K_k(y1, y2) given two samples y1 = Y1[k], y2 = Y2[k] at different time points.
  ## This function is NOT vectorized.
  ## `y1` and `y2` are two **scalar** samples (of the same variable) at different time points.
  ## Returns a scalar.
  ## `param` is a single list containing the parameters for the specified kernel for this particular variable.
  ## For kernel="gaussian", it should include `bandwidth` as the bandwidth.
  ## For kernel="linear", no parameter is needed,
  ## For kernel="polynomial", it should include `intercept` and `degree`.
  ## For kernel="matern", it should include `lengthscale`.
  if ((length(y1) != 1) | (length(y2) != 1)) {stop("Incompatible input size")}

  kernel <- match.arg(kernel)  # if no kernel is specified, use gaussian kernel
  if (kernel == "gaussian"){
    if (is.null(param$bandwidth)){stop("Please specify bandwidth for gaussian kernel")}
    res <- gaussian_kernel(y1 = y1,
                           y2 = y2,
                           bandwidth=param$bandwidth)
  } else if (kernel == "linear"){
    res <- linear_kernel(y1 = y1,
                         y2 = y2)
  } else if (kernel == "polynomial"){
    if (is.null(param$intercept) | is.null(param$degree)){stop("Please specify intercept and degree for polynomial kernel")}
    res <- polynomial_kernel(y1 = y1,
                             y2 = y2,
                             intercept = param$intercept,
                             degree = param$degree)
  } else if (kernel == "matern"){
    if (is.null(param$lengthscale)){stop("Please specify lengthscale for matern kernel")}
    res <- matern_kernel(y1 = y1,
                         y2 = y2,
                         lengthscale = param$lengthscale)
  } else {
    stop("Does not support the specified kernel.")
    res <- NA
  }

  if (any(is.na(res) | is.nan(res))) {stop("NaN produced in kernel_main. Please check kernel computation.")}

  return (res)
}


kernel_inter <- function(Y1, Y2, k, l, kernel = c("gaussian", "linear", "polynomial", "matern"), kernel_params = NULL){
  ## Computes K_{kl}(Y1, Y2) (k is not equal to l) given two samples.
  ## This function is NOT vectorized.
  ## `Y1` and `Y2` are two **vector** samples of length p at different time points.
  ## Returns a scalar.
  ## `kernel_params` is a list of p lists. `kernel_params[[k]]` contains the parameters of variable `k` for the specified kernel. Check `kernel_main` to see the parameters needed for each kernel. The same kernel_main is used for all variables.
  if (length(Y1) != length(Y2)) {stop("Y1 and Y2 should contain the same variables.")}
  if (length(kernel_params) != length(Y1)) {stop("kernel_params should be of length p.")}
  if (k == l) {stop("Interaction terms are only computed between different variables.")}

  kernel <- match.arg(kernel)  # if no kernel is specified, use gaussian kernel
  val_k <- kernel_main(y1=Y1[k], y2=Y2[k], kernel=kernel, param=kernel_params[[k]])
  val_l <- kernel_main(y1=Y1[l], y2=Y2[l], kernel=kernel, param=kernel_params[[l]])
  res <- val_k * val_l

  if (any(is.na(res) | is.nan(res))) {stop("NaN produced in kernel_inter. Please check kernel computation.")}

  return(res)
}


kernel_theta <- function(theta_j, Y1, Y2, interaction_term, kernel = c("gaussian", "linear", "polynomial", "matern"), kernel_params = NULL){
  ## Computes K_{theta_j}(Y1, Y2) in the paper that generates the RKHS of F_j (which we estimate) given two samples. Formula specified by paragraph below Eq. (12).
  ## This function is NOT vectorized.
  ## `theta_j` is a vector with length p^2 if interaction_term = TRUE OR with length p if interaction_term = FALSE.
  ## `Y1` and `Y2` are two **vector** samples of length p at different time points.
  ## Returns a scalar.
  ## `kernel_params` is a list of p lists. `kernel_params[[k]]` contains the parameters of variable k for the specified kernel. Check `kernel_main` to see the parameters needed for each kernel.
  if (length(Y1) != length(Y2)) {stop("Y1 and Y2 should contain the same variables.")}
  if (length(kernel_params) != length(Y1)) {stop("kernel_params should be of length p.")}
  if (interaction_term & (length(theta_j) != length(Y1)^2)) {stop("interaction_term = TRUE: theta_j should have length p^2.")}
  if ((!interaction_term) & (length(theta_j) != length(Y1))) {stop("interaction_term = FALSE: theta_j should have length p.")}
  kernel <- match.arg(kernel)  # if no kernel is specified, use gaussian kernel

  p <- length(Y1)  # number of variables
  res_main <- sapply(1:p, function(k){
    kernel_main(y1=Y1[k], y2=Y2[k], kernel=kernel, param=kernel_params[[k]])
  })  # kernel values of the main effects

  if (!interaction_term){
    res <- sum(res_main * theta_j)
  } else {
    # compute kernel for interaction effect
    ## Note K_{kl} = K_k \times K_l for interaction between variable k and l
    mat <- res_main %o% res_main
    diag(mat) <- sqrt(diag(mat))  # recover main effects on the diagonal K_k
    res <- sum(mat * matrix(theta_j, nrow = p, ncol = p, byrow = TRUE))
  }

  if (any(is.na(res) | is.nan(res))) {stop("NaN produced in kernel_theta. Please check kernel computation.")}

  return (res)
}

