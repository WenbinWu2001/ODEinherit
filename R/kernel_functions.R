#' Kernel Functions
#'
#' Computes common kernel functions *elementwise* for two equal-length numeric
#' vectors. Each function is vectorized: the result at position *i* depends only
#' on `y1[i]` and `y2[i]`.
#'
#' @param y1 A numeric vector.
#' @param y2 A numeric vector of the same length as `y1`.
#'
#' @return A numeric vector of length `length(y1)` giving the kernel values for
#'   each pair of corresponding elements in `y1` and `y2`.
#'
#' @details The following kernels are provided:
#' \itemize{
#'   \item \strong{Linear}: \deqn{K(y_1, y_2) = y_1 y_2}
#'   \item \strong{Polynomial}: \deqn{K(y_1, y_2) = (y_1 y_2 + c)^{d},}
#'         with intercept \eqn{c \in \mathbb{R}} and degree \eqn{d \in \mathbb{N}_0}.
#'   \item \strong{Gaussian (RBF)}: \deqn{K(y_1, y_2) = \exp\left(-\frac{(y_1 - y_2)^2}{2 \sigma^2}\right),}
#'         where \eqn{\sigma} is the bandwidth.
#'   \item \strong{Matern (\eqn{\nu = 3/2})}:
#'         \deqn{K(y_1, y_2) = \left(1 + \frac{\sqrt{3} |y_1 - y_2|}{\ell}\right) \exp\left(-\frac{\sqrt{3} |y_1 - y_2|}{\ell}\right),}
#'         with lengthscale \eqn{\ell > 0}.
#' }
#'
#' Inputs must have the same length; kernel-specific parameters must be valid
#' (e.g., positive bandwidth/lengthscale, non-negative integer degree). Inputs
#' containing `NA` will yield `NA` in the corresponding positions.
#'
#' @examples
#' linear_kernel(1:3, 4:6)
#' polynomial_kernel(1:3, 4:6, intercept = 1, degree = 2)
#' gaussian_kernel(1:3, 4:6, bandwidth = 1)
#' matern_kernel(1:3, 4:6, lengthscale = 0.5)
#'
#' @name kernel_functions
NULL


#' @rdname kernel_functions
#' @export
linear_kernel <- function(y1, y2){
  if (!is.numeric(y1) || !is.numeric(y2)) {stop("`y1` and `y2` must be numeric vectors.")}
  if (length(y1) != length(y2)) {stop("`y1` and `y2` must have the same length.")}

  res <- y1*y2
  return (res)
}


#' @rdname kernel_functions
#' @param intercept Numeric intercept \eqn{c}.
#' @param degree Non-negative integer degree \eqn{d}.
#' @export
polynomial_kernel <- function(y1, y2, intercept, degree){
  if (!is.numeric(y1) || !is.numeric(y2)) {stop("`y1` and `y2` must be numeric vectors.")}
  if (length(y1) != length(y2)) {stop("`y1` and `y2` must have the same length.")}
  if (is.na(intercept) || is.null(intercept)) {stop("`intercept` is NA or NULL.")}
  if (is.na(degree) || is.null(degree)) {stop("`degree` is NA or NULL.")}
  if (degree < 0 || degree %% 1 != 0) {stop("`degree` must be a non-negative integer.")}

  res <- (y1*y2 + intercept)^degree

  if (any(is.na(res) | is.nan(res))) {stop("Computing polynomial kernel: result contains NA or NaN values.")}

  return (res)
}


#' @rdname kernel_functions
#' @param bandwidth Positive numeric bandwidth \eqn{\sigma}.
#' @export
gaussian_kernel <- function(y1, y2, bandwidth){
  if (!is.numeric(y1) || !is.numeric(y2)) {stop("`y1` and `y2` must be numeric vectors.")}
  if (length(y1) != length(y2)) {stop("`y1` and `y2` must have the same length.")}
  if (is.na(bandwidth) || is.null(bandwidth)) {stop("`bandwidth` is NA or NULL.")}
  if (!is.numeric(bandwidth) || bandwidth <= 0) {stop("`bandwidth` must be a positive number.")}

  res <- exp(- (y1-y2)^2 / (2 * bandwidth^2))

  if (any(is.na(res) | is.nan(res))) {stop("Computing Gaussian kernel: result contains NA or NaN values.")}

  return (res)
}


#' @rdname kernel_functions
#' @param lengthscale Positive numeric lengthscale \eqn{\ell}.
#' @export
matern_kernel <- function(y1, y2, lengthscale){

  if (!is.numeric(y1) || !is.numeric(y2)) {stop("`y1` and `y2` must be numeric vectors.")}
  if (length(y1) != length(y2)) {stop("`y1` and `y2` must have the same length.")}
  if (is.na(lengthscale) || is.null(lengthscale)) {stop("`lengthscale` is NA or NULL.")}
  if (lengthscale <= 0) {stop("`lengthscale` must be positive.")}

  res <- (1 + sqrt(3)*abs(y1-y2)/lengthscale) * exp(- sqrt(3)*abs(y1-y2)/lengthscale)

  if (any(is.na(res) | is.nan(res))) {stop("Computing Matern kernel: result contains NA or NaN values.")}

  return (res)
}
