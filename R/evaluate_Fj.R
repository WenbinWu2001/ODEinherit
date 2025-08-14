#' Evaluate \eqn{F_j} and Recovery the Trajectory for a Single Variable
#'
#' Evaluates the derivative function \eqn{F_j} of a single variable on the
#' integration grid `tt` and recovers the trajectory by numerical integration,
#' given estimates of \eqn{b_j}, \eqn{c_j}, and \eqn{\theta_j}.
#'
#' @inheritParams kernelODE_step1
#' @inheritParams kernelODE_step2
#'
#' @param bj A numeric scalar, giving the estimated \eqn{b_j} from Kernel ODE
#'   (i.e., `res_bj[j]`).
#' @param cj A numeric vector of length `n`, giving the estimated \eqn{c_j} from
#'   Kernel ODE (i.e., `res_cj[,j]`).
#' @param kk_array Optional precomputed kernel array on `tt`. An array of
#'   dimension (`len`, `len`, `p`) when `interaction_term = FALSE`, or (`len`,
#'   `len`, `p^2`) when `interaction_term = TRUE`, where `len = length(tt)`.
#'   Providing `kk_array` enables reuse across variables and can greatly reduce
#'   computation.
#' @param theta_j A numeric vector of length `p` (if `interaction = FALSE`) or
#'   `p^2` (if `interaction = TRUE`), giving the estimated \eqn{\theta_j}
#'   coefficients for variable \eqn{j} from Kernel ODE (i.e., `res_theta[,j]`).
#' @param Yj A numeric vector of length `n`, giving the observed trajectory for variable
#'   \eqn{j} (i.e., `Y[, j]`).
#' @param yy_smth Numeric matrix of dimension (`len`, `p`); smoothed
#'   trajectories evaluated on `tt` (i.e., output of `kernelODE_step1()`).
#'
#'
#' @returns A list with components:
#' \describe{
#'   \item{`theta_j0`}{A numeric scalar giving estimated initial condition for variable \eqn{j}.}
#'   \item{`Fj_est`}{A numeric vector (length `len`) giving the evaluated \eqn{F_j} on `tt`.}
#'   \item{`yy_est`}{A numeric vector (length `len`) giving the recovered trajectory on `tt`.}
#'   \item{`TV_est`}{A numeric scalar giving the total variation \eqn{\int |F_j(t)| \, dt}
#'     approximated on `tt`.}
#' }
#'
#' @details Given \eqn{b_j}, \eqn{c_j}, and \eqn{\theta_j}, the function
#' constructs the kernel-weighted integral operator on the grid `tt` to evaluate
#' \eqn{F_j}, estimates the initial condition \eqn{\theta_{j0}}, then recovers
#' the trajectory via cumulative summation on `tt` (first-order approximation).
#' When provided, `kk_array` is reused to avoid recomputing kernel blocks.
#'
#' @seealso [kernelODE_step1()], [kernelODE_step2()], [assess_recov_traj()]
#'
#' @examples
#' set.seed(1)
#' obs_time <- seq(0, 1, length.out = 10)
#' Y <- cbind(sin(2 * pi * obs_time), cos(4 * pi * obs_time)) + 0.1 * matrix(rnorm(20), 10, 2)  # each col is a variable
#' tt <- seq(0, 1, length.out = 100)
#' res_step1 <- kernelODE_step1(Y = Y, obs_time = obs_time, tt = tt)
#'
#' kernel <- "gaussian"
#' kernel_params <- auto_select_kernel_params(kernel = kernel, Y = Y)
#' res_step2 <- kernelODE_step2(Y = Y, obs_time = obs_time, yy_smth = res_step1$yy_smth, tt = tt, kernel = kernel, kernel_params = kernel_params)
#'
#' j <- 1  # evaluate Fj for the first variable
#' res_eval <- evaluate_Fj(bj = res_step2$res_bj[j],
#'                         cj = res_step2$res_cj[,j],
#'                         interaction_term = F,
#'                         kernel = kernel,
#'                         kernel_params = kernel_params,
#'                         obs_time = obs_time,
#'                         theta_j = res_step2$res_theta[,j],
#'                         tt = tt,
#'                         Yj = Y[,j],
#'                         yy_smth = res_step1$yy_smth)
#' yy_j_est <- res_eval$yy_est
#'
#' # plot the evaluated traj
#' plot(NA, type = "n",
#'      xlab = "Time index", ylab = "Value",
#'      xlim = c(0,1), ylim = range(c(yy_j_est, Y[,j]), na.rm = T))
#' lines(obs_time, Y[,j], lty = 1)
#' lines(tt, yy_j_est, lty = 2)
#' legend("topright",
#'        legend = c("obs.", "eval."),
#'        lty = c(1,2),
#'        col = "black")
#'
#' @export
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
  n <- length(obs_time)
  p <- ncol(yy_smth)
  len <- length(tt)
  delta <- 1/len
  tt_mean <- .construct_tt_mean(obs_time, tt)

  if ((length(bj) != 1)) {stop("bj should be a scalar.")}
  if (length(cj) != n) {stop("cj should be a vector of length n.")}
  if (length(theta_j) != ifelse(interaction_term, p^2, p)) {stop("theta_j should be a vector of length p or p^2.")}
  if ((!is.null(kk_array)) & any(dim(kk_array) != c(len, len, length(theta_j)))) {stop("kk_array should be an array of dimension (length(tt), length(tt), p) or (length(tt), length(tt), p^2).")}

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
