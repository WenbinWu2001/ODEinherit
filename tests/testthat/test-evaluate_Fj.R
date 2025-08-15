# Helpers
expect_numeric_scalar <- function(x) {
  expect_true(is.numeric(x))
  expect_length(x, 1)
  expect_true(is.finite(x))
}
expect_numeric_vec_len <- function(x, len) {
  expect_true(is.numeric(x))
  expect_length(x, len)
  expect_true(all(is.finite(x)))
}

test_that("evaluate_Fj: res_eval has the correct structure, types, and lengths", {
  set.seed(1)

  # example in the doc
  obs_time <- seq(0, 1, length.out = 10)
  Y <- cbind(sin(2 * pi * obs_time), cos(4 * pi * obs_time)) + 0.1 * matrix(rnorm(20), 10, 2)  # each col is a variable
  tt <- seq(0, 1, length.out = 100)
  res_step1 <- kernelODE_step1(Y = Y, obs_time = obs_time, tt = tt)

  kernel <- "gaussian"
  kernel_params <- auto_select_kernel_params(kernel = kernel, Y = Y)
  res_step2 <- kernelODE_step2(Y = Y, obs_time = obs_time, yy_smth = res_step1$yy_smth, tt = tt, kernel = kernel, kernel_params = kernel_params)

  j <- 1  # evaluate Fj for the first variable
  res_eval <- evaluate_Fj(bj = res_step2$res_bj[j],
                          cj = res_step2$res_cj[,j],
                          interaction_term = F,
                          kernel = kernel,
                          kernel_params = kernel_params,
                          obs_time = obs_time,
                          theta_j = res_step2$res_theta[,j],
                          tt = tt,
                          Yj = Y[,j],
                          yy_smth = res_step1$yy_smth)
  yy_j_est <- res_eval$yy_est

  # Structure
  expect_type(res_eval, "list")
  expect_named(res_eval, c("theta_j0", "Fj_est", "yy_est", "TV_est", "tt"))

  # Types & lengths (formats)
  len <- length(tt)
  expect_numeric_scalar(res_eval$theta_j0)
  expect_numeric_vec_len(res_eval$Fj_est, len)
  expect_numeric_vec_len(res_eval$yy_est, len)
  expect_numeric_scalar(res_eval$TV_est)

  # Internal consistency
  expect_length(res_eval$Fj_est, length(res_eval$yy_est))
})
