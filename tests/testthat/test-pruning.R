# Helpers
expect_numeric_matrix_dim <- function(x, nr, nc) {
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(nr, nc))
  expect_true(is.numeric(c(x)) || is.logical(c(x))) # allow 0/1 numeric or logical
}
expect_binary01_matrix <- function(x) {
  expect_true(all(na.omit(c(x)) %in% c(0, 1)))
}
expect_numeric_array_dim <- function(x, d1, d2, d3) {
  expect_true(is.array(x))
  expect_equal(dim(x), c(d1, d2, d3))
  expect_true(is.numeric(c(x)) || is.logical(c(x)))
}

test_that("prune_network (K=1): formats for network_pruned, R2_avg_mat_pruned, R2_multi_arr_pruned", {
  set.seed(1)
  # Example pipeline
  obs_time <- seq(0, 1, length.out = 10)
  Y <- cbind(
    sin(2 * pi * obs_time),
    cos(4 * pi * obs_time)
  ) + 0.1 * matrix(rnorm(20), 10, 2)
  tt <- seq(0, 1, length.out = 100)

  res_step1 <- kernelODE_step1(Y = Y, obs_time = obs_time, tt = tt)

  kernel <- "gaussian"
  kernel_params <- auto_select_kernel_params(kernel = kernel, Y = Y)

  res_step2 <- kernelODE_step2(
    Y = Y, obs_time = obs_time, yy_smth = res_step1$yy_smth, tt = tt,
    kernel = kernel, kernel_params = kernel_params,
    interaction_term = FALSE, max_iter = 1, tol = 1e-3
  )

  p <- ncol(Y)

  res_prune <- prune_network(
    network_original   = res_step2$network_est,
    prune_thres        = 0.05,
    depth              = 1,
    eval_edge_R2_pruned = TRUE,
    Y_list             = list(Y),
    yy_smth_list       = list(res_step1$yy_smth),
    obs_time_list      = list(obs_time),
    tt                 = tt,
    kernel             = kernel,
    kernel_params_list = list(kernel_params),
    interaction_term   = FALSE,
    theta_initial_list = list(matrix(1, nrow = p, ncol = p)),
    max_iter           = 1,
    tol                = 1e-3,
    parallel           = FALSE,
    verbose            = 0
  )

  # network_pruned
  expect_numeric_matrix_dim(res_prune$network_pruned, p, p)
  expect_binary01_matrix(res_prune$network_pruned)

  # R2_avg_mat_pruned: numeric matrix p x p (NA allowed)
  expect_numeric_matrix_dim(res_prune$R2_avg_mat_pruned, p, p)

  # R2_multi_arr_pruned: numeric array p x p x K (here K=1)
  expect_numeric_array_dim(res_prune$R2_multi_arr_pruned, p, p, 1)
})

test_that("prune_network (K=2): formats for network_pruned, R2_avg_mat_pruned, R2_multi_arr_pruned", {
  set.seed(2)
  # Two samples (possibly different n), same p
  obs_time1 <- seq(0, 1, length.out = 10)
  obs_time2 <- seq(0, 1, length.out = 12)
  p <- 2

  Y1 <- cbind(
    sin(2 * pi * obs_time1),
    cos(4 * pi * obs_time1)
  ) + 0.1 * matrix(rnorm(length(obs_time1) * p), length(obs_time1), p)
  Y2 <- cbind(
    sin(2 * pi * obs_time2),
    cos(4 * pi * obs_time2)
  ) + 0.1 * matrix(rnorm(length(obs_time2) * p), length(obs_time2), p)

  tt <- seq(0, 1, length.out = 100)

  res1 <- kernelODE_step1(Y = Y1, obs_time = obs_time1, tt = tt)
  res2 <- kernelODE_step1(Y = Y2, obs_time = obs_time2, tt = tt)

  kernel <- "gaussian"
  kp1 <- auto_select_kernel_params(kernel = kernel, Y = Y1)
  kp2 <- auto_select_kernel_params(kernel = kernel, Y = Y2)

  # Use a reasonable starting network from sample 1
  step2_1 <- kernelODE_step2(
    Y = Y1, obs_time = obs_time1, yy_smth = res1$yy_smth, tt = tt,
    kernel = kernel, kernel_params = kp1,
    interaction_term = FALSE, max_iter = 1, tol = 1e-3
  )

  # Provide theta_initial_list for both samples
  theta_initial_list <- list(matrix(1, nrow = p, ncol = p),
                             matrix(1, nrow = p, ncol = p))

  res_prune <- prune_network(
    network_original   = step2_1$network_est,
    prune_thres        = 0.05,
    depth              = 1,
    eval_edge_R2_pruned = TRUE,
    Y_list             = list(Y1, Y2),
    yy_smth_list       = list(res1$yy_smth, res2$yy_smth),
    obs_time_list      = list(obs_time1, obs_time2),
    tt                 = tt,
    kernel             = kernel,
    kernel_params_list = list(kp1, kp2),
    interaction_term   = FALSE,
    theta_initial_list = theta_initial_list,
    max_iter           = 1,
    tol                = 1e-3,
    parallel           = FALSE,
    verbose            = 0
  )

  # network_pruned
  expect_numeric_matrix_dim(res_prune$network_pruned, p, p)
  expect_binary01_matrix(res_prune$network_pruned)

  # R2_avg_mat_pruned: numeric matrix p x p
  expect_numeric_matrix_dim(res_prune$R2_avg_mat_pruned, p, p)

  # R2_multi_arr_pruned: numeric array p x p x K (K=2)
  expect_numeric_array_dim(res_prune$R2_multi_arr_pruned, p, p, 2)
})
