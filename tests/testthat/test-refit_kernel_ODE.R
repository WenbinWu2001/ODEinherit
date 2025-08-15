
test_that("refit_kernel_ODE: structure, types, and dimensions are correct", {
  set.seed(123)

  obs_time <- seq(0, 1, length.out = 20)
  Y <- cbind(
    sin(2 * pi * obs_time),
    cos(4 * pi * obs_time)
  ) + 0.05 * matrix(rnorm(20 * 2), 20, 2)
  tt <- seq(0, 1, length.out = 50)  # fine grid containing obs_time points

  res1 <- kernelODE_step1(Y = Y, obs_time = obs_time, tt = tt)

  kernel <- "gaussian"
  # Keep params simple and valid for p=3
  kernel_params <- replicate(ncol(Y), list(bandwidth = 0.5), simplify = FALSE)

  # Provide a valid adjacency
  adj <- matrix(c(0, 1,
                  1, 1),
                nrow = 2, ncol = 2, byrow = T)

  out <- refit_kernel_ODE(
    Y = Y,
    obs_time = obs_time,
    yy_smth = res1$yy_smth,
    tt = tt,
    kernel = kernel,
    kernel_params = kernel_params,
    interaction_term = FALSE,
    adj_matrix = adj,
    max_iter = 1,  # keep fast
    verbose = 0
  )

  # Top-level structure
  expect_type(out, "list")
  expect_named(out, c("metrics", "Y_refit"))

  # Y_refit checks
  expect_true(is.matrix(out$Y_refit))
  expect_true(is.numeric(out$Y_refit))
  expect_equal(dim(out$Y_refit), dim(Y))

  # metrics structure and types
  expect_type(out$metrics, "list")
  expect_named(out$metrics, c("R2", "R2_per_var_vec"))
  expect_true(is.numeric(out$metrics$R2))
  expect_true(is.numeric(out$metrics$R2_per_var_vec))
  expect_length(out$metrics$R2_per_var_vec, ncol(Y))

  # Sanity: R^2 in [0,1]
  expect_true(all(out$metrics$R2_per_var_vec >= 0 & out$metrics$R2_per_var_vec <= 1))
  expect_true(out$metrics$R2 >= 0 && out$metrics$R2 <= 1)

  # Consistency: recompute metrics and compare
  metrics_manual <- assess_recov_traj(Y = Y, Y_est = out$Y_refit)
  expect_equal(out$metrics$R2, metrics_manual$R2, tolerance = 1e-12)
  expect_equal(out$metrics$R2_per_var_vec, metrics_manual$R2_per_var_vec, tolerance = 1e-12)
})
