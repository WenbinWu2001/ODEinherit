test_that("kernelODE_step2: runs on small dataset and output structure is correct", {
  set.seed(123)
  obs_time <- seq(0, 1, length.out = 10)
  Y <- cbind(sin(2 * pi * obs_time), cos(4 * pi * obs_time)) + 0.1 * matrix(rnorm(20), 10, 2)  # each col is a variable
  tt <- seq(0, 1, length.out = 20)
  yy_smth <- kernelODE_step1(Y = Y, obs_time = obs_time, tt = tt)$yy_smth
  kernel_params <- list(list(bandwidth = 1), list(bandwidth = 1))

  res <- kernelODE_step2(Y = Y,
                         obs_time = obs_time,
                         yy_smth = yy_smth,
                         tt = tt,
                         kernel = "gaussian",
                         kernel_params = kernel_params,
                         max_iter = 1,
                         tol = 1e-3,
                         eval_loss = FALSE)

  # Structure
  expect_type(res, "list")
  expect_named(res, c("res_theta", "res_best_kappa", "res_bj", "res_cj",
                      "res_best_eta", "res_loss_path", "network_est",
                      "num_iter", "config"))

  # Dimensions
  p <- ncol(Y)
  expect_equal(dim(res$res_theta)[2], p)
  expect_length(res$res_best_kappa, p)
  expect_length(res$res_bj, p)
  expect_equal(dim(res$res_cj), c(nrow(Y), p))
  expect_length(res$res_best_eta, p)
  expect_equal(dim(res$network_est), c(p, p))

  # Types
  expect_true(is.numeric(res$res_theta))
  expect_true(is.numeric(res$res_best_kappa))
  expect_true(is.numeric(res$res_bj))
  expect_true(is.numeric(res$res_best_eta))
  expect_true(is.matrix(res$res_cj))
  expect_true(is.matrix(res$network_est))
})

test_that("kernelODE_step2: respects adj_matrix input", {
  set.seed(321)
  obs_time <- seq(0, 1, length.out = 10)
  Y <- cbind(sin(2 * pi * obs_time), cos(4 * pi * obs_time)) + 0.1 * matrix(rnorm(20), 10, 2)  # each col is a variable
  tt <- seq(0, 1, length.out = 20)
  yy_smth <- kernelODE_step1(Y = Y, obs_time = obs_time, tt = tt)$yy_smth
  kernel_params <- list(list(bandwidth = 1), list(bandwidth = 1))
  adj_matrix <- matrix(c(1,0,0,1), nrow = 2)

  res <- kernelODE_step2(Y = Y,
                         obs_time = obs_time,
                         yy_smth = yy_smth,
                         tt = tt,
                         kernel = "gaussian",
                         kernel_params = kernel_params,
                         adj_matrix = adj_matrix,
                         max_iter = 1)

  expect_equal(dim(res$network_est), c(2, 2))
  expect_true(all(res$network_est %in% c(0, 1)))
})
