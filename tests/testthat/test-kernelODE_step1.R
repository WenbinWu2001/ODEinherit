test_that("kernelODE_step1: output structure, dimensions, and types are correct", {
  set.seed(123)
  obs_time <- seq(0, 1, length.out = 10)
  p <- 3
  Y <- cbind(
    sin(2 * pi * obs_time),
    cos(4 * pi * obs_time),
    obs_time^2
  ) + 0.05 * matrix(rnorm(10 * p), 10, p)
  tt <- seq(0, 1, length.out = 50)

  result <- kernelODE_step1(Y = Y, obs_time = obs_time, tt = tt)

  # Structure check
  expect_type(result, "list")
  expect_named(result, c("yy_smth", "init_vals_smth", "deriv_smth"))

  # yy_smth
  expect_true(is.matrix(result$yy_smth))
  expect_type(result$yy_smth, "double")
  expect_equal(dim(result$yy_smth), c(length(tt), p))

  # init_vals_smth
  expect_true(is.numeric(result$init_vals_smth))
  expect_length(result$init_vals_smth, p)

  # deriv_smth
  expect_true(is.matrix(result$deriv_smth))
  expect_type(result$deriv_smth, "double")
  expect_equal(dim(result$deriv_smth), c(length(tt), p))
})

test_that("kernelODE_step1: errors on invalid inputs", {
  obs_time <- seq(0, 1, length.out = 10)
  tt <- seq(0, 1, length.out = 50)
  Y <- matrix(rnorm(20), 10, 2)

  # Wrong dimensions
  expect_error(kernelODE_step1(Y[1:5, ], obs_time, tt))  # nrow mismatch
})
