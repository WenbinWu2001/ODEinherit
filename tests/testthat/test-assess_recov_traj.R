test_that("assess_recov_traj: structure, types, and lengths are correct", {
  set.seed(1)
  n <- 20; p <- 3
  Y <- cbind(
    sin(seq(0, 1, length.out = n)),
    rnorm(n),
    seq_len(n)
  )
  Y_est <- Y + matrix(rnorm(n * p, sd = 0.1), n, p)

  res <- assess_recov_traj(Y, Y_est)

  expect_type(res, "list")
  expect_named(res, c("R2", "R2_per_var_vec"))
  expect_true(is.numeric(res$R2))
  expect_true(is.numeric(res$R2_per_var_vec))
  expect_length(res$R2_per_var_vec, p)
})

test_that("assess_recov_traj: perfect recovery gives R2 = 1 per variable and overall", {
  set.seed(2)
  n <- 15; p <- 2
  Y <- cbind(rnorm(n), runif(n))
  Y_est <- Y

  res <- assess_recov_traj(Y, Y_est)

  expect_equal(res$R2_per_var_vec, rep(1, p))
  expect_equal(res$R2, 1)
})

test_that("assess_recov_traj: predicting column means gives R2 = 0 (by truncation)", {
  set.seed(3)
  n <- 30; p <- 2
  Y <- cbind(rnorm(n), rnorm(n, mean = 3, sd = 2))
  Y_est <- sapply(seq_len(p), function(j) rep(mean(Y[, j]), n))

  res <- assess_recov_traj(Y, Y_est)

  expect_equal(res$R2_per_var_vec, rep(0, p))
  expect_equal(res$R2, 0)
})

test_that("assess_recov_traj: worse-than-mean predictions are truncated at 0", {
  set.seed(4)
  n <- 25; p <- 2
  Y <- cbind(rnorm(n), rnorm(n))
  # Very bad predictions: invert and add noise
  Y_est <- -(Y) + matrix(rnorm(n * p, sd = 5), n, p)

  res <- assess_recov_traj(Y, Y_est)

  expect_true(all(res$R2_per_var_vec >= 0))
  expect_true(res$R2 >= 0)
})

test_that("assess_recov_traj: numeric range is within [0, 1]", {
  set.seed(5)
  n <- 40; p <- 3
  Y <- matrix(rnorm(n * p), n, p)
  Y_est <- Y + matrix(rnorm(n * p, sd = 0.5), n, p)

  res <- assess_recov_traj(Y, Y_est)

  expect_true(all(res$R2_per_var_vec >= 0 & res$R2_per_var_vec <= 1))
  expect_true(res$R2 >= 0 && res$R2 <= 1)
})

test_that("assess_recov_traj: overall R2 equals mean of per-variable R2", {
  set.seed(6)
  n <- 18; p <- 4
  Y <- matrix(rnorm(n * p), n, p)
  Y_est <- Y + matrix(rnorm(n * p, sd = 0.3), n, p)

  res <- assess_recov_traj(Y, Y_est)
  expect_equal(res$R2, mean(res$R2_per_var_vec))
})


test_that("assess_recov_traj: vector inputs are reshaped to 1-column matrices", {
  set.seed(8)
  n <- 12
  Y <- rnorm(n)
  Y_est <- Y + rnorm(n, sd = 0.1)

  res <- assess_recov_traj(Y, Y_est)

  expect_true(is.numeric(res$R2))
  expect_length(res$R2_per_var_vec, 1)

  # Consistency with matrix version
  res2 <- assess_recov_traj(matrix(Y, ncol = 1), matrix(Y_est, ncol = 1))
  expect_equal(res$R2, res2$R2, tolerance = 1e-12)
  expect_equal(res$R2_per_var_vec, res2$R2_per_var_vec, tolerance = 1e-12)
})

test_that("assess_recov_traj: errors when dimensions differ", {
  set.seed(9)
  Y <- matrix(rnorm(20), 10, 2)
  Y_est_bad1 <- matrix(rnorm(15), 5, 3)
  Y_est_bad2 <- matrix(rnorm(30), 10, 3)
  Y_est_bad3 <- matrix(rnorm(10), 5, 2)

  expect_error(assess_recov_traj(Y, Y_est_bad1), "same dimension")
  expect_error(assess_recov_traj(Y, Y_est_bad2), "same dimension")
  expect_error(assess_recov_traj(Y, Y_est_bad3), "same dimension")
})
