test_that("linear_kernel: correctness, type/length, NA warning, errors", {
  y1 <- c(1, 2, -3)
  y2 <- c(4, 0.5, 10)

  out <- linear_kernel(y1, y2)
  expect_equal(out, y1 * y2, tolerance = 1e-12)
  expect_true(is.numeric(out))
  expect_length(out, length(y1))

  # symmetry
  expect_equal(linear_kernel(y2, y1), out, tolerance = 1e-12)

  # NA in inputs -> NA in result and a warning (per new code)
  y1_na <- c(1, NA, 3)
  y2_na <- c(2,  5, 7)
  expect_warning(out_na <- linear_kernel(y1_na, y2_na), "result contains NA|NaN")
  expect_true(is.na(out_na[2]))
  expect_equal(out_na[c(1,3)], y1_na[c(1,3)] * y2_na[c(1,3)])

  # errors
  expect_error(linear_kernel(1:3, 1:2), "same length")
  expect_error(linear_kernel("a", 1:3), "numeric")
  expect_error(linear_kernel(1:3, "b"), "numeric")
})

test_that("polynomial_kernel: correctness, degree rules, type/length, NA warning, errors", {
  y1 <- c(1, 2, -3)
  y2 <- c(4, 0.5, 10)

  # basic
  out <- polynomial_kernel(y1, y2, intercept = 1, degree = 2)
  expect_equal(out, (y1 * y2 + 1)^2, tolerance = 1e-12)
  expect_true(is.numeric(out))
  expect_length(out, length(y1))

  # degree = 0 => ones
  out0 <- polynomial_kernel(y1, y2, intercept = 3, degree = 0)
  expect_equal(out0, rep(1, length(y1)))

  # NA propagation with warning
  expect_warning(
    out_na <- polynomial_kernel(c(1, NA, 3), c(2, 5, 7), intercept = 0, degree = 2),
    "result contains NA|NaN"
  )
  expect_true(is.na(out_na[2]))

  # errors on bad args
  expect_error(polynomial_kernel(1:3, 1:2, 1, 2), "same length")
  expect_error(polynomial_kernel("a", 1:3, 1, 2), "numeric")
  expect_error(polynomial_kernel(1:3, "b", 1, 2), "numeric")
  expect_error(polynomial_kernel(1:3, 1:3, intercept = NA, degree = 2), "intercept")
  expect_error(polynomial_kernel(1:3, 1:3, intercept = 1, degree = NA), "degree")
  expect_error(polynomial_kernel(1:3, 1:3, intercept = 1, degree = -1), "non-negative")
  expect_error(polynomial_kernel(1:3, 1:3, intercept = 1, degree = 2.5), "non-negative")
})

test_that("gaussian_kernel: correctness, type/length, monotonicity, NA warning, errors", {
  y1 <- c(0, 1, 2)
  y2 <- c(0, 3, 2)
  bw <- 1

  out <- gaussian_kernel(y1, y2, bandwidth = bw)
  manual <- exp(- (y1 - y2)^2 / (2 * bw^2))
  expect_equal(out, manual, tolerance = 1e-12)
  expect_true(all(out > 0 & out <= 1))
  expect_true(is.numeric(out))
  expect_length(out, length(y1))

  # symmetry
  expect_equal(gaussian_kernel(y2, y1, bw), out, tolerance = 1e-12)

  # monotonicity (farther => smaller)
  expect_true(gaussian_kernel(0, 0, 1) > gaussian_kernel(0, 2, 1))

  # NA propagation with warning
  expect_warning(
    out_na <- gaussian_kernel(c(0, NA), c(0, 1), 1),
    "result contains NA|NaN"
  )
  expect_true(is.na(out_na[2]))

  # errors
  expect_error(gaussian_kernel(1:3, 1:2, 1), "same length")
  expect_error(gaussian_kernel("a", 1:3, 1), "numeric")
  expect_error(gaussian_kernel(1:3, "b", 1), "numeric")
  expect_error(gaussian_kernel(1:3, 1:3, bandwidth = NA), "bandwidth")
  expect_error(gaussian_kernel(1:3, 1:3, bandwidth = 0), "positive")
  expect_error(gaussian_kernel(1:3, 1:3, bandwidth = -1), "positive")
})

test_that("matern_kernel (nu = 3/2): correctness, type/length, monotonicity, NA warning, errors", {
  y1 <- c(0, 1, 2)
  y2 <- c(0, 3, 2)
  ls <- 2

  out <- matern_kernel(y1, y2, lengthscale = ls)
  manual <- (1 + sqrt(3) * abs(y1 - y2) / ls) * exp(- sqrt(3) * abs(y1 - y2) / ls)
  expect_equal(out, manual, tolerance = 1e-12)
  expect_true(all(out > 0 & out <= 1))
  expect_true(is.numeric(out))
  expect_length(out, length(y1))

  # symmetry
  expect_equal(matern_kernel(y2, y1, ls), out, tolerance = 1e-12)

  # monotonicity
  expect_true(matern_kernel(0, 0, ls) > matern_kernel(0, 3, ls))

  # NA propagation with warning
  expect_warning(
    out_na <- matern_kernel(c(0, NA), c(0, 1), ls),
    "result contains NA|NaN"
  )
  expect_true(is.na(out_na[2]))

  # errors
  expect_error(matern_kernel(1:3, 1:2, ls), "same length")
  expect_error(matern_kernel("a", 1:3, ls), "numeric")
  expect_error(matern_kernel(1:3, "b", ls), "numeric")
  expect_error(matern_kernel(1:3, 1:3, lengthscale = NA), "lengthscale")
  expect_error(matern_kernel(1:3, 1:3, lengthscale = 0), "positive")
  expect_error(matern_kernel(1:3, 1:3, lengthscale = -1), "positive")
})
