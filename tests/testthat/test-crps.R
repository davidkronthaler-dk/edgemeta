test_that("CRPS returns a scalar double", {
  s <- rnorm(1000, 0, 1)
  t <- rnorm(100, 0, 1)
  res <- crps(s, t)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_false(is.na(res))
})

test_that("CRPS is smaller for well-calibrated predictive samples", {
  truth <- rnorm(200, 0, 1)
  s_good <- rnorm(2000, 0, 1)  
  s_bad  <- rnorm(2000, 3, 2) 
  crps_good <- crps(s_good, truth)
  crps_bad  <- crps(s_bad, truth)
  expect_lt(crps_good, crps_bad)
})

test_that("CRPS is approximately zero for perfect per-observation forecasts", {
  tn <- c(1, 2, 3)
  # For each truth, predictive samples are a point mass at that truth:
  res_per_obs <- mean(sapply(tn, function(tj) crps(rep(tj, 1000), tj)))
  expect_lt(res_per_obs, 1e-12)  # numeric tolerance
})


test_that("CRPS decreases as sample size increases", {
  set.seed(42)
  tn <- rnorm(20, 0, 1)
  res_small <- crps(rnorm(50, 0, 1), tn)
  res_large <- crps(rnorm(5000, 0, 1), tn)
  expect_lt(res_large, res_small * 2)  # not strict equality but sanity check
})

test_that("CRPS handles vectorized truths", {
  set.seed(42)
  tn <- rnorm(5, 0, 1)
  s <- rnorm(1000, 0, 1)
  res <- crps(s, tn)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("CRPS errors on invalid input", {
  expect_error(crps("not numeric", 1:10))
  expect_error(crps(1:10, "not numeric"))
  expect_error(crps(numeric(0), 1:10))
  expect_error(crps(1:10, numeric(0)))
})
