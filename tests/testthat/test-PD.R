# Valid estimates and standard errors
es = c(0.17,  1.20,  1.10, -0.0019, -2.33)
se = c(0.52, 0.93, 0.63, 0.3, 0.28)

# Correct run
expect_no_error(PredDist(es = es, se  = se, method = "FullCD"))
expect_no_error(PredDist(es = es, se  = se, method = "SimplifiedCD"))
expect_no_error(PredDist(es = es, se  = se, method = "FixedTau2"))

# Test missing argument
expect_error(PredDist(es = es))
expect_error(PredDist(se = se))
expect_error(PredDist(es = es,  se = se))
expect_error(PredDist(method = "FullCD"))

# Invalid estimated and standard errors
xes0 = c(0.2)
xse0 = c(0.1)
xes1 = c(0.17,  1.20,  1.10, -0.0019, -2.33)
xse1 = c(0.52, 0.93, 0.63, 0.3)
xes2 = c(-0.52, 0.93, 0.63, 0.3, 0.28)
xse2 = c(-0.52, 0.93, 0.63, 0.3, 0.28)

# Test wrong specifications
expect_error(PredDist(es = xes0, se  = xse0, method = "FullCD"))
expect_error(PredDist(es = xes1, se  = xse1, method = "FullCD"))
expect_error(PredDist(es = xes2, se  = xse2, method = "FullCD"))

# Choice of method
test_that("method argument must be one of the allowed choices", {
  expect_error(PredDist(es = es, se = se, method = "invalid"))
  expect_no_error(PredDist(es = es, se = se, method = "FullCD"))
  expect_no_error(PredDist(es = es, se = se, method = "SimplifiedCD"))
  expect_no_error(PredDist(es = es, se = se, method = "FixedTau2"))
})

# Level of prediction interval
test_that("level.pi must be between 0 and 1", {
  expect_error(PredDist(es = es, se = se, method = "FullCD", level.pi = 1.1))
  expect_error(PredDist(es = es, se = se, method = "FullCD", level.pi = -0.1))
  expect_no_error(PredDist(es = es, se = se, method = "FullCD", level.pi = 0.95))
})

# Number of samples
test_that("n_samples must be positive integer for FullCD", {
  expect_error(PredDist(es = es, se = se, method = "FullCD", n_samples = -100))
  expect_error(PredDist(es = es, se = se, method = "FullCD", n_samples = "ten"))
  expect_no_error(PredDist(es = es, se = se, method = "FullCD", n_samples = 100))
})

# Theta New
test_that("theta_new input works with FixedTau2 and SimplifiedCD", {
  theta = seq(-5, 5, length.out = 100)
  expect_no_error(PredDist(es = es, se = se, method = "FixedTau2", theta_new = theta))
  expect_no_error(PredDist(es = es, se = se, method = "SimplifiedCD", theta_new = theta))
})

# Different tau2 methods
test_that("method.tau2 variations for FixedTau2 work", {
  expect_no_error(PredDist(es = es, se = se, method = "FixedTau2", method.tau2 = "REML"))
  expect_no_error(PredDist(es = es, se = se, method = "FixedTau2", method.tau2 = "DL"))
  expect_no_error(PredDist(es = es, se = se, method = "FixedTau2", method.tau2 = "SJ"))
  expect_error(PredDist(es = es, se = se, method = "FixedTau2", method.tau2 = "invalid"))
})

# Return
test_that("returned object has expected structure", {
  result1 <- PredDist(es = es, se = se, method = "FullCD")
  expect_named(result1, c("PI", "samples"))
  result2 <- PredDist(es = es, se = se, method = "SimplifiedCD")
  expect_named(result2, c("PI", "PD.theta_new", "fPD"))
  result3 <- PredDist(es = es, se = se, method = "FixedTau2")
  expect_named(result3, c("PI", "PD.theta_new", "fPD"))
})

# Subdivisions
test_that("integration subdivisions argument is respected", {
  expect_no_error(PredDist(es = es, se = se, method = "FixedTau2", subdivisions = 100))
  expect_error(PredDist(es = es, se = se, method = "FixedTau2", subdivisions = -10))
})











