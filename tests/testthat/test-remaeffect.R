# Valid estimates and standard errors
es = c(0.17,  1.20,  1.10, -0.0019, -2.33)
se = c(0.52, 0.93, 0.63, 0.3, 0.28)

# Correct run
expect_no_error(remaeffect(es = es, se  = se))

# Test missing argument
expect_error(remaeffect(es = es))
expect_error(remaeffect(se = se))
expect_error(remaeffect())

# Invalid estimated and standard errors
xes0 = c(0.2)
xse0 = c(0.1)
xes1 = c(0.17,  1.20,  1.10, -0.0019, -2.33)
xse1 = c(0.52, 0.93, 0.63, 0.3)
xes2 = c(-0.52, 0.93, 0.63, 0.3, 0.28)
xse2 = c(-0.52, 0.93, 0.63, 0.3, 0.28)

# Test wrong specifications
expect_error(remaeffect(es = xes0, se  = xse0))
expect_error(remaeffect(es = xes1, se  = xse1))
expect_error(remaeffect(es = xes2, se  = xse2))

# Level of prediction interval
test_that("level.ci must be between 0 and 1", {
  expect_error(remaeffect(es = es, se = se, level.ci = 1.1))
  expect_error(remaeffect(es = es, se = se, level.ci = -0.1))
  expect_no_error(remaeffect(es = es, se = se, level.ci = 0.95))
})

# Number of samples
test_that("n_samples must be positive integer for FullCD", {
  expect_error(remaeffect(es = es, se = se, n_samples = -100))
  expect_error(remaeffect(es = es, se = se, n_samples = "ten"))
  expect_no_error(remaeffect(es = es, se = se, n_samples = 100))
})

# Return
test_that("returned object has expected structure", {
  result1 <- remaeffect(es = es, se = se)
  expect_named(result1, c("estimate", "CI", "pval", "cd_mu", "cd_tau2"))
})
