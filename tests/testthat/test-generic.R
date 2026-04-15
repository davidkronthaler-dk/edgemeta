es = c(0.17,  1.20,  1.10, -0.0019, -2.33)
se = c(0.52, 0.93, 0.63, 0.3, 0.28)

# Correct run
p1 <- PredDist(es = es, se  = se, method = "FullCD")
p2 <- PredDist(es = es, se  = se, method = "SimplifiedCD")
p3 <- PredDist(es = es, se  = se, method = "FixedTau2")

# Plot
test_that("Generic plotting works", {
  expect_no_error(plot(p1, type = "theta_new"))
  expect_no_error(plot(p1, type = "mu"))
  expect_no_error(plot(p1, type = "tau2"))
  expect_no_error(plot(p2, type = "theta_new"))
  expect_no_error(plot(p2, type = "mu"))
  expect_no_error(plot(p2, type = "tau2"))
  expect_no_error(plot(p3, type = "theta_new"))
  expect_no_error(plot(p3, type = "mu"))
  expect_warning(plot(p3, type = "tau2"))
})

# Summary
test_that("Generic summary works",{
  expect_no_error(summary(p1))
  expect_no_error(summary(p2))
  expect_no_error(summary(p3))
  expect_no_error(print(summary(p1)))
  expect_no_error(print(summary(p2)))
  expect_no_error(print(summary(p3)))
})

# Prob
test_that("Conf works", {
  expect_no_error(conf(p1, 0, 1))
  expect_no_error(conf(p2, -Inf, Inf))
  expect_error(conf(p3, a, b))
})
