es = c(0.17,  1.20,  1.10, -0.0019, -2.33)
se = c(0.52, 0.93, 0.63, 0.3, 0.28)

test_that("correct execution for key configurations", {
  conditions <- list(
    list(method = "FullCD",       level.pi = 0.05, n_samples = 50,   seed = NULL,     method.tau2 = "REML"),
    list(method = "SimplifiedCD", level.pi = 0.5,  n_samples = 1000, seed = 5,        method.tau2 = "PM"),
    list(method = "FixedTau2",    level.pi = 0.99, n_samples = 50,   seed = NULL,  method.tau2 = "REML"),
    list(method = "FullCD",       level.pi = 0.5,  n_samples = 1000, seed = 5,     method.tau2 = "PM"),
    list(method = "FixedTau2",    level.pi = 0.01, n_samples = 50,   seed = NULL,    method.tau2 = "REML")
  )
  
  for (cond in conditions) {
    expect_no_error(
      PredDist(
        es = es,
        se = se,
        method = cond$method,
        level.pi = cond$level.pi,
        n_samples = cond$n_samples,
        seed = cond$seed,
        method.tau2 = cond$method.tau2
      )
    )
  }
  
  r1id <- PredDist(es = es, se = se, method = "FullCD", seed = 999)
  r2id <- PredDist(es = es, se = se, method = "FullCD", seed = 999)
  expect_identical(r1id, r2id)
})

test_that("handles invalid input types", {
  expect_error(PredDist(es = "not numeric", se = 1:5))
  expect_error(PredDist(es = 1:5, se = "not numeric"))
  expect_error(PredDist(es = 1:5, se = 1:5, method = "INVALID"))
  expect_error(PredDist(es = 1:5, se = 1:5, method = "FullCD", method.tau2 = "INVALID"))
  expect_error(PredDist(es = numeric(0), se = 1:5))
  expect_error(PredDist(es = 1:5, se = numeric(0)))
  expect_error(PredDist(c(1, NA), c(1, 0.5)))
  expect_error(PredDist(c(1, "a"), c(1,1)))
  expect_error(PredDist(es, se, NULL))
  expect_error(PredDist(es, se, NA))
  expect_error(PredDist(es, se, method.tau2 = NA))
  expect_error(PredDist(es, se, method.tau2 = NULL))
  expect_error(PredDist(es, se, n_samples = 1.5))
  expect_error(PredDist(es, se, n_samples = -Inf))
  expect_warning(PredDist(es, se, seed = -100))
  expect_warning(expect_warning(expect_error(PredDist(es, se, seed = Inf))))
  expect_warning(expect_warning(expect_error(PredDist(es, se, seed = "hundred"))))
  expect_error(PredDist(es = c(0.1, Inf, 0.3), se = c(0.1, 0.2, 0.3)))
  expect_error(PredDist(es = c(0.1, 0.2, 0.3), se = c(0.1, Inf, 0.3)))
})

test_that("tolerates small sample sizes", {
  set.seed(1)
  es <- rnorm(2)
  se <- runif(2, 0.1, 1.5)
  expect_no_error(PredDist(es, se, method = "FullCD", n_samples = 5))
  expect_no_error(PredDist(es, se, method = "SimplifiedCD"))
  expect_error(PredDist(es[1], se[1], method = "FullCD"))
})

test_that("works with different prediction levels and seeds", {
  es <- rnorm(5)
  se <- runif(5, 0.1, 2)
  levels <- c(0.01, 0.5, 0.99)
  seeds <- list(NULL, 123, 999)
  for (lvl in levels) {
    for (s in seeds) {
      expect_no_error(
        PredDist(es, se, level.pi = lvl, seed = s, method = "FullCD")
      )
    }
  }
})

test_that("fails when es and se have different lengths", {
  es <- rnorm(5)
  se <- runif(6, 0.1, 2)
  expect_error(PredDist(es, se))
})

test_that("handles method.tau2 choices", {
  es <- rnorm(5)
  se <- runif(5, 0.1, 2)
  expect_no_error(PredDist(es, se, method.tau2 = "REML"))
  expect_no_error(PredDist(es, se, method = "FixedTau2", method.tau2 = "PM"))
})

# Return
test_that("returned object has expected structure", {
  result1 <- PredDist(es = es, se = se, method = "FullCD")
  expect_named(result1, c("PI", "samples"))
  result2 <- PredDist(es = es, se = se, method = "SimplifiedCD")
  expect_named(result2, c("PI", "samples"))
  result3 <- PredDist(es = es, se = se, method = "FixedTau2")
  expect_named(result3, c("PI", "samples"))
})

