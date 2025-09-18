es <- c(0.17,  1.20,  1.10, -0.0019, -2.33)
se <- c(0.52, 0.93, 0.63, 0.3, 0.28)

test_that("correct execution for key configurations", {
  conditions <- list(
    list(method = "MC",    level.ci = 0.05, n_samples = 50,   seed = NULL, mu0 = 0,    method.tau2 = "REML"),
    list(method = "GAQ",   level.ci = 0.5,  n_samples = 1000, seed = 5,   mu0 = 1,    method.tau2 = "PM"),
    list(method = "NHEU",  level.ci = 0.99, n_samples = 50,   seed = NULL, mu0 = -Inf, method.tau2 = "REML"),
    list(method = "MC",    level.ci = 0.5,  n_samples = 1000, seed = 5,   mu0 = Inf,  method.tau2 = "PM"),
    list(method = "GAQ",   level.ci = 0.01, n_samples = 50,   seed = NULL, mu0 = 0,    method.tau2 = "REML")
  )
  
  for (cond in conditions) {
    expect_no_error(
      remaeffect(
        es = es,
        se = se,
        method = cond$method,
        level.ci = cond$level.ci,
        n_samples = cond$n_samples,
        seed = cond$seed,
        mu0 = cond$mu0,
        method.tau2 = cond$method.tau2
      )
    )
  }
  
  r1id <- remaeffect(es = es, se = se, method = "MC", seed = 999)
  r2id <- remaeffect(es = es, se = se, method = "MC", seed = 999)
  expect_identical(r1id, r2id)
})

test_that("handles invalid input types", {
  expect_error(remaeffect(es = "not numeric", se = 1:5))
  expect_error(remaeffect(es = 1:5, se = "not numeric"))
  expect_error(remaeffect(es = 1:5, se = 1:5, method = "INVALID"))
  expect_error(remaeffect(es = 1:5, se = 1:5, method = "GAQ", method.tau2 = "INVALID"))
  expect_error(remaeffect(es = numeric(0), se = 1:5))
  expect_error(remaeffect(es = 1:5, se = numeric(0)))
  expect_error(remaeffect(c(1, NA), c(1, 0.5), "MC"))
  expect_error(remaeffect(c(1, "a"), c(1,1), "MC"))
  expect_error(remaeffect(es, se, NULL))
  expect_error(remaeffect(es, se, NA))
  expect_error(remaeffect(es, se, "NHEU", method.tau2 = NA))
  expect_error(remaeffect(es, se, "NHEU", method.tau2 = NULL))
  expect_error(remaeffect(es, se, "MC", n_samples = 1.5))
  expect_error(remaeffect(es, se, "MC", n_samples = -Inf))
  expect_warning(remaeffect(es, se, "MC", seed = -100))
  expect_warning(expect_warning(expect_error(remaeffect(es, se, "MC", seed = Inf))))
  expect_warning(expect_warning(expect_error(remaeffect(es, se, "MC", seed = "hundred"))))
  expect_error(remaeffect(es = c(0.1, Inf, 0.3), se = c(0.1, 0.2, 0.3)))
  expect_error(remaeffect(es = c(0.1, 0.2, 0.3), se = c(0.1, Inf, 0.3)))
})

test_that("handles mu0 at boundary values", {
  es <- rnorm(5)
  se <- runif(5, 0.1, 2)
  expect_no_error(remaeffect(es, se, method = "GAQ", mu0 = -1e7))
  expect_no_error(remaeffect(es, se, method = "MC", mu0 = Inf))
  expect_no_error(remaeffect(es, se, method = "NHEU", mu0 = 0))
})

test_that("tolerates small sample sizes", {
  es <- rnorm(2)
  se <- runif(2, 0.1, 0.5)
  expect_no_error(remaeffect(es, se, method = "MC", n_samples = 5))
  expect_no_error(remaeffect(es, se, method = "GAQ"))
  expect_no_error(remaeffect(es, se, method = "NHEU"))
  expect_error(remaeffect(es[1], se[1], method = "MC"))
})

test_that("works with different confidence levels and seeds", {
  es <- rnorm(5)
  se <- runif(5, 0.1, 2)
  levels <- c(0.01, 0.5, 0.99)
  seeds <- list(NULL, 123, 999)
  for (lvl in levels) {
    for (s in seeds) {
      expect_no_error(
        remaeffect(es, se, level.ci = lvl, seed = s, method = "MC")
      )
    }
  }
})

test_that("fails when es and se have different lengths", {
  es <- rnorm(5)
  se <- runif(6, 0.1, 2)
  expect_error(remaeffect(es, se))
})

test_that("handles method.tau2 choices", {
  es <- rnorm(5)
  se <- runif(5, 0.1, 2)
  expect_no_error(remaeffect(es, se, method = "NHEU", method.tau2 = "REML"))
  expect_no_error(remaeffect(es, se, method = "NHEU", method.tau2 = "PM"))
})

