test_that("Gelman, Hill & Vehtari Figure 16.1 values", {
  r <- ps_retrodesign(2, 8.1)
  expect_equal(r$power, 0.05701, tolerance = 1e-4)
  expect_equal(r$type_s, 0.2396, tolerance = 1e-3)
  expect_equal(r$exaggeration, 9.537, tolerance = 1e-3)
  expect_equal(r$min_sig_ratio, qnorm(0.975) * 8.1 / 2)
})

test_that("normal closed form equals the t integration at large df", {
  a <- ps_retrodesign(c(1, -0.5, 3), c(1, 0.4, 2))
  b <- ps_retrodesign(c(1, -0.5, 3), c(1, 0.4, 2), df = 1e7)
  expect_equal(a$power, b$power, tolerance = 1e-6)
  expect_equal(a$type_s, b$type_s, tolerance = 1e-5)
  expect_equal(a$exaggeration, b$exaggeration, tolerance = 1e-5)
})

test_that("engine agrees with its own simulation, including finite df", {
  for (d in list(c(2, 8.1, Inf), c(-1, 0.6, Inf), c(0.5, 0.4, 6))) {
    a <- ps_retrodesign(d[1], d[2], df = d[3])
    s <- ps_retrodesign_sim(d[1], d[2], df = d[3], n_sims = 4e5, seed = 11)
    expect_lt(abs(a$power - s$sim_power), 4 * sqrt(a$power * (1 - a$power) / 4e5))
    expect_equal(a$exaggeration, s$sim_exaggeration, tolerance = 0.02)
  }
})

test_that("sign symmetry", {
  a <- ps_retrodesign(1.3, 0.9); b <- ps_retrodesign(-1.3, 0.9)
  expect_equal(a[, .(power, type_s, exaggeration)], b[, .(power, type_s, exaggeration)])
})

test_that("input validation", {
  expect_error(ps_retrodesign(0, 1))
  expect_error(ps_retrodesign(1, 0))
  expect_error(ps_retrodesign(1, 1, alpha = 1))
})

test_that("ps_retrodesign_sim restores RNG state", {
  set.seed(5); x <- runif(1); set.seed(5)
  ps_retrodesign_sim(1, 1, n_sims = 10, seed = 99)
  expect_equal(runif(1), x)
})
