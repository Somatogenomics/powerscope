test_that("t-test power equals power.t.test(strict = TRUE)", {
  for (n in c(3, 10, 64)) for (d in c(0.2, 0.8)) {
    expect_equal(ps_means(d, 1.3, n)$power,
                 power.t.test(n = n, delta = d, sd = 1.3, strict = TRUE)$power, tolerance = 1e-10)
    expect_equal(ps_means(d, 1.3, n, paired = TRUE)$power,
                 power.t.test(n = n, delta = d, sd = 1.3, type = "paired", strict = TRUE)$power,
                 tolerance = 1e-10)
  }
})

test_that("t-test exaggeration and Type S match simulated t-tests", {
  for (p in list(list(delta = 0.5, sd = 1, n1 = 5, n2 = 5),
                 list(delta = 0.2, sd = 1, n1 = 10, n2 = 15))) {
    a <- do.call(ps_means, p)
    s <- ps_simulate("means", p, n_sims = 1e5, seed = 3)
    expect_lt(abs(a$power - s$sim_power), 4 * s$mc_se_power)
    expect_equal(a$exaggeration, s$sim_exaggeration, tolerance = 0.03)
    expect_equal(a$type_s, s$sim_type_s, tolerance = 0.15)
  }
})

test_that("normal-approximation pooled props equal power.prop.test(strict = TRUE)", {
  for (n in c(20, 80, 300)) {
    expect_equal(ps_props(0.3, 0.45, n, exact = FALSE)$power,
                 power.prop.test(n = n, p1 = 0.3, p2 = 0.45, strict = TRUE)$power, tolerance = 1e-10)
  }
})

test_that("exact props match simulated z-tests, unequal groups", {
  p <- list(p1 = 0.1, p2 = 0.2, n1 = 100, n2 = 150)
  a <- do.call(ps_props, p); s <- ps_simulate("props", p, n_sims = 1e5, seed = 4)
  expect_match(a$method, "exact")
  expect_lt(abs(a$power - s$sim_power), 4 * s$mc_se_power)
  expect_equal(a$exaggeration, s$sim_exaggeration, tolerance = 0.01)
})

test_that("grid crossing and defaults", {
  r <- ps_means(c(0.2, 0.5), 1, c(10, 20, 30))
  expect_equal(nrow(r), 6L)
  expect_true(all(r$n1 == r$n2))
  expect_true(all(diff(r[delta == 0.5]$power) > 0))
})
