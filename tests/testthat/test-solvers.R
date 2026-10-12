test_that("sample size: textbook n = 64 per group for d = 0.5", {
  r <- ps_sample_size(ps_means, "n1", 0.8, delta = 0.5, sd = 1, link = list(n2 = 1))
  expect_equal(r$n1, 64); expect_equal(r$n2, 64)
  below <- ps_means(0.5, 1, 63)$power
  expect_lt(below, 0.8); expect_gte(r$power, 0.8)
})

test_that("sample size for RNA-seq and unreachable target", {
  r <- ps_sample_size(ps_rnaseq, "n1", 0.9, log2fc = 1, mean_count = 50, dispersion = 0.1,
                      link = list(n2 = 1), n_max = 1000)
  expect_gte(r$power, 0.9)
  expect_lt(ps_rnaseq(1, 50, 0.1, r$n1 - 1)$power, 0.9)
  expect_warning(x <- ps_sample_size(ps_means, "n1", 0.8, delta = 1e-4, sd = 1, n_max = 100))
  expect_null(x)
})

test_that("minimum detectable effect reaches target power", {
  r <- ps_mde(ps_means, "delta", 0.8, sd = 1, n1 = 20, n2 = 20, interval = c(1e-4, 10))
  expect_equal(r$power, 0.8, tolerance = 1e-8)
  r2 <- ps_mde(ps_casecontrol, "odds_ratio", 0.8, maf = 0.3, n_cases = 2000, interval = c(1.0001, 5))
  expect_equal(r2$power, 0.8, tolerance = 1e-8)
})

test_that("FDR per-test alpha solves Jung's equation", {
  pf <- function(a) ps_means(delta = 1.5, sd = 1, n1 = 20, n2 = 20, alpha = a)$power
  a <- ps_fdr_alpha(pf, q = 0.05, m = 10000, pi1 = 0.05)
  et <- attr(a, "expected_true"); ef <- attr(a, "expected_false")
  expect_equal(ef / (et + ef), 0.05, tolerance = 1e-6)
  expect_warning(ps_fdr_alpha(function(a) a, q = 0.05, m = 1000, pi1 = 0.01))
})
