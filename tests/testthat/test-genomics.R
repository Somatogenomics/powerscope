close_power <- function(a, s, k = 4, abs_tol = 0) {
  expect_lt(abs(a$power - s$sim_power), k * s$mc_se_power + abs_tol)
}

test_that("RNA-seq NB power matches simulated Wald tests", {
  for (p in list(list(log2fc = 1, mean_count = 50, dispersion = 0.1, n1 = 3, n2 = 3),
                 list(log2fc = -0.5, mean_count = 200, dispersion = 0.05, n1 = 6, n2 = 6),
                 list(log2fc = 0.6, mean_count = 10, dispersion = 0, n1 = 4, n2 = 4))) {
    a <- do.call(ps_rnaseq, p); s <- ps_simulate("rnaseq", p, n_sims = 5e4, seed = 7)
    close_power(a, s)
    expect_equal(a$exaggeration, s$sim_exaggeration, tolerance = 0.02)
  }
})

test_that("RNA-seq depth_factor is equivalent to scaling mean_count", {
  expect_equal(ps_rnaseq(1, 20, 0.1, 4, depth_factor = 2)$power,
               ps_rnaseq(1, 40, 0.1, 4)$power)
})

test_that("ps_rnaseq_genes FDR solution gives FDR q in expectation", {
  set.seed(2)
  genes <- data.frame(mean_count = rlnorm(1500, log(100), 1.5),
                      dispersion = rlnorm(1500, log(0.1), 0.6))
  r <- ps_rnaseq_genes(genes, log2fc = 1, n1 = 5, q = 0.05, pi1 = 0.1)
  s <- r$summary
  expect_equal(s$expected_false_discoveries /
                 (s$expected_false_discoveries + s$expected_true_discoveries), 0.05, tolerance = 1e-6)
  expect_equal(nrow(r$per_gene), 1500L)
  expect_error(ps_rnaseq_genes(genes, 1, 5, q = 0.05))
})

test_that("pool-seq variance and simulation", {
  r <- ps_poolseq(0.3, 0.5, pool_size1 = 20, depth1 = 50)
  v <- function(p) p * (1 - p) * (1 / 40 + 1 / 50 - 1 / 2000)
  expect_equal(r$se, sqrt(v(0.3) + v(0.5)))
  expect_equal(ps_poolseq(0.3, 0.5, 20, depth1 = 50, target = "pool")$se,
               sqrt(0.21 / 50 + 0.25 / 50))
  expect_gt(ps_poolseq(0.3, 0.5, 20, depth1 = 50, overdispersion = 2)$se, r$se)
  p <- list(p1 = 0.2, p2 = 0.4, pool_size1 = 30, depth1 = 40, n_pools1 = 3)
  a <- do.call(ps_poolseq, p); s <- ps_simulate("poolseq", p, n_sims = 5e4, seed = 8)
  # plug-in Wald differs slightly from the analytic expectation (see validation/)
  close_power(a, s, abs_tol = 0.015)
})

test_that("case-control and QTL against simulation", {
  p <- list(maf = 0.2, odds_ratio = 1.3, n_cases = 1000, n_controls = 1000, alpha = 0.05)
  a <- do.call(ps_casecontrol, p); s <- ps_simulate("casecontrol", p, n_sims = 5e4, seed = 9)
  close_power(a, s, abs_tol = 0.01)
  q <- list(maf = 0.3, beta_sd = 0.15, n = 300, alpha = 0.05)
  a <- do.call(ps_qtl, q); s <- ps_simulate("qtl", q, n_sims = 4000, seed = 10)
  close_power(a, s)
  expect_equal(ps_qtl(0.5, 0.2, 100)$r2, 0.5 * 0.04)
  expect_equal(ps_casecontrol(0.2, 2, 100)$case_freq, 0.4 / 1.2)
})

test_that("genome-wide default alpha for association", {
  expect_equal(ps_casecontrol(0.2, 1.2, 1000)$crit, qnorm(1 - 2.5e-8))
})
