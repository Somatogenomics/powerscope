#' Monte Carlo check of a design by simulating the actual test
#'
#' Simulates data under a design, applies the test a user would run (with
#' estimated, not known, nuisance parameters), and returns empirical power,
#' Type S rate and exaggeration ratio. Used to validate the analytic results.
#'
#' Supported designs and the test simulated:
#' * `"means"`: two-sample pooled-variance t-test on normal data.
#' * `"props"`: pooled two-proportion z-test on binomial counts.
#' * `"rnaseq"`: NB counts; group means are the MLEs under known dispersion;
#'   Wald z on log fold change with plug-in variance.
#' * `"poolseq"`: binomial chromosomes into each pool, binomial reads; Wald z
#'   with plug-in variance (`target = "population"`, `overdispersion = 1`).
#' * `"casecontrol"`: binomial allele counts; Wald z on log odds ratio
#'   (0.5 added to all cells when any is zero).
#' * `"qtl"`: genotypes from HWE, normal residuals; `lm` slope t-test.
#'
#' @param design One of the names above.
#' @param params Named list of scalar design parameters (same names as the
#'   corresponding `ps_*` function).
#' @param n_sims Number of simulated studies.
#' @param seed Seed; the RNG state is restored afterwards.
#' @return One-row data.table: `sim_power`, `sim_type_s`, `sim_exaggeration`,
#'   `sim_n_sig`, `n_sims`, `mc_se_power`.
#' @examples
#' ps_simulate("means", list(delta = 0.5, sd = 1, n1 = 20, n2 = 20), n_sims = 2000, seed = 1)
#' @export
ps_simulate <- function(design = c("means", "props", "rnaseq", "poolseq", "casecontrol", "qtl"),
                        params, n_sims = 1e4, seed = 1) {
  design <- match.arg(design)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  P <- params
  a <- if (is.null(P$alpha)) (if (design %in% c("casecontrol", "qtl")) 5e-8 else 0.05) else P$alpha
  `%||%` <- function(x, y) if (is.null(x)) y else x

  r <- switch(design,
    means = {
      n1 <- P$n1; n2 <- P$n2 %||% n1
      x1 <- matrix(stats::rnorm(n_sims * n1, 0, P$sd), n_sims)
      x2 <- matrix(stats::rnorm(n_sims * n2, P$delta, P$sd), n_sims)
      m1 <- rowMeans(x1); m2 <- rowMeans(x2)
      sp <- sqrt(((n1 - 1) * apply(x1, 1, stats::var) + (n2 - 1) * apply(x2, 1, stats::var)) / (n1 + n2 - 2))
      est <- m2 - m1; stat <- est / (sp * sqrt(1 / n1 + 1 / n2))
      list(est = est, sig = abs(stat) >= stats::qt(1 - a / 2, n1 + n2 - 2), true = P$delta)
    },
    props = {
      n1 <- P$n1; n2 <- P$n2 %||% n1
      x1 <- stats::rbinom(n_sims, n1, P$p1); x2 <- stats::rbinom(n_sims, n2, P$p2)
      est <- x2 / n2 - x1 / n1; pb <- (x1 + x2) / (n1 + n2)
      se0 <- sqrt(pb * (1 - pb) * (1 / n1 + 1 / n2))
      list(est = est, sig = se0 > 0 & abs(est) >= stats::qnorm(1 - a / 2) * se0, true = P$p2 - P$p1)
    },
    rnaseq = {
      n1 <- P$n1; n2 <- P$n2 %||% n1; phi <- P$dispersion
      mu1 <- P$mean_count * (P$depth_factor %||% 1); mu2 <- mu1 * 2^P$log2fc
      rnb <- function(n, mu) if (phi > 0) stats::rnbinom(n, size = 1 / phi, mu = mu) else stats::rpois(n, mu)
      y1 <- rowMeans(matrix(rnb(n_sims * n1, mu1), n_sims))
      y2 <- rowMeans(matrix(rnb(n_sims * n2, mu2), n_sims))
      ok <- y1 > 0 & y2 > 0
      est <- ifelse(ok, log2(y2 / y1), 0)
      se <- sqrt((1 / pmax(y1, 1e-12) + phi) / n1 + (1 / pmax(y2, 1e-12) + phi) / n2) / log(2)
      list(est = est, sig = ok & abs(est) >= stats::qnorm(1 - a / 2) * se, true = P$log2fc)
    },
    poolseq = {
      k1 <- (P$ploidy %||% 2) * P$pool_size1; k2 <- (P$ploidy %||% 2) * (P$pool_size2 %||% P$pool_size1)
      C1 <- P$depth1; C2 <- P$depth2 %||% C1
      r1 <- P$n_pools1 %||% 1; r2 <- P$n_pools2 %||% r1
      draw <- function(p, k, C, r) {
        f <- matrix(stats::rbinom(n_sims * r, k, p) / k, n_sims)
        rowMeans(matrix(stats::rbinom(n_sims * r, C, f) / C, n_sims))
      }
      f1 <- draw(P$p1, k1, C1, r1); f2 <- draw(P$p2, k2, C2, r2)
      v <- function(p, k, C, r) p * (1 - p) * (1 / k + 1 / C - 1 / (k * C)) / r
      se <- sqrt(v(f1, k1, C1, r1) + v(f2, k2, C2, r2))
      est <- f2 - f1
      list(est = est, sig = se > 0 & abs(est) >= stats::qnorm(1 - a / 2) * se, true = P$p2 - P$p1)
    },
    casecontrol = {
      nca <- 2 * P$n_cases; nco <- 2 * (P$n_controls %||% P$n_cases)
      p0 <- P$maf; p1 <- P$odds_ratio * p0 / (1 + p0 * (P$odds_ratio - 1))
      a1 <- stats::rbinom(n_sims, nca, p1); a0 <- stats::rbinom(n_sims, nco, p0)
      cells <- cbind(a1, nca - a1, a0, nco - a0)
      cells <- cells + 0.5 * (apply(cells, 1, min) == 0)
      est <- log(cells[, 1] * cells[, 4] / (cells[, 2] * cells[, 3]))
      se <- sqrt(rowSums(1 / cells))
      list(est = est, sig = abs(est) >= stats::qnorm(1 - a / 2) * se, true = log(P$odds_ratio))
    },
    qtl = {
      n <- P$n; het <- 2 * P$maf * (1 - P$maf); r2 <- het * P$beta_sd^2
      est <- sig <- numeric(n_sims)
      for (i in seq_len(n_sims)) {
        g <- stats::rbinom(n, 2, P$maf)
        y <- P$beta_sd * g + stats::rnorm(n, 0, sqrt(1 - r2))
        cf <- summary(stats::lm(y ~ g))$coefficients
        if (nrow(cf) < 2) { est[i] <- 0; sig[i] <- FALSE; next }
        est[i] <- cf[2, 1]; sig[i] <- abs(cf[2, 3]) >= stats::qt(1 - a / 2, n - 2)
      }
      list(est = est, sig = as.logical(sig), true = P$beta_sd)
    }
  )
  sig <- r$sig; est <- r$est
  pw <- mean(sig)
  data.table::data.table(
    sim_power = pw,
    sim_type_s = if (any(sig)) mean(sign(est[sig]) != sign(r$true)) else NA_real_,
    sim_exaggeration = if (any(sig)) mean(abs(est[sig])) / abs(r$true) else NA_real_,
    sim_n_sig = sum(sig), n_sims = n_sims,
    mc_se_power = sqrt(pw * (1 - pw) / n_sims)
  )
}
