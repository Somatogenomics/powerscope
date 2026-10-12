#' Case-control allelic association power
#'
#' Allelic test on `2 * n` chromosomes per group under Hardy-Weinberg
#' equilibrium and a multiplicative (per-allele) model. The effect is the
#' log allelic odds ratio, with Wald SE
#' `sqrt(1/(2 n_cases p1) + 1/(2 n_cases (1-p1)) + 1/(2 n_controls p0) + 1/(2 n_controls (1-p0)))`,
#' where `p0` is the control risk-allele frequency and `p1` the case
#' frequency implied by the odds ratio.
#'
#' The default `alpha = 5e-8` is the human genome-wide convention (Pe'er et
#' al. 2008). It depends on the number of effectively independent tests, so
#' for other species or panels derive it from your marker set (e.g.
#' Bonferroni on the effective number of tests after LD pruning).
#'
#' @param maf Risk-allele frequency in controls.
#' @param odds_ratio Allelic odds ratio (not 1).
#' @param n_cases,n_controls Numbers of diploid individuals.
#' @param alpha Two-sided significance level.
#' @return A `ps_result` data.table; estimate scale is log odds ratio.
#' @references Pe'er I, Yelensky R, Altshuler D, Daly MJ (2008). Estimation of
#'   the multiple testing burden for genomewide association studies of nearly
#'   all common variants. *Genetic Epidemiology* 32:381-385.
#'   \doi{10.1002/gepi.20303}
#' @examples
#' ps_casecontrol(maf = 0.2, odds_ratio = 1.3, n_cases = c(1000, 3000), n_controls = 3000)
#' @export
ps_casecontrol <- function(maf, odds_ratio, n_cases, n_controls = NULL, alpha = 5e-8) {
  .check_prob(maf, "maf"); .check_pos(odds_ratio, "odds_ratio"); .check_pos(n_cases, "n_cases")
  if (any(odds_ratio == 1)) stop("`odds_ratio` must differ from 1")
  g <- .ps_grid(list(maf = maf, odds_ratio = odds_ratio, n_cases = n_cases,
                     n_controls = n_controls, alpha = alpha))
  if (is.null(n_controls)) g[, n_controls := n_cases]
  p0 <- g$maf
  p1 <- g$odds_ratio * p0 / (1 + p0 * (g$odds_ratio - 1))
  g[, case_freq := p1]
  se <- sqrt(1 / (2 * g$n_cases * p1) + 1 / (2 * g$n_cases * (1 - p1)) +
             1 / (2 * g$n_controls * p0) + 1 / (2 * g$n_controls * (1 - p0)))
  out <- .ps_finish(g, "casecontrol_allelic", "log odds ratio", log(g$odds_ratio), se, Inf, g$alpha)
  out[, method := "allelic Wald test, HWE, multiplicative model"]
  out[]
}

#' Quantitative-trait association power (additive SNP)
#'
#' Linear regression of a trait on allele count (0, 1, 2) at a biallelic SNP
#' in Hardy-Weinberg equilibrium. The effect `beta_sd` is the per-allele
#' change in trait units of the trait's total SD, so the variance explained is
#' `r2 = 2 p (1-p) beta_sd^2` and
#' `SE(beta) = sqrt((1 - r2) / (n * 2 p (1-p)))`, with t reference on
#' `n - 2` df. Default `alpha` as in [ps_casecontrol()].
#'
#' @param maf Minor (or effect) allele frequency.
#' @param beta_sd Per-allele effect in trait-SD units.
#' @param n Number of individuals.
#' @param alpha Two-sided significance level.
#' @return A `ps_result` data.table with `r2`; estimate scale is per-allele
#'   effect in trait SD.
#' @examples
#' ps_qtl(maf = 0.3, beta_sd = 0.1, n = c(1000, 5000, 20000))
#' @export
ps_qtl <- function(maf, beta_sd, n, alpha = 5e-8) {
  .check_prob(maf, "maf"); .check_pos(n, "n")
  if (any(n < 3)) stop("`n` must be >= 3")
  if (any(beta_sd == 0)) stop("`beta_sd` must be nonzero")
  g <- .ps_grid(list(maf = maf, beta_sd = beta_sd, n = n, alpha = alpha))
  het <- 2 * g$maf * (1 - g$maf)
  g[, r2 := het * beta_sd^2]
  if (any(g$r2 >= 1)) stop("implied variance explained >= 1; reduce `beta_sd`")
  se <- sqrt((1 - g$r2) / (g$n * het))
  out <- .ps_finish(g, "qtl_additive", "per-allele effect (trait SD)", g$beta_sd, se, g$n - 2, g$alpha)
  out[, method := "additive linear regression, shifted t"]
  out[]
}
