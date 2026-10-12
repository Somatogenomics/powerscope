#' Smallest sample size reaching a target power
#'
#' Integer bisection over a sample-size argument of any powerscope design
#' function. Power is assumed non-decreasing in that argument; the result is
#' checked at `n` and `n - 1`.
#'
#' @param design A design function, e.g. [ps_means()].
#' @param n_arg Name of the sample-size argument of `design` (e.g. `"n1"`,
#'   `"n_per_group"`, `"n_cases"`).
#' @param target_power Target power.
#' @param ... Other (scalar) arguments passed to `design`.
#' @param n_min,n_max Search range.
#' @param link Optional named list of other arguments that should equal
#'   `n_arg` times a ratio, e.g. `list(n2 = 1)` for equal group sizes.
#' @return A one-row result for the smallest `n` with power >= `target_power`,
#'   with attribute `"solved"` naming the argument; `NULL` with a warning if
#'   `n_max` is not enough.
#' @examples
#' ps_sample_size(ps_means, "n1", 0.8, delta = 0.5, sd = 1, link = list(n2 = 1))
#' @export
ps_sample_size <- function(design, n_arg, target_power = 0.8, ..., n_min = 2L,
                           n_max = 1e6L, link = NULL) {
  dots <- list(...)
  pw <- function(n) {
    a <- dots; a[[n_arg]] <- n
    for (nm in names(link)) a[[nm]] <- ceiling(n * link[[nm]])
    r <- do.call(design, a)
    if (nrow(r) != 1L) stop("design arguments must be scalars when solving")
    r
  }
  if (pw(n_max)$power < target_power) {
    warning("target power not reached at ", n_arg, " = ", n_max)
    return(NULL)
  }
  lo <- as.integer(n_min); hi <- as.integer(n_max)
  if (pw(lo)$power >= target_power) {
    res <- pw(lo)
  } else {
    while (hi - lo > 1L) {
      mid <- (lo + hi) %/% 2L
      if (pw(mid)$power >= target_power) hi <- mid else lo <- mid
    }
    res <- pw(hi)
  }
  data.table::setattr(res, "solved", n_arg)
  res
}

#' Minimum detectable effect at a target power
#'
#' Finds the effect size for which a design reaches `target_power`, by root
#' finding on a continuous effect argument.
#'
#' @inheritParams ps_sample_size
#' @param effect_arg Name of the effect argument of `design` (e.g. `"delta"`,
#'   `"log2fc"`, `"odds_ratio"`).
#' @param interval Search interval for the effect, on one side of the null.
#' @return One-row result at the minimum detectable effect, with attribute
#'   `"solved"`.
#' @examples
#' ps_mde(ps_means, "delta", 0.8, interval = c(1e-4, 10), sd = 1, n1 = 20, n2 = 20)
#' @export
ps_mde <- function(design, effect_arg, target_power = 0.8, ..., interval) {
  dots <- list(...)
  pw <- function(e) { a <- dots; a[[effect_arg]] <- e; do.call(design, a) }
  f <- function(e) pw(e)$power - target_power
  fl <- f(interval[1]); fu <- f(interval[2])
  if (sign(fl) == sign(fu)) {
    warning("target power not bracketed by `interval`")
    return(NULL)
  }
  root <- stats::uniroot(f, interval, tol = 1e-10)$root
  res <- pw(root)
  data.table::setattr(res, "solved", effect_arg)
  res
}

#' Per-test significance level that controls the false discovery rate
#'
#' For `m` tests of which a fraction `pi1` carry a true effect, the
#' Benjamini-Hochberg procedure at level `q` behaves approximately like a
#' fixed per-test level `alpha*` solving
#' `alpha* = q * m1 * power(alpha*) / (m0 * (1 - q))` (Jung 2005).
#' The per-test level is then used in place of `alpha` in any design.
#'
#' @param power_fun Function of `alpha` returning the average power across
#'   the truly non-null tests.
#' @param q Target FDR.
#' @param m Number of tests.
#' @param pi1 Fraction of tests with a true effect. This must be estimated
#'   from your data or prior studies; it has no default.
#' @return The per-test level `alpha*`, with attributes `expected_true` and
#'   `expected_false` (expected discoveries). `NA` if no positive solution
#'   exists (power too low for any expected discoveries).
#' @references Jung SH (2005). Sample size for FDR-control in microarray data
#'   analysis. *Bioinformatics* 21:3097-3104. \doi{10.1093/bioinformatics/bti456}
#' @examples
#' pf <- function(a) ps_means(delta = 1.5, sd = 1, n1 = 20, n2 = 20, alpha = a)$power
#' ps_fdr_alpha(pf, q = 0.05, m = 10000, pi1 = 0.05)
#' @export
ps_fdr_alpha <- function(power_fun, q, m, pi1) {
  .check_prob(q, "q"); .check_prob(pi1, "pi1"); .check_pos(m, "m")
  m1 <- m * pi1; m0 <- m - m1
  k <- q * m1 / (m0 * (1 - q))
  f <- function(la) { a <- exp(la); a - k * power_fun(a) }
  grid <- seq(log(1e-14), log(q), length.out = 200)
  fv <- vapply(grid, f, numeric(1))
  neg <- which(fv < 0)
  if (!length(neg)) {
    warning("no positive per-test alpha satisfies the FDR equation; power is too low")
    return(NA_real_)
  }
  # largest root: last sign change from negative to positive
  i <- max(neg)
  if (i == length(grid)) {
    a <- q  # FDR satisfied even at alpha = q
  } else {
    a <- exp(stats::uniroot(f, c(grid[i], grid[i + 1]), tol = 1e-12)$root)
  }
  pw <- power_fun(a)
  structure(a, expected_true = m1 * pw, expected_false = m0 * a, avg_power = pw)
}
