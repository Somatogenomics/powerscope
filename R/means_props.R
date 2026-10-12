#' Two-group or paired comparison of means (t-test)
#'
#' All three quantities are exact for the t-test on normal data. Power and
#' Type S use the noncentral t distribution (both rejection tails, as
#' `stats::power.t.test(strict = TRUE)`). The exaggeration ratio
#' E(|mean difference| | significant) / |delta| is computed by integrating
#' over the sampling distribution of the pooled SD, which is independent of
#' the mean difference under normality. The `sig_cutoff` and
#' `min_sig_ratio` columns refer to a sample SD equal to `sd`.
#'
#' @param delta True difference in means (group 2 minus group 1), or mean
#'   within-pair difference when `paired = TRUE`.
#' @param sd Within-group SD (two-sample) or SD of the paired differences.
#'   Must be estimated from pilot data or literature.
#' @param n1 Group 1 size, or number of pairs when `paired = TRUE`.
#' @param n2 Group 2 size; defaults to `n1`. Ignored when `paired = TRUE`.
#' @param alpha Two-sided significance level.
#' @param paired Logical; paired design.
#' @return A `ps_result` data.table; all vector arguments are fully crossed.
#' @examples
#' ps_means(delta = 0.5, sd = 1, n1 = c(10, 20, 64))
#' @export
ps_means <- function(delta, sd, n1, n2 = NULL, alpha = 0.05, paired = FALSE) {
  if (any(delta == 0)) stop("`delta` must be nonzero")
  .check_pos(sd, "sd"); .check_pos(n1, "n1")
  g <- .ps_grid(list(delta = delta, sd = sd, n1 = n1, n2 = if (!paired) n2, alpha = alpha))
  if (paired) {
    g[, n2 := NA_real_]
    se <- g$sd / sqrt(g$n1); df <- g$n1 - 1
  } else {
    if (is.null(n2)) g[, n2 := n1]
    se <- g$sd * sqrt(1 / g$n1 + 1 / g$n2); df <- g$n1 + g$n2 - 2
  }
  if (any(df < 1)) stop("too few observations: degrees of freedom < 1")
  out <- .ps_finish(g, if (paired) "means_paired" else "means_two_group",
                    "difference in means", g$delta, se, df, g$alpha)
  q <- stats::qt(1 - g$alpha / 2, df); ncp <- g$delta / se
  p_hi <- stats::pt(q, df, ncp, lower.tail = FALSE)
  p_lo <- stats::pt(-q, df, ncp)
  ex <- mapply(.t_exaggeration, g$delta, se, df, q)
  out[, power := p_hi + p_lo]
  out[, type_s := ifelse(delta > 0, p_lo, p_hi) / power]
  out[, exaggeration := ex]
  out[, method := "exact t-test (noncentral t; SD integrated out)"]
  out[]
}

# E(|d| | significant)/|D| for the t-test: d ~ N(D, tau^2) independent of
# w = S^2/sigma^2 ~ chi^2_df / df; reject when |d| >= q * tau * sqrt(w).
.t_exaggeration <- function(D, tau, df, q) {
  inner <- function(w) {
    cut <- q * tau * sqrt(w)
    a <- (cut - D) / tau; b <- (-cut - D) / tau
    p <- stats::pnorm(a, lower.tail = FALSE) + stats::pnorm(b)
    e <- D * stats::pnorm(a, lower.tail = FALSE) + tau * stats::dnorm(a) -
         D * stats::pnorm(b) + tau * stats::dnorm(b)
    dw <- stats::dchisq(w * df, df) * df
    cbind(p * dw, e * dw)
  }
  P <- stats::integrate(function(w) inner(w)[, 1], 0, Inf, rel.tol = 1e-9)$value
  E <- stats::integrate(function(w) inner(w)[, 2], 0, Inf, rel.tol = 1e-9)$value
  E / P / abs(D)
}

# Normal estimate N(D, s^2) with an arbitrary rejection cutoff `cut`.
.normal_cut <- function(D, s, cut) {
  a <- (cut - D) / s; b <- (-cut - D) / s
  p_hi <- stats::pnorm(a, lower.tail = FALSE); p_lo <- stats::pnorm(b)
  power <- p_hi + p_lo
  e_hi <- D * p_hi + s * stats::dnorm(a)
  e_lo <- -D * p_lo + s * stats::dnorm(b)
  list(power = power, type_s = ifelse(D > 0, p_lo, p_hi) / power,
       exaggeration = (e_hi + e_lo) / power / abs(D))
}

# Exact enumeration over binomial outcomes for the two-proportion z-test.
.props_exact <- function(p1, p2, n1, n2, alpha, method) {
  x1 <- 0:n1; x2 <- 0:n2
  pr <- outer(stats::dbinom(x1, n1, p1), stats::dbinom(x2, n2, p2))
  est <- outer(x1 / n1, x2 / n2, function(a, b) b - a)
  if (method == "pooled") {
    pb <- outer(x1, x2, "+") / (n1 + n2)
    se <- sqrt(pb * (1 - pb) * (1 / n1 + 1 / n2))
  } else {
    se <- sqrt(outer(x1 / n1 * (1 - x1 / n1) / n1, x2 / n2 * (1 - x2 / n2) / n2, "+"))
  }
  sig <- se > 0 & abs(est) >= stats::qnorm(1 - alpha / 2) * se
  P <- sum(pr[sig]); D <- p2 - p1
  list(power = P,
       type_s = sum(pr[sig & sign(est) != sign(D)]) / P,
       exaggeration = sum(pr[sig] * abs(est[sig])) / P / abs(D))
}

#' Two-group comparison of proportions
#'
#' `method = "pooled"` is the usual two-sample z-test with the pooled
#' proportion in the null SE; `method = "wald"` uses the unpooled SE.
#' When `exact = TRUE` (default when `n1 * n2 <= 4e6`), power, Type S and the
#' exaggeration ratio are computed exactly by enumerating all binomial
#' outcomes. Otherwise a normal approximation is used; for `"pooled"` with
#' equal groups its power equals `stats::power.prop.test(strict = TRUE)`.
#' The normal approximation can be off by a few percentage points with
#' unequal groups or proportions near 0 or 1.
#'
#' @param p1,p2 True proportions in groups 1 and 2.
#' @param n1 Group 1 size.
#' @param n2 Group 2 size; defaults to `n1`.
#' @param alpha Two-sided significance level.
#' @param method `"pooled"` or `"wald"`.
#' @param exact `NULL` (automatic), `TRUE` or `FALSE`.
#' @return A `ps_result` data.table.
#' @examples
#' ps_props(p1 = 0.5, p2 = 0.75, n1 = c(30, 60, 90))
#' @export
ps_props <- function(p1, p2, n1, n2 = NULL, alpha = 0.05, method = c("pooled", "wald"),
                     exact = NULL) {
  method <- match.arg(method)
  .check_prob(p1, "p1"); .check_prob(p2, "p2"); .check_pos(n1, "n1")
  g <- .ps_grid(list(p1 = p1, p2 = p2, n1 = n1, n2 = n2, alpha = alpha))
  if (is.null(n2)) g[, n2 := n1]
  if (any(g$p1 == g$p2)) stop("`p1` and `p2` must differ")
  if (any(g$n1 != round(g$n1) | g$n2 != round(g$n2))) stop("`n1`, `n2` must be integers")
  D <- g$p2 - g$p1
  se1 <- sqrt(g$p1 * (1 - g$p1) / g$n1 + g$p2 * (1 - g$p2) / g$n2)
  out <- .ps_finish(g, paste0("props_", method), "difference in proportions",
                    D, se1, Inf, g$alpha)
  if (method == "pooled") {
    pbar <- (g$n1 * g$p1 + g$n2 * g$p2) / (g$n1 + g$n2)
    se0 <- sqrt(pbar * (1 - pbar) * (1 / g$n1 + 1 / g$n2))
    cut <- stats::qnorm(1 - g$alpha / 2) * se0
    r <- .normal_cut(D, se1, cut)
    out[, `:=`(sig_cutoff = cut, power = r$power, type_s = r$type_s,
               exaggeration = r$exaggeration, min_sig_ratio = cut / abs(D))]
  }
  use_exact <- if (is.null(exact)) g$n1 * g$n2 <= 4e6 else rep(exact, nrow(g))
  out[, method := paste0(method, ifelse(use_exact, ", exact enumeration", ", normal approximation"))]
  for (i in which(use_exact)) {
    r <- .props_exact(g$p1[i], g$p2[i], g$n1[i], g$n2[i], g$alpha[i], method)
    data.table::set(out, i, c("power", "type_s", "exaggeration"),
                    list(r$power, r$type_s, r$exaggeration))
  }
  out[]
}
